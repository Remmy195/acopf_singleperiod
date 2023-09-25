from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math



def goqp(log,all_data):

    log.joint(' creating ampl object ...\n')
    ampl = AMPL()

    all_data['ampl_object'] = ampl
        
    modfile = all_data['modfile']
    solver  = all_data['solver']

    log.joint(" reading modfile ...\n")
    
    t0 = time.time()
    ampl.read(modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    
    ampl.eval("option display_precision 0;")
    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)    
    ampl.setOption('presolve',0)
    log.joint(" solver set to " + solver + "\n")
    ampl.eval("option show_stats 2;")
    #ampl.eval("option presolve_eps=1e-5;")

    if all_data['solver'] == 'gurobi_ampl':
        ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 iisfind=1 writemodel=jabr.lp resultfile=jabr.ilp';")

    if all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'convex=1 feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=7';")
        else:
            ampl.eval("option knitro_options 'convex=1 blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';") #try linsolver_numthreads=1 and choose linear solver... #number of threads?';")

    #bar_conic_enable=1 --> too slow for some problems ...
    # bar_refinement=1
    
    #if all_data['solver'] == 'knitroampl' and all_data['knitropresolveoff']:
    #    ampl.eval("option knitro_options 'presolve=0';")
        
    #if all_data['solver'] == 'knitroampl' and all_data['conic']:
    #    ampl.eval("option knitro_options 'bar_conic_enable=1';")

    #if all_data['solver'] == 'knitroampl' and all_data['conic'] and all_data['knitropresolveoff']:
    #    ampl.eval("option knitro_options 'bar_conic_enable=1 presolve=0';")        
    #if all_data['solver'] == 'knitroampl' and all_data['multistart']:
    #    ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=3 ms_terminate =1';")


    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']
	#mp_angle = all_data['mp_angle']
        print(mp_vm)

    IDtoCountmap = all_data['IDtoCountmap']
        
    #buses
    buses = {}
    Pd = {}
    Qd = {}
    Vmax = {}
    Vmin = {}    
    branches_f = {}
    branches_t = {}
    bus_gens = {}
    bus_Gs = {}
    bus_Bs = {}
    Vinit  = {}

    for bus in all_data['buses'].values():
        buscount = bus.count
        buses[buscount] = bus.nodeID 
        Pd[buscount] = bus.Pd
        Qd[buscount] = bus.Qd
        if all_data['fix_point']:
            log.joint(' buscount ' + str(buscount) + ' VM ' + str(mp_vm[buscount]) + '\n')
            Vmin[buscount] = mp_vm[buscount] - tolerance                                                                                          
            Vmax[buscount] = mp_vm[buscount] + tolerance
        else:
            Vmax[buscount] = bus.Vmax
            Vmin[buscount] = bus.Vmin

        Vinit[buscount]      = 1

        branches_f[buscount] = []
        branches_t[buscount] = []

        bus_gens[buscount] = bus.genidsbycount #load this as a sparse ds

        if ( bus.Gs != 0 ) and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
            bus_Gs[buscount] = bus.Gs
        if ( bus.Bs != 0 ) and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
            bus_Bs[buscount] = bus.Bs
        
        for branchid in bus.frombranchids.values():
            branches_f[buscount].append(branchid) 
            
        for branchid in bus.tobranchids.values():
            branches_t[buscount].append(branchid)      

    
    ampl.getSet('buses').setValues(list(buses))
    ampl.getSet('bus_Gs').setValues(list(bus_Gs))
    ampl.getSet('bus_Bs').setValues(list(bus_Bs))

    ampl.get_parameter('Gs').setValues(list(bus_Gs.values()))
    ampl.get_parameter('Bs').setValues(list(bus_Bs.values()))
    ampl.get_parameter('Pd').set_values(Pd)
    ampl.get_parameter('Qd').set_values(Qd)
    ampl.get_parameter('Vmax').set_values(Vmax)
    ampl.get_parameter('Vmin').set_values(Vmin)
    ampl.get_parameter('Vinit').set_values(Vinit)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set   = ampl.getSet('bus_gens')
    
    for bus in all_data['buses'].values():
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
        bus_gens_set[bus.count].setValues(bus_gens[bus.count])
    
    #branches
    branches = {}
    U = {}
    Gtt = {}
    Btt = {}
    Gff = {}
    Bff = {}
    Gtf = {}
    Btf = {}
    Gft = {}
    Bft = {}
    bus_f = {}
    bus_t = {}
    CSmax = {}

    i2max       = {}
    g           = {}
    b           = {}
    bshunt      = {}
    ratio       = {}
    phase_angle = {}
    
    counter = 0 
    for branch in all_data['branches'].values():
        branchcount = branch.count
        branches[branchcount] = (branch.id_f,branch.id_t)
        U[branchcount] = branch.limit
        Gtt[branchcount] = branch.Gtt
        Btt[branchcount] = branch.Btt
        Gff[branchcount] = branch.Gff
        Bff[branchcount] = branch.Bff
        Gtf[branchcount] = branch.Gtf
        Btf[branchcount] = branch.Btf
        Gft[branchcount] = branch.Gft
        Bft[branchcount] = branch.Bft
        CSmax[branchcount] = Vmax[branch.id_f] * Vmax[branch.id_t]
        bus_f[branchcount] = branch.id_f
        bus_t[branchcount] = branch.id_t

        y                        = branch.y
        i2max[branchcount]       = branch.limit**2 / Vmin[branch.id_f] * Vmin[branch.id_f]
        g[branchcount]           = y.real
        b[branchcount]           = y.imag
        bshunt[branchcount]      = branch.bc
        ratio[branchcount]       = branch.ratio
        phase_angle[branchcount] = branch.angle_rad
        

    ampl.getSet('branches').setValues(list(branches))
    ampl.get_parameter('Gtt').setValues(Gtt)
    ampl.get_parameter('Btt').setValues(Btt)
    ampl.get_parameter('Gff').setValues(Gff)
    ampl.get_parameter('Bff').setValues(Bff)
    ampl.get_parameter('Gtf').setValues(Gtf)
    ampl.get_parameter('Btf').setValues(Btf)
    ampl.get_parameter('Gft').setValues(Gft)
    ampl.get_parameter('Bft').setValues(Bft) 
    ampl.get_parameter('U').setValues(U)
    ampl.get_parameter('bus_f').setValues(bus_f)
    ampl.get_parameter('bus_t').setValues(bus_t)
    ampl.get_parameter('CSmax').setValues(CSmax)

    ampl.get_parameter('i2max').setValues(i2max)    
    ampl.get_parameter('g').setValues(g) 
    ampl.get_parameter('b').setValues(b)
    ampl.get_parameter('bshunt').setValues(bshunt)    
    ampl.get_parameter('ratio').setValues(ratio)
    ampl.get_parameter('phase_angle').setValues(phase_angle)

    #gens
    gens      = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} #we assume at most quadratic cost (for the moment)

    for gen in all_data['gens'].values():
        gencount = gen.count #this corresponds to bus count
        gens[gencount] = gen.nodeID

        #if gen.status != 0 and gen.status != 1:
        #    breakexit('check')
        
        Pmax[gencount]      = gen.Pmax * gen.status
        Pmin[gencount]      = gen.Pmin * gen.status
        Qmax[gencount]      = gen.Qmax * gen.status
        Qmin[gencount]      = gen.Qmin * gen.status

        buscount = IDtoCountmap[gen.nodeID]
        bus      = all_data['buses'][buscount]

        if bus.nodetype == 3:
            Qmax[gencount] = 2 * all_data['summaxgenQ']
            Qmin[gencount] = - 2 * all_data['summaxgenQ']            

        #print(" cost vector",gen.costvector) #check small case5.m
        #print(" cost degree",gen.costdegree)

        fixedcost[gencount] = gen.costvector[2] * gen.status
        lincost[gencount] = gen.costvector[1]
        quadcost[gencount] = gen.costvector[0]
    
    ampl.getSet('gens').setValues(list(gens))
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)


    #cuts
    add_cuts(log,all_data)
    
    ampl.get_parameter('num_jabr').set(all_data['num_jabr'])
    ampl.getSet('jabr_cuts').setValues(list(all_data['jabr_cuts'].keys()))
    ampl.get_parameter('c_jabr').setValues(all_data['c_jabr'])
    ampl.get_parameter('s_jabr').setValues(all_data['s_jabr'])
    ampl.get_parameter('vf_jabr').setValues(all_data['vf_jabr'])
    ampl.get_parameter('vt_jabr').setValues(all_data['vt_jabr'])

    ampl.get_parameter('num_i2').set(all_data['num_i2'])
    ampl.getSet('i2_cuts').setValues(list(all_data['i2_cuts'].keys()))
    ampl.get_parameter('Pf_i2').setValues(all_data['Pf_i2'])
    ampl.get_parameter('Qf_i2').setValues(all_data['Qf_i2'])
    ampl.get_parameter('vf_i2').setValues(all_data['vf_i2'])
    ampl.get_parameter('i2f_i2').setValues(all_data['i2f_i2'])

    ampl.get_parameter('num_limf').set(all_data['num_limf'])
    ampl.getSet('limf_cuts').setValues(list(all_data['limf_cuts'].keys()))
    ampl.get_parameter('Pf_lim').setValues(all_data['Pf_lim'])
    ampl.get_parameter('Qf_lim').setValues(all_data['Qf_lim'])

    ampl.get_parameter('num_limt').set(all_data['num_limt'])
    ampl.getSet('limt_cuts').setValues(list(all_data['limt_cuts'].keys()))
    ampl.get_parameter('Pt_lim').setValues(all_data['Pt_lim'])
    ampl.get_parameter('Qt_lim').setValues(all_data['Qt_lim'])
    

    log.joint(" sets, parameters and cuts loaded\n")

    log.joint(" saving processed data to all_data\n")
    
    expand = True
    if expand:
        time1 = time.time()
        filename = 'basemodel.out'
        doNLP = 1
        if doNLP:
            filename = 'NLP.out'
        log.joint('Now expanding to %s.\n'%(filename))
        amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb,_var[j].ub);' #shows full mod    
        modelout = ampl.getOutput(amplstate)
        outfile = open(filename,"w")
        outfile.write("model = " + str(modelout) + "\n")
        outfile.close()
        
    #SOLVE
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()

    log.joint(" ------------------------\n")


    #GET SOLUTION
    total_cost = ampl.get_objective("total_cost")
    v = ampl.get_variable("v")
    c = ampl.get_variable("c")
    s = ampl.get_variable("s")
    Pf = ampl.get_variable("Pf")
    Pt = ampl.get_variable("Pt")
    Pg = ampl.get_variable("Pg")
    Qf = ampl.get_variable("Qf")
    Qt = ampl.get_variable("Qt")
    if all_data['modfile'] == 'jabr_i2.mod':
        i2f    = ampl.get_variable("i2f")
        dic_i2f = i2f.get_values().to_dict()
    dic_v = v.get_values().to_dict()
    dic_c = c.get_values().to_dict()
    dic_s = s.get_values().to_dict()
    dic_Pg = Pg.get_values().to_dict()
    dic_Pf = Pf.get_values().to_dict()
    dic_Pt = Pt.get_values().to_dict()
    dic_Qf = Qf.get_values().to_dict()
    dic_Qt = Qt.get_values().to_dict()
    
    dicPLoss = {}
    dicQLoss = {}
    
    for branch in dic_Pf.keys():
        dicPLoss[branch] = dic_Pf[branch] + dic_Pt[branch] 
        dicQLoss[branch] = dic_Qf[branch] + dic_Qt[branch] 

    objval    = total_cost.get().value()
    timesofar = time.time() - all_data['T0']
    
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" objective " + str(objval) + '\n')
    log.joint(" total active power generation " + str(sum(dic_Pg.values())) + '\n')
    log.joint(" total active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" total active power loss " + str(sum(dicPLoss.values())) + '\n')
    log.joint(" total net reactive power loss " + str(sum(dicQLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')

    log.joint(" writing casename, modfile, obj and runtime to summary_qp.log\n")

    summary_qp = open("summary_qp.log","a+")

    summary_qp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile'] +  ' obj ' + str(objval) + ' runtime ' + str(timesofar) + '\n')

    summary_qp.close()
    
    log.joint(" ------------------------\n")
    
    breakexit("write down solution?")
    
    
    branches_data = all_data['branches']
    buses_data    = all_data['buses']
    IDtoCountmap  = all_data['IDtoCountmap']
    tolerance     = 1e-05 #default 1e-5
    
    casefilename  = all_data['casefilename'] 
    filename      = 'sols_qp/ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.txt'
    filenamelp    = 'sols_qp/ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.lp'
    
    thefile       = open(filename,'w+')
    thefilelp     = open(filenamelp,'w+')
    
    log.joint(' writing solution to ' + filename + ' and to ' + filenamelp + '\n')

    thefile.write('value ' + str(total_cost.get().value()) + '\n')
    
    thefile.write('voltages:\n')

    for buscount in buses.keys():
        v2val      = dic_v[buscount]
        vval       = (v2val)**(0.5)
        f          = buses_data[buscount].nodeID
        line = 'bus ' + str(buscount) + ' M ' + str(vval) + '\n'
        thefile.write(line)
        
        linelp_v = str(v2val - tolerance) + ' <= c_' + str(f) + '_' + str(f) + ' <= ' + str(v2val + tolerance) + '\n'
        thefilelp.write(linelp_v)
    
    thefile.write('power flows and cs:\n')
    
    for branchid in branches.keys():

        branch     = branches_data[branchid]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        Pfval      = dic_Pf[branchid]
        Ptval      = dic_Pt[branchid]
        Qfval      = dic_Qf[branchid]
        Qtval      = dic_Qt[branchid]
        cftval     = dic_c[branchid]
        sftval     = dic_s[branchid]

        if all_data['modfile'] == 'jabr_i2.mod':
            i2fval = dic_i2f[branchid]
            line   = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' i2ft ' + str(i2fval) + '\n'
        else:
            line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + '\n'

        thefile.write(line)

        linelp_Pf = str(Pfval-tolerance) + ' <= P_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Pfval+tolerance) + '\n'
        thefilelp.write(linelp_Pf)
        linelp_Pt = str(Ptval-tolerance) + ' <= P_' + str(branchid) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Ptval+tolerance) + '\n'
        thefilelp.write(linelp_Pt)
        linelp_Qf = str(Qfval-tolerance) + ' <= Q_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Qfval+tolerance) + '\n'
        thefilelp.write(linelp_Qf)

        linelp_Qt = str(Qtval-tolerance) + ' <= Q_' + str(branchid) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Qtval+tolerance) + '\n'

        thefilelp.write(linelp_Qt)
        
        linelp_cft = str(cftval-tolerance) + ' <= c_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(cftval+tolerance) + '\n'

        thefilelp.write(linelp_cft)

        linelp_sft = str(sftval-tolerance) + ' <= s_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(sftval+tolerance) + '\n' 

        thefilelp.write(linelp_sft)

        if all_data['modfile'] == 'jabr_i2.mod' or all_data['modfile'] == 'qp.mod':
            linelp_i2f = str(i2fval-tolerance) + ' <= i2_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(i2fval+tolerance) + '\n' 

            thefilelp.write(linelp_i2f)
            
    thefile.close()
    thefilelp.close()

    log.joint(' done writing knitro solution to file\n')

def add_cuts(log,all_data):

    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif '_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    #elif '_line' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 5]
    elif '_line5' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = 'cuts/cuts_' + original_casename + '.txt' 
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    theround  = lines[0].split()[3]
    firstline = lines[1].split()
    jabr      = 1

    if firstline[0] == '#Jabr-envelope':
        numjabrcuts = int(firstline[3])
        log.joint(' number of Jabr-envelope cuts = ' + str(numjabrcuts) + ' at round ' + str(theround) + '\n')
    elif firstline[0] == '#i2-envelope':
        jabr      = 0
        i2        = 1
        numi2cuts = int(firstline[3])
        log.joint(' no Jabr-envelope cuts to add\n')
        log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + ' at round ' + str(theround) + '\n')
    elif firstline[0] == '#limit-envelope':
        jabr         = 0
        i2           = 0
        numlimitcuts = firstline[3]
        log.joint(' no Jabr nor i2-envelope cuts to add\n')
        log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) + ' at round ' + str(theround) + '\n')
    else:
        log.joint(' no cuts added\n')
        return None

    linenum = 2

    buses          = all_data['buses']
    branches       = all_data['branches']
    IDtoCountmap   = all_data['IDtoCountmap']

    jabr_cuts      = {}
    i2_cuts        = {}
    limf_cuts      = {}
    limt_cuts      = {}
    
    c_jabr         = {}
    s_jabr         = {}
    vf_jabr        = {}
    vt_jabr        = {}

    vf_i2          = {}
    i2f_i2         = {}
    Pf_i2          = {}
    Qf_i2          = {}

    Pf_lim         = {}
    Qf_lim         = {}

    Pt_lim         = {}
    Qt_lim         = {}

    num_jabr       = 0
    num_i2         = 0
    num_limf       = 0
    num_limt       = 0
    
    if jabr:
        log.joint(' adding Jabr-envelope cuts ...\n')
    elif i2:
        log.joint(' adding i2-envelope cuts ...\n')
    else:
        log.joint(' adding limit-envelope cuts ...\n')

    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            numi2cuts = int(thisline[3])
            log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
            linenum += 1
            jabr = 0
            i2   = 1
            continue

        elif thisline[0] == '#limit-envelope' and i2:
            numlimcuts = int(thisline[3])
            log.joint(' number of limit-envelope cuts = ' + str(numlimcuts) + '\n')
            linenum += 1
            i2       = 0
            continue

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            cutid     = int(thisline[7])
            rnd       = int(thisline[9])
            coeff_cft = float(thisline[15])
            coeff_sft = float(thisline[17])
            coeff_cff = float(thisline[19])
            coeff_ctt = float(thisline[21])


            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                log.joint(' branchid ' + str(branchid) + ' branch.count ' + str(branch.count) + ' f ' + str(f) + ' branch.f ' + str(branch.f) + ' t ' + str(t) + ' branch.t ' + str(branch.t) + '\n')
                breakexit('bug')

            log.joint(' --> new Jabr-envelope cut\n')
            log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + '\n' )
            log.joint(' LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' + str(coeff_sft) + ' cff ' + str(coeff_cff) + ' ctt ' + str(coeff_ctt) + '\n' )
            
            num_jabr += 1
            
            jabr_cuts[(num_jabr,branchid)] = (coeff_cft,coeff_sft,coeff_cff,coeff_ctt)
            c_jabr[num_jabr]   = coeff_cft
            s_jabr[num_jabr]   = coeff_sft
            vf_jabr[num_jabr]  = coeff_cff
            vt_jabr[num_jabr]  = coeff_ctt
            
            linenum += 1

        elif (jabr == 0) and i2:
            
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            coeff_Pft  = float(thisline[15])
            coeff_Qft  = float(thisline[17])
            coeff_cff  = float(thisline[19])
            coeff_i2ft = float(thisline[21])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')

            log.joint(' --> new i2-envelope cut\n')
            log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + '\n' )
            log.joint(' LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' + str(coeff_Qft) + ' cff ' + str(coeff_cff) + ' i2ft ' + str(coeff_i2ft) + '\n' )
            
            num_i2 += 1
            
            i2_cuts[(num_i2,branchid)] = (coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft)
            Pf_i2[num_i2]   = coeff_Pft
            Qf_i2[num_i2]   = coeff_Qft
            vf_i2[num_i2]   = coeff_cff
            i2f_i2[num_i2]  = coeff_i2ft
            
            linenum += 1

        elif (jabr == 0) and (i2 == 0):
                
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')
            
            if thisline[14] == 'Pft':
                coeff_P    = float(thisline[15])
                coeff_Q    = float(thisline[17])

                num_limf += 1
                
                log.joint(' --> new limit-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' Pft ' + str(coeff_P) + ' Qft ' + str(coeff_Q) + '\n')
                
                limf_cuts[(num_limf,branchid)] = (coeff_P,coeff_Q)
                Pf_lim[num_limf]   = coeff_P
                Qf_lim[num_limf]   = coeff_Q

                linenum += 1
            elif thisline[14] == 'Ptf':
                coeff_P    = float(thisline[15])
                coeff_Q    = float(thisline[17])

                num_limt += 1

                log.joint(' --> new limit-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' t ' + str(t) + ' f ' + str(f) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' Ptf ' + str(coeff_P) + ' Qtf ' + str(coeff_Q) + '\n')

                limt_cuts[(num_limt,branchid)] = (coeff_P,coeff_Q)
                Pt_lim[num_limt]   = coeff_P
                Qt_lim[num_limt]   = coeff_Q

                linenum += 1                
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')


    if num_jabr != numjabrcuts:
        log.joint(' num jabr-cuts .txt file ' + str(numjabrcuts) + ' num jabr-cuts added ' + str(num_jabr) + '\n')
        log.joint(' potential bug\n')

    if num_i2 != numi2cuts:
        log.joint(' num i2-cuts .txt file ' + str(numi2cuts) + ' num i2-cuts added ' + str(num_i2) + '\n')
        log.joint(' potential bug\n')
                
    if (num_limf + num_limt) != numlimcuts:
        log.joint(' num lim-cuts .txt file ' + str(numlimcuts) + ' num i2-cuts added ' + str((num_limf + num_limt)) + '\n')
        log.joint(' potential bug\n')
            
    all_data['jabr_cuts']      = jabr_cuts
    all_data['i2_cuts']        = i2_cuts
    all_data['limf_cuts']      = limf_cuts
    all_data['limt_cuts']      = limt_cuts
    
    all_data['c_jabr']         = c_jabr
    all_data['s_jabr']         = s_jabr
    all_data['vf_jabr']        = vf_jabr
    all_data['vt_jabr']        = vt_jabr

    all_data['Pf_i2']          = Pf_i2
    all_data['Qf_i2']          = Qf_i2
    all_data['vf_i2']          = vf_i2
    all_data['i2f_i2']         = i2f_i2

    all_data['Pf_lim']         = Pf_lim
    all_data['Qf_lim']         = Qf_lim

    all_data['Pt_lim']         = Pt_lim
    all_data['Qt_lim']         = Qt_lim

    all_data['num_jabr']       = num_jabr
    all_data['num_i2']         = num_i2
    all_data['num_limf']       = num_limf
    all_data['num_limt']       = num_limt 
            


def getsol_knitro(log,all_data):

    casefilename = all_data['casefilename']
    casename = casefilename.split('/')[2].split('.')[0]
    #filename = 'knitro_sols/'+ casename +'.txt'
    filename = 'ksol_'+ casename +'.txt'
    try:
        thefile = open(filename, "r")
        lines = thefile.readlines()
        lenlines = len(lines)
        thefile.close()
    except:
        log.stateandquit("cannot open file " + datafilename)
        sys.exit("failure")

    branches     = all_data['branches']
    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    mp_vm      = {}
    mp_angle   = {}
    mp_cvalues = {}
    mp_svalues = {}

    mp_Pfvalues = {}
    mp_Ptvalues = {}
    mp_Qfvalues = {}
    mp_Qtvalues = {}
    linenum = 2

    log.joint(' reading file\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()

        if thisline[0] == 'bus':
            buscount        = int(thisline[1])
            mp_vm[buscount] = float(thisline[3])
            #mp_cvalues[bus] = mp_vm[bus]**2
        elif thisline[0] == 'branch':
            branchcount     = int(thisline[1])
            mp_Pfvalues[branchcount] = float(thisline[7])
            mp_Ptvalues[branchcount] = float(thisline[9])
            mp_Qfvalues[branchcount] = float(thisline[11])
            mp_Qtvalues[branchcount] = float(thisline[13])
            #mp_cvalues[branchcount] = float(thisline[15])
            #mp_svalues[branchcount] = float(thisline[17])
        linenum += 1

    all_data['mp_Pfvalues'] = mp_Pfvalues
    all_data['mp_Ptvalues'] = mp_Ptvalues
    all_data['mp_Qfvalues'] = mp_Qfvalues
    all_data['mp_Qtvalues'] = mp_Qtvalues

    all_data['mp_vm']      = mp_vm
    #all_data['mp_cvalues'] = mp_cvalues
    #all_data['mp_svalues'] = mp_svalues
    

def getsol(log,all_data):
    casefilename = all_data['casefilename']
    casename = casefilename.split('/')[2].split('.')[0]
    filename = 'mp_sols/solution_va_'+ casename +'.txt'
    
    log.joint(" reading file matpower solution volt magnitudes and angles " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    lenlines    = len(lines)

    mp_vm       = {}
    mp_angle    = {}
    
    buses                  = all_data['buses']
    IDtoCountmap           = all_data['IDtoCountmap']

    linenum  = 1

    while linenum < lenlines:
        thisline = lines[linenum].split(',')
        bus_id   = int(thisline[0])
        buscount = IDtoCountmap[bus_id]
        bus      = buses[buscount]

        mp_vm[buscount]        = float(thisline[1])
        mp_angle[buscount]     = float(thisline[2]) * math.pi / 180

        log.joint('bus ' + str(bus_id) + ' v ' + str(mp_vm[buscount]) + ' angle ' + str(mp_angle[buscount]) + '\n')

        linenum += 1
    
    all_data['mp_vm']    = mp_vm
    all_data['mp_angle'] = mp_angle    
