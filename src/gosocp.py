from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def gosocp(log,all_data):


    log.joint(' creating ampl object ...\n')
    ampl = AMPL()

    all_data['ampl_object'] = ampl
        
    modfile = all_data['modfile']
    solver  = all_data['solver']

    log.joint(" reading modfile ...\n")
    
    t0 = time.time()
    ampl.read('../modfiles/' + modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    
    ampl.eval("option display_precision 0;")
    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)    

    if True:
        ampl.setOption('presolve',0)
        log.joint(' AMPL presolve off\n')

    log.joint(" solver set to " + solver + "\n")
    ampl.eval("option show_stats 2;")

    if all_data['solver'] == 'gurobi_ampl':
        ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 TimeLimit=1000';")
        #ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 iisfind=1 writemodel=jabr.lp resultfile=jabr.ilp';")

    if all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=5 maxtime_real=1000 convex=1';") #bar_conic_enable=1  honorbnds=1
        elif all_data['multistart']:
            ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';")
        elif all_data['knitropresolveoff']:
            ampl.eval("option knitro_options 'presolve=0';")       
        else:
            ampl.eval("option knitro_options 'blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';") 

    if all_data['solver'] == 'mosek':
        ampl.eval("option mosek_options 'outlev=1 sol:chk:feastol=1e-08 chk:mode=2';")
        #cvt:socp=0     
    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']


    IDtoCountmap = all_data['IDtoCountmap']
        
    #buses
    buses        = {}
    Pd           = {}
    Qd           = {}
    Vmax         = {}
    Vmin         = {}    
    branches_f   = {}
    branches_t   = {}
    bus_gens     = {}
    bus_Gs       = {}
    bus_Bs       = {}
    Vinit        = {}


    for bus in all_data['buses'].values():
        buscount = bus.count

        if ( len(bus.frombranchids) == 0 ) and ( len(bus.tobranchids) == 0 ): #else we could have just skipped bus.Gs and bus.Bs
            log.joint(' flawed bus ' + str(bus.nodeID) + ' busid ' + str(buscount) + '\n') #we are not skipping it
        #    breakexit('check')
        #    continue

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
    
    for bus in all_data['buses'].values(): #if REMOVING flawed buses, change this
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
        bus_gens_set[bus.count].setValues(bus_gens[bus.count])
    
    #branches 
    branches  = {}
    U         = {}
    Gtt       = {}
    Btt       = {}
    Gff       = {}
    Bff       = {}
    Gtf       = {}
    Btf       = {}
    Gft       = {}
    Bft       = {}
    bus_f     = {}
    bus_t     = {}
    CSmax     = {}
    c_ubound  = {}
    c_lbound  = {}
    s_ubound  = {}
    s_lbound  = {}
    Cinit     = {}
    Sinit     = {}
    
    for branch in all_data['branches'].values():
        branchcount = branch.count
        if branch.status != 1: #the reader does not consider branches whose status is 0, hence this should never be true
            log.joint(' branch ' + str(branchcount) + ' off, we skip it\n')
            breakexit('check')
            continue
        
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

        Cinit[branchcount] = 1
        Sinit[branchcount] = 0

        # cs bounds
        # Assumption: zero angle difference is always allowed
        maxprod     = Vmax[branch.id_f] * Vmax[branch.id_t]
        minprod     = Vmin[branch.id_f] * Vmin[branch.id_t]

        ubound      = maxprod
        lbound      = -maxprod
        maxanglerad = branch.maxangle_rad
        minanglerad = branch.minangle_rad
        
        # Cosine
        if maxanglerad <= 0.5*math.pi:
            # In this case minangle <= 0                                              
            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)  # Which is negative          
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod
                
        elif maxanglerad <= math.pi:
            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*math.cos(maxanglerad)
            elif minanglerad >= -math.pi:
                lbound = maxprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5*math.pi:
            lbound = -maxprod    

        elif maxanglerad <= 2*math.pi:
            lbound = -maxprod

        else:
            ubound = maxprod
            lbound = -maxprod

        c_ubound[branchcount] = ubound
        c_lbound[branchcount] = lbound
        
        # Sine                                                                         
        if maxanglerad <= 0.5*math.pi:
            ubound = maxprod*math.sin(maxanglerad)

            if  minanglerad >= -0.5*math.pi:
                lbound = maxprod*math.sin(minanglerad)
            elif  minanglerad >= -math.pi:
                lbound = -maxprod
            elif  minanglerad >= -1.5*math.pi:
                ubound = maxprod*max( math.sin(maxanglerad), math.sin(minanglerad))
                lbound = -maxprod
            else:
                ubound = maxprod
                lbound = -maxprod

        elif maxanglerad <= math.pi:
            ubound = maxprod

            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*math.sin(minanglerad)
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5*math.pi:
            ubound = maxprod

            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*min(math.sin(maxanglerad), math.sin(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod
        else:
            ubound = maxprod
            lbound = -maxprod

        s_ubound[branchcount] = ubound
        s_lbound[branchcount] = lbound
    
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
    ampl.get_parameter('Cinit').setValues(Cinit)
    ampl.get_parameter('Sinit').setValues(Sinit)
    ampl.get_parameter('c_ubound').setValues(c_ubound)
    ampl.get_parameter('c_lbound').setValues(c_lbound)
    ampl.get_parameter('s_ubound').setValues(s_ubound)
    ampl.get_parameter('s_lbound').setValues(s_lbound)
    
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

    log.joint(" sets and parameters loaded\n")

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
    Qf = ampl.get_variable("Qf")
    Qt = ampl.get_variable("Qt")
    Pg = ampl.get_variable("Pg")
    Qg = ampl.get_variable("Qg")
    dic_v = v.get_values().to_dict()
    dic_c = c.get_values().to_dict()
    dic_s = s.get_values().to_dict()
    dic_Pg = Pg.get_values().to_dict()
    dic_Qg = Qg.get_values().to_dict()
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

    log.joint(" writing casename, modfile, obj and runtime to summary_socp.log\n")

    summary_socp = open("summary_socp.log","a+")

    summary_socp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile'] +  ' obj ' + str(objval) + ' runtime ' + str(timesofar) + '\n')

    summary_socp.close()
    
    log.joint(" ------------------------\n")
    
    #breakexit("write down solution?")
    
    branches_data = all_data['branches']
    buses_data    = all_data['buses']
    IDtoCountmap  = all_data['IDtoCountmap']
    tolerance     = 1e-05
    
    casefilename  = all_data['casefilename']
    filename      = 'ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.txt'
    filenamelp    = 'ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.lp'
    #filename      = 'sols_socp/ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.txt'
    #filenamelp    = 'sols_socp/ksol_' + all_data['modfile'] + '_' + all_data['casename'] + '.lp'
    
    thefile       = open(filename,'w+')
    thefilelp     = open(filenamelp,'w+')
    
    log.joint(' writing solution to ' + filename + ' and to ' + filenamelp + '\n')

    thefile.write('value ' + str(objval) + '\n')
    
    thefile.write('voltages:\n') # add power injections?

    for buscount in buses.keys():
        v2val      = dic_v[buscount]
        vval       = (v2val)**(0.5)
        f          = buses_data[buscount].nodeID
        line = 'bus ' + str(buscount) + ' M ' + str(vval) + '\n'
        thefile.write(line)
        
        linelp_v = str(v2val - tolerance) + ' <= c_' + str(f) + '_' + str(f) + ' <= ' + str(v2val + tolerance) + '\n'
        thefilelp.write(linelp_v)
    
    thefile.write('power flows and cs variables:\n')
    
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

    thefile.write('generation:\n')
        
    for genid in gens.keys():
        gen_nodeID = gens[genid]
        line_gen = 'genid ' + str(genid) + ' bus ' + str(gen_nodeID) + ' GP ' + str(dic_Pg[genid]) + ' GQ ' + str(dic_Qg[genid]) + '\n'
        thefile.write(line_gen)
        
    thefile.close()
    thefilelp.close()

    log.joint(' done writing knitro solution to file\n')


def getsol_knitro(log,all_data):

    
    casename = all_data['casename']
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

    
    casename = all_data['casename']
    filename = 'solution_va_'+ casename +'.txt'
    #filename = 'mp_sols/solution_va_'+ casename +'.txt'
    
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
