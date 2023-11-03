from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

        
def goac(log,all_data):


    log.joint(' creating ampl object ...\n')
    ampl = AMPL()

    all_data['ampl_object'] = ampl
        
    modfile = all_data['modfile']
    solver  = all_data['solver']

    log.joint(" reading modfile ...\n")
    
    t0 = time.time()
    ampl.read('../modfiles/'+modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    ampl.eval("option display_precision 0;")
    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)
    
    if True:
        ampl.setOption('presolve',0)
        log.joint(' AMPL presolve off\n')
        
    log.joint(" solver set to " + solver + "\n")
    
    if all_data['solver'] == 'gurobi_ampl':
        ampl.eval("options option gurobi_options 'method=2';")
        ampl.eval("options option gurobi_options 'bar=1';")

    if all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';")
        elif all_data['multistart']:
            ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';")
        elif all_data['knitropresolveoff']:
            ampl.eval("option knitro_options 'presolve=0';")            
        else:
            ampl.eval("option knitro_options 'blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';") 
                

    if all_data['fix'] or all_data['initial_point']:
        getsol(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']
        mp_angle = all_data['mp_angle']
        

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
    theta_min    = {}
    theta_max    = {}
    Vinit        = {}
    thetainit    = {}


    for bus in all_data['buses'].values():
        buscount = bus.count

        if ( len(bus.frombranchids) == 0 ) and ( len(bus.tobranchids) == 0 ): #else we could have just skipped bus.Gs and bus.Bs
            log.joint(' isolated bus ' + str(bus.nodeID) + ' busid ' + str(buscount) + '\n') #we are not skipping it
        #    breakexit('check')
        #    continue
        
        buses[buscount] = bus.nodeID 
        Pd[buscount] = bus.Pd
        Qd[buscount] = bus.Qd

        if all_data['fix']:
            Vmin[buscount] = mp_vm[buscount] - tolerance
            Vmax[buscount] = mp_vm[buscount] + tolerance
        else:
            Vmin[buscount] = bus.Vmin
            Vmax[buscount] = bus.Vmax

        
        if all_data['initial_point'] == 0:
            Vinit[buscount]      = 1
            thetainit[buscount]  = 0

        branches_f[buscount] = []
        branches_t[buscount] = []

        bus_gens[buscount] = bus.genidsbycount #load this as a sparse ds

        if buscount == all_data['refbus']:
            theta_min[buscount] = 0
            theta_max[buscount] = 0
        else:
            if all_data['fix']:
                theta_min[buscount] = mp_angle[buscount] - tolerance
                theta_max[buscount] = mp_angle[buscount] + tolerance
            else:
                theta_min[buscount] = -2 * math.pi
                theta_max[buscount] = 2 * math.pi


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
    ampl.get_parameter('theta_min').setValues(theta_min)
    ampl.get_parameter('theta_max').setValues(theta_max)

    if all_data['initial_point']:
        ampl.get_parameter('Vinit').set_values(mp_vm)
        ampl.get_parameter('thetainit').set_values(mp_angle)
    else:
        ampl.get_parameter('Vinit').set_values(Vinit)
        ampl.get_parameter('thetainit').set_values(thetainit)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set = ampl.getSet('bus_gens')
    
    for bus in all_data['buses'].values():
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

    
    for branch in all_data['branches'].values():
        branchcount            = branch.count

        if branch.status != 1: #the reader does not consider branches whose status is 0, hence this should never be true
            log.joint(' branch ' + str(branchcount) + ' OFF, we skip it\n')
            breakexit('check')
            continue

        branches[branchcount]  = (branch.id_f,branch.id_t)
        U[branchcount]         = branch.limit
        Gtt[branchcount]       = branch.Gtt
        Btt[branchcount]       = branch.Btt
        Gff[branchcount]       = branch.Gff
        Bff[branchcount]       = branch.Bff
        Gtf[branchcount]       = branch.Gtf
        Btf[branchcount]       = branch.Btf
        Gft[branchcount]       = branch.Gft
        Bft[branchcount]       = branch.Bft
        bus_f[branchcount]     = branch.id_f
        bus_t[branchcount]     = branch.id_t

            
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
        gencount            = gen.count 
        gens[gencount]      = gen.nodeID
        Pmax[gencount]      = gen.Pmax * gen.status
        Pmin[gencount]      = gen.Pmin * gen.status
        Qmax[gencount]      = gen.Qmax * gen.status
        Qmin[gencount]      = gen.Qmin * gen.status

        buscount = IDtoCountmap[gen.nodeID]
        bus      = all_data['buses'][buscount]
        if bus.nodetype == 3:
            Qmax[gencount] = 2 * all_data['summaxgenQ']
            Qmin[gencount] = - 2 * all_data['summaxgenQ']
    
        
        fixedcost[gencount] = gen.costvector[2] * gen.status
        lincost[gencount]   = gen.costvector[1]
        quadcost[gencount]  = gen.costvector[0]
        
        #print(" cost vector",gen.costvector) #check small case5.m
        #print(" cost degree",gen.costdegree)
        
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
    
    expand = False
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
    v          = ampl.get_variable("v")
    theta      = ampl.get_variable("theta")
    Pf         = ampl.get_variable("Pf")
    Pt         = ampl.get_variable("Pt")
    Qf         = ampl.get_variable("Qf")
    Qt         = ampl.get_variable("Qt")
    Pg         = ampl.get_variable("Pg")
    Qg         = ampl.get_variable("Qg")
    dic_v      = v.get_values().to_dict()
    dic_theta  = theta.get_values().to_dict()
    dic_Pg     = Pg.get_values().to_dict()
    dic_Qg     = Qg.get_values().to_dict()
    dic_Pf     = Pf.get_values().to_dict()
    dic_Pt     = Pt.get_values().to_dict()
    dic_Qf     = Qf.get_values().to_dict()
    dic_Qt     = Qt.get_values().to_dict()
    dicPLoss   = {}
    dicQLoss   = {}
    
    for branch in dic_Pf.keys():
        dicPLoss[branch] = dic_Pf[branch] + dic_Pt[branch] 
        dicQLoss[branch] = dic_Qf[branch] + dic_Qt[branch] 

    objval    = total_cost.get().value()
    timesofar = time.time() - all_data['T0']
        
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" objective " + str(objval) + '\n')
    log.joint(" active power generation " + str(sum(dic_Pg.values())) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(dicPLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(dic_Qg.values())) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(dicQLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')
    
    log.joint(" writing casename, modfile, obj and runtime to summary_ac.log\n")

    summary_ac = open("summary_ac.log","a+") #later add feasibility errors, etc

    summary_ac.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile'] +  ' obj ' + str(objval) + ' runtime ' + str(timesofar) + '\n')

    summary_ac.close()


    log.joint(" ------------------------\n")

    if False:
        topfg_name = 'top_flows_gens/' + all_data['casename'] + '_fg.txt' 
        topfg      = open(topfg_name,"w")

        topflows_f = sorted(dic_Pf.items(), key = lambda x: math.fabs(x[1]), reverse = True)[:10]
        topflows_g = sorted(dic_Pt.items(), key = lambda x: math.fabs(x[1]), reverse = True)[:10]
        topflows   = {}

        branchid_f   = int(topflows_f[0][0])
        flow_f       = topflows_f[0][1]
        branchid_t   = int(topflows_g[0][0])
        flow_t       = topflows_g[0][1]
        i_f          = 0
        i_t          = 0

        while True:
            if math.fabs(flow_f) >= math.fabs(flow_t):
                topflows[branchid_f] = flow_f
                if i_f == 9:
                    break
                i_f         += 1
                branchid_f   = int(topflows_f[i_f][0])
                flow_f       = topflows_f[i_f][1]
            else:
                topflows[branchid_t] = flow_t
                if i_t == 9:
                    break
                i_t         += 1
                branchid_t   = int(topflows_g[i_t][0])
                flow_t       = topflows_g[i_t][1]


        log.joint(' top 10 flows are:\n')
        topfg.write(' top 10 flows are:\n')

        for branchid in topflows.keys():
            f          = all_data['branches'][branchid].f
            t          = all_data['branches'][branchid].t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            bus_f      = all_data['buses'][count_of_f]
            bus_t      = all_data['buses'][count_of_t]

            log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' flow ' + str(topflows[branchid]) + ' fdegree ' + str(bus_f.degree) + ' tdegree ' + str(bus_t.degree) + '\n')
            topfg.write(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' flow ' + str(topflows[branchid]) + ' fdegree ' + str(bus_f.degree) + ' tdegree ' + str(bus_t.degree) + '\n')

        top_active_gens  = dict(sorted(dic_Pg.items(), key = lambda x: x[1], reverse = True)[:10]) 

        log.joint(' top 10 active power gens are:\n')
        topfg.write(' top 10 active power gens are:\n')

        for genid in top_active_gens.keys():
            log.joint(' genid ' + str(int(genid)) + ' at bus ' + str(gens[genid]) + ' active power generation ' + str(top_active_gens[genid]) + '\n')
            topfg.write(' genid ' + str(int(genid)) + ' at bus ' + str(gens[genid]) + ' active power generation ' + str(top_active_gens[genid]) + '\n')

        top_reactive_gens  = dict(sorted(dic_Qg.items(), key = lambda x: x[1], reverse = True)[:10]) 

        log.joint(' top 10 reactive power gens are:\n')
        topfg.write(' top 10 reactive power gens are:\n')

        for genid in top_reactive_gens.keys():
            log.joint(' genid ' + str(int(genid)) + ' at bus ' + str(gens[genid]) + ' active power generation ' + str(top_reactive_gens[genid]) + '\n')
            topfg.write(' genid ' + str(int(genid)) + ' at bus ' + str(gens[genid]) + ' active power generation ' + str(top_reactive_gens[genid]) + '\n')

        topfg.close()

        log.joint(" ------------------------\n")
    
    #breakexit("write solution?")

    branches_data = all_data['branches']
    buses_data    = all_data['buses']
    IDtoCountmap  = all_data['IDtoCountmap']
    tolerance     = 1e-05
    
    filename      = 'ksol_' + all_data['casename'] + '.txt'
    filenamelp    = 'ksol_' + all_data['casename'] + '.lp'
    #filename      = 'sols/ksol_' + all_data['casename'] + '.txt'
    #filenamelp    = 'sols/ksol_' + all_data['casename'] + '.lp'    
    thefile       = open(filename,'w+')
    thefilelp     = open(filenamelp,'w+')
    
    log.joint(' writing solution to ' + filename + '\n')

    thefile.write('value ' + str(total_cost.get().value()) + '\n')
    
    thefile.write('voltages and angles:\n')

    for buscount in buses.keys():
        vval       = dic_v[buscount]
        f          = buses_data[buscount].nodeID
        thetaval   = dic_theta[buscount]
        line = 'bus ' + str(buscount) + ' M ' + str(vval) + ' A ' + str(thetaval) + '\n'
        thefile.write(line)
        
        v2val = vval*vval
        linelp_v = str(v2val - tolerance) + ' <= c_' + str(f) + '_' + str(f) + ' <= ' + str(v2val + tolerance) + '\n'
        thefilelp.write(linelp_v)
    
    thefile.write('power flows:\n')
    
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
        thetaval   = dic_theta[count_of_f] - dic_theta[count_of_t]                                            
        
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' theta ' + str(thetaval) + '\n'
        thefile.write(line)

        linelp_Pf = str(Pfval-tolerance) + ' <= P_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Pfval+tolerance) + '\n'
        thefilelp.write(linelp_Pf)
        linelp_Pt = str(Ptval-tolerance) + ' <= P_' + str(branchid) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Ptval+tolerance) + '\n'
        thefilelp.write(linelp_Pt)
        linelp_Qf = str(Qfval-tolerance) + ' <= Q_' + str(branchid) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Qfval+tolerance) + '\n'
        thefilelp.write(linelp_Qf)
        linelp_Qt = str(Qtval-tolerance) + ' <= Q_' + str(branchid) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Qtval+tolerance) + '\n'
        thefilelp.write(linelp_Qt)
        

    thefile.close()
    thefilelp.close()

    log.joint(' done writing knitro solution to file\n\n')

def getsol_knitro(log,all_data):

    casefilename  = all_data['casename'] 
    #filename      = 'ksol_' + casename + '.txt'
    filename      = 'ksol_' + casename + '.txt'
    print(filename)
    
    try:
        thefile  = open(filename, "r")
        lines    = thefile.readlines()
        lenlines = len(lines)
        thefile.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")
    
    log.joint(' reading voltages\n')

    buses_data    = all_data['buses']
    branches_data = all_data['branches']
    IDtoCountmap  = all_data['IDtoCountmap']
    
    mp_vm    = {}
    mp_angle = {}
    mp_Pfvalues = {}
    mp_Ptvalues = {}
    mp_Qfvalues = {}
    mp_Qtvalues = {}

    linenum  = 2

    while linenum < lenlines:
        thisline = lines[linenum].split()

        if thisline[0] == 'bus':
            busid              = int(thisline[1])
            buscount           = IDtoCountmap[busid]
            mp_vm[buscount]    = float(thisline[3])
            mp_angle[buscount] = float(thisline[5]) * math.pi / 180
        elif thisline[0] == 'branch':
            branchcount     = int(thisline[1])
            mp_Pfvalues[branchcount] = float(thisline[7])/all_data['baseMVA']
            mp_Ptvalues[branchcount] = float(thisline[9])/all_data['baseMVA']
            mp_Qfvalues[branchcount] = float(thisline[11])/all_data['baseMVA']
            mp_Qtvalues[branchcount] = float(thisline[13])/all_data['baseMVA']

        linenum += 1

    print(mp_vm)
    all_data['mp_vm']      = mp_vm
    all_data['mp_angle']   = mp_angle          
        
    all_data['mp_Pfvalues'] = mp_Pfvalues
    all_data['mp_Ptvalues'] = mp_Ptvalues
    all_data['mp_Qfvalues'] = mp_Qfvalues
    all_data['mp_Qtvalues'] = mp_Qtvalues


def getsol(log,all_data):
    
    casename = all_data['casename']
    filename = 'mp_sols/solution_va_'+ casename +'.txt'
    
    log.joint(" reading file matpower solution volt magnitudes and angles " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    lenlines      = len(lines)
    mp_vm         = {}
    mp_angle      = {}
    IDtoCountmap  = all_data['IDtoCountmap']
    linenum       = 1

    while linenum < lenlines:
        thisline = lines[linenum].split(',')
        bus_id   = int(thisline[0])
        buscount = IDtoCountmap[bus_id]

        mp_vm[buscount]        = float(thisline[1])
        mp_angle[buscount]     = float(thisline[2]) * math.pi / 180

        log.joint('bus ' + str(bus_id) + ' v ' + str(mp_vm[buscount]) + ' angle ' + str(mp_angle[buscount]) + '\n')

        linenum += 1
    
    all_data['mp_vm']    = mp_vm
    all_data['mp_angle'] = mp_angle
