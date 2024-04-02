###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                        mjv2153@columbia.edu                               ##
##                                                                           ##
## Oct 2023                                                                  ##
###############################################################################

from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math
import os
import platform

def goac_mtp(log,all_data):

    #readcuts(log,all_data)
    #breakexit('check')
    
    log.joint(' creating ampl object ...\n')
    
    ampl         = AMPL()
    casename     = all_data['casename']
    casetype     = all_data['casetype']
    modfile      = all_data['modfile']
    solver       = all_data['solver']
    IDtoCountmap = all_data['IDtoCountmap']
    T            = all_data['T']
    

    log.joint(" reading modfile ...\n")
    t0 = time.time()
    ampl.read('../modfiles/'+modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    
    # Setting up AMPL options
    ampl.setOption('display_precision',0)
    ampl.setOption('expand_precision',0)
    ampl.setOption('show_stats',2)
    ampl.setOption('solver',solver)    

    log.joint(" solver set to " + solver + "\n")
    
    if all_data['AMPL_presolve']:
        ampl.setOption('presolve',1)
        log.joint(' AMPL presolve on\n')
    else:
        ampl.setOption('presolve',0)        
        log.joint(' AMPL presolve off\n')
        
    # Solver options
    feastol        = all_data['feastol']
    opttol         = all_data['opttol']
    linsolver      = all_data['linsolver']
    max_time       = all_data['max_time']
    wstart         = all_data['wstart']
    bar_initmu     = all_data['bar_initmu']
    
    solver_options = "option knitro_options 'feastol_abs=" + feastol + " opttol_abs=" + opttol + " blasoptionlib=1 numthreads=20 linsolver=" + linsolver + " maxtime_real=" + max_time + " strat_warm_start=" + wstart + " bar_initmu=" + bar_initmu + "';"
    
    ampl.eval(solver_options)


    if all_data['multistart']:
        solver_options = "option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=5 maxtime_real=2400 ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';"
        ampl.eval(solver_options)

    # Initial solution (flat-start or previously computed AC solutions)
    Vinit_mtp, thetadiffinit_mtp = initial_solution_mtp(log,all_data)        


    # Setting up buses
    buses        = {}
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

    
    # Getting multi-time period loads and ramp rates from file
    if all_data['T'] == 1:
        Pd             = all_data['Pd'] = getloads_check(log,all_data)
    else:
        Pd             = all_data['Pd'] = getloads(log,all_data)
        
    rampru, ramprd = getrampr(log,all_data)

    
    for bus in all_data['buses'].values():
        buscount = bus.count

        if ( len(bus.frombranchids) == 0 ) and ( len(bus.tobranchids) == 0 ): 
            log.joint(' isolated bus ' + str(bus.nodeID) + ' busid ' + str(buscount) + '\n')
        
        buses[buscount]      = bus.nodeID
        Qd[buscount]         = bus.Qd
        Vmin[buscount]       = bus.Vmin
        Vmax[buscount]       = bus.Vmax
        branches_f[buscount] = []
        branches_t[buscount] = []
        bus_gens[buscount]   = bus.genidsbycount 
        
        if buscount == all_data['refbus']:
            theta_min[buscount] = 0
            theta_max[buscount] = 0
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

    ampl.get_parameter('T').set(T)            
    ampl.getSet('buses').setValues(list(buses))
    ampl.getSet('bus_Gs').setValues(list(bus_Gs))
    ampl.getSet('bus_Bs').setValues(list(bus_Bs))    
    ampl.get_parameter('Gs').setValues(list(bus_Gs.values()))
    ampl.get_parameter('Bs').setValues(list(bus_Bs.values()))
    ampl.get_parameter('Vmax').set_values(Vmax)
    ampl.get_parameter('Vmin').set_values(Vmin)
    ampl.get_parameter('theta_min').setValues(theta_min)
    ampl.get_parameter('theta_max').setValues(theta_max)
    ampl.get_parameter('Pd').set_values(Pd)    
    ampl.get_parameter('Qd').set_values(Qd)
    ampl.get_parameter('Vinit').set_values(Vinit_mtp)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set   = ampl.getSet('bus_gens')
    
    for bus in all_data['buses'].values():
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
        bus_gens_set[bus.count].setValues(bus_gens[bus.count])

        
    # Setting up branches
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
    maxangle  = {}
    minangle  = {}
    
    for branch in all_data['branches'].values():
        branchcount            = branch.count

        if branch.status != 1:
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
        maxangle[branchcount]  = branch.maxangle_rad
        minangle[branchcount]  = branch.minangle_rad
                    
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
    ampl.get_parameter('maxangle').setValues(maxangle)
    ampl.get_parameter('minangle').setValues(minangle)
    ampl.get_parameter('thetadiffinit').setValues(thetadiffinit_mtp)

    # Setting up generators
    gens      = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} 

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
        
    ampl.getSet('gens').setValues(list(gens))
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)
    ampl.get_parameter('rampru').setValues(rampru)
    ampl.get_parameter('ramprd').setValues(ramprd)

    log.joint(" sets and parameters loaded\n")

    
    # Expand
    if all_data['expand']:
        filename = casename + "_" + casetype + "_NLP.out"
        log.joint('Now expanding to %s.\n'%(filename))
        amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb,_var[j].ub);' 
        modelout = ampl.getOutput(amplstate)
        outfile = open(filename,"w")
        outfile.write("model = " + str(modelout) + "\n")
        outfile.close()

        
    # Solve
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()


    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")
    
    # Get solution
    solver_status = int(ampl.getOutput("display solve_result_num;").split(' ')[2].strip())
    
    total_costvar           = ampl.get_objective("total_cost")
    vvar                    = ampl.get_variable("v")
    thetavar                = ampl.get_variable("theta")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['vvalues']     = vvar.get_values().to_dict()
    all_data['thetavalues'] = thetavar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    sumPd = 0
    for buscount,h in Pd.keys():
        sumPd += Pd[buscount,h]

    sumQd = 0
    for buscount in Qd.keys():
        sumQd += Qd[buscount]

    sumQd = sumQd * T
    
    PLoss    = {}
    QLoss    = {}
    sumPLoss = 0
    sumQLoss = 0

    for branchcount,h in all_data['Pfvalues'].keys():
        PLoss[branchcount,h] = all_data['Pfvalues'][branchcount,h] + all_data['Ptvalues'][branchcount,h]
        QLoss[branchcount,h] = all_data['Qfvalues'][branchcount,h] + all_data['Qtvalues'][branchcount,h]
        sumPLoss += PLoss[branchcount,h]
        sumQLoss += QLoss[branchcount,h]

    sumGenP = 0
    sumGenQ = 0
    for gencount,h in all_data['GenPvalues'].keys():
        sumGenP += all_data['GenPvalues'][gencount,h]
        sumGenP += all_data['GenQvalues'][gencount,h]
    
    timesofar = time.time() - all_data['T0']        
        
    log.joint(" case " + casename + " time-periods " + str(T)
              + " " + casetype + "\n")
    log.joint(' solver ' + all_data['solver'] + ' status '
              + str(solver_status) + ' warmstart '
              + str(all_data['wstart']) + '\n')
    log.joint(" objective " + str(all_data['objvalue']) + '\n')    
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" active power generation " + str(sumGenP) + '\n')
    log.joint(" active power demand " + str(sumPd) + '\n')
    log.joint(" active power loss " + str(sumPLoss) + '\n')
    log.joint(" reactive power generation " + str(sumGenQ) + '\n')
    log.joint(" reactive power demand " + str(sumQd) + '\n')
    log.joint(" reactive power loss " + str(sumQLoss) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')    
    
    log.joint(" writing casename, modfile, obj and runtime to summary_ac.log\n")

    summary_ac = open("summary_ac.log","a+")

    summary_ac.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile']
                     + ' T ' + str(all_data['T']) + ' solver ' + all_data['solver'] + ' obj '
                     + str(all_data['objvalue']) + ' runtime ' + str(timesofar) + '\n')

    summary_ac.close()

    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")
    
    if all_data['writesol']:
        writesol(log,all_data)
        writesol_qcqp_allvars(log,all_data)

    log.joint('\n solver status ' + str(solver_status) + ' Heuristic 1\n')        
    if solver_status:        
        return 1
    else:
        return 0
        
def writesol(log,all_data):

    branches     = all_data['branches']
    buses        = all_data['buses']
    gens         = all_data['gens']
    IDtoCountmap = all_data['IDtoCountmap']
    T            = all_data['T']
    tolerance    = 1e-05

    objvalue     = all_data['objvalue']
    vvalues      = all_data['vvalues']
    thetavalues  = all_data['thetavalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']
    

    datenow       = '_01_28_24'
    
    if all_data['uniform']:
        filename      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_u' + str(all_data['uniform_drift']) + '.txt'
    elif all_data['uniform2']:
        filename      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_u2.txt'        
    elif all_data['nperturb']:
        filename      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_n1.txt'

    thefile       = open(filename,'w+')

    log.joint(' writing solution to ' + filename + '\n')

    machinename    = "cool1"
    now            = time.time()
    AMPL_version   = 'Version 20231012'
    solver_version = 'Artelys Knitro 13.2.0'
    opsystem       = 'Fedora 34 (Workstation Edition)'
    processor      = 'Intel(R) Xeon(R) Linux64 CPU E5-2687W v3 3.10GHz'
    cores          = '20 physical cores, 40 logical processors'
    ram            = '256 GB RAM'
    
    
    thefile.write('/ACsolution : ' + all_data['casename'] + '\n')
    thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefile.write('/MachineName : ' + machinename + '\n')
    thefile.write('/Processor : ' + processor + '\n')
    thefile.write('/OS : ' + opsystem + '\n')
    thefile.write('/Cores : ' + cores + '\n')
    thefile.write('/RAM : ' + ram + '\n')
    thefile.write('/AMPL : ' + AMPL_version + '\n')
    thefile.write('/Solver : ' + solver_version + '\n')
    thefile.write('objvalue ' + str(all_data['objvalue']) + '\n')
    thefile.write('time-periods ' + str(T) + '\n')
    
    thefile.write('voltages and angles:\n')

    for buscount in buses.keys():
        for k in range(T):
            vval       = vvalues[buscount,k]
            f          = buses[buscount].nodeID
            thetaval   = thetavalues[buscount,k]
            line = 'bus ' + str(buscount) + ' M ' + str(vval) + ' A ' + str(thetaval) + ' k ' + str(k) + '\n'
            thefile.write(line)

    thefile.write('power flows and cs variables:\n')
    
    for branchid in branches.keys():
        for k in range(T):
            branch     = branches[branchid]
            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            Pfval      = Pfvalues[branchid,k]
            Ptval      = Ptvalues[branchid,k]
            Qfval      = Qfvalues[branchid,k]
            Qtval      = Qtvalues[branchid,k]
            thetaval   = thetavalues[count_of_f,k] - thetavalues[count_of_t,k]
            cftval     = vvalues[count_of_f,k] * vvalues[count_of_t,k] * math.cos(thetaval)
            sftval     = vvalues[count_of_f,k] * vvalues[count_of_t,k] * math.sin(thetaval)        
        
            line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' k ' + str(k) + '\n'
            thefile.write(line)

    thefile.write('generation:\n')
        
    for genid in gens.keys():
        for k in range(T):
            gen     = gens[genid] 
            nodeID  = gen.nodeID
            line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPvalues[genid,k]) + ' GQ ' + str(GenQvalues[genid,k]) + ' k ' + str(k) + '\n'
            thefile.write(line_gen)
        
        
    thefile.close()

    log.joint(' done writing knitro solution to .txt file\n\n')

def writesol_qcqp_allvars(log,all_data):

    branches      = all_data['branches']
    buses         = all_data['buses']
    gens          = all_data['gens']
    IDtoCountmap  = all_data['IDtoCountmap']
    T             = all_data['T']
    tolerance     = 1e-05

    objvalue     = all_data['objvalue']
    vvalues      = all_data['vvalues']
    thetavalues  = all_data['thetavalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']
    Pd           = all_data['Pd']


    datenow       = '_01_28_24'

    if all_data['uniform']:
        filenamevars      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_u' + str(all_data['uniform_drift']) + '.sol'
    elif all_data['uniform2']:
        filenamevars      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_u2.sol'        
    elif all_data['nperturb']:
        filenamevars      = 'ACsol_' + all_data['casename'] + '_T' + str(T) + '_n1.sol'                                                                                     
    thefilevars   = open(filenamevars,'w+')
    
    log.joint(' writing solution to ' + filenamevars + '\n')

    machinename    = "cool1"
    now            = time.time()
    AMPL_version   = 'Version 20231012'
    solver_version = 'Artelys Knitro 13.2.0'
    opsystem       = 'Fedora 34 (Workstation Edition)'
    processor      = 'Intel(R) Xeon(R) Linux64 CPU E5-2687W v3 3.10GHz'
    cores          = '20 physical cores, 40 logical processors'
    ram            = '256 GB RAM'
    
    thefilevars.write('/ACsolution : ' + all_data['casename'] + '\n')
    thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefilevars.write('/MachineName : ' + machinename + '\n')
    thefilevars.write('/Processor : ' + processor + '\n')
    thefilevars.write('/OS : ' + opsystem + '\n')
    thefilevars.write('/Cores : ' + cores + '\n')
    thefilevars.write('/RAM : ' + ram + '\n')
    thefilevars.write('/AMPL : ' + AMPL_version + '\n')
    thefilevars.write('/Solver : ' + solver_version + '\n')
    thefilevars.write('/Objvalue : ' + str(all_data['objvalue']) + '\n')
    thefilevars.write('/Time-periods ' + str(T) + '\n')
    
    for buscount in buses.keys():
        for k in range(T):
            bus        = buses[buscount]
            vval       = vvalues[buscount,k]
            f          = bus.nodeID
            thetaval   = thetavalues[buscount,k]
            v2value    = vval**2
            evalue     = vval * math.cos(thetaval)
            fvalue     = vval * math.sin(thetaval)

            v2name     = 'c_' + str(f) + '_' + str(f) + '_' + str(k)
            ename      = 'e_' + str(f) + '_' + str(k)
            fname      = 'f_' + str(f) + '_' + str(k)

            v2line     = v2name + ' = ' + str(v2value) + '\n'
            thefilevars.write(v2line)

            eline     = ename + ' = ' + str(evalue) + '\n'
            thefilevars.write(eline)

            fline     = fname + ' = ' + str(fvalue) + '\n'
            thefilevars.write(fline)        

            IPvalue = - Pd[buscount,k]
            IQvalue = - bus.Qd
            IPname  = 'IP_' + str(bus.nodeID) + '_' + str(k)
            IQname  = 'IQ_' + str(bus.nodeID) + '_' + str(k)

            for gencounter in bus.genidsbycount:
                if gens[gencounter].status:
                    IPvalue += GenPvalues[gencounter,k]
                    IQvalue += GenQvalues[gencounter,k]

            IPline = IPname + ' = ' + str(IPvalue) + '\n'
            thefilevars.write(IPline)

            IQline = IQname + ' = ' + str(IQvalue) + '\n'
            thefilevars.write(IQline)
        
        
    for branchid in branches.keys():
        for k in range(T):
            branch     = branches[branchid]
            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            thetaftval = thetavalues[count_of_f,k] - thetavalues[count_of_t,k]
            Pfval      = Pfvalues[branchid,k]
            Ptval      = Ptvalues[branchid,k]
            Qfval      = Qfvalues[branchid,k]
            Qtval      = Qtvalues[branchid,k]
            cftval     = vvalues[count_of_f,k] * vvalues[count_of_t,k] * math.cos(thetaftval)
            sftval     = vvalues[count_of_f,k] * vvalues[count_of_t,k] * math.sin(thetaftval)
            
            Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
            Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
            Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k) 
            Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
            cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
            sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)

            Pfline  = Pfname + ' = ' + str(Pfval) + '\n'
            thefilevars.write(Pfline)

            Ptline  = Ptname + ' = ' + str(Ptval) + '\n'
            thefilevars.write(Ptline)

            Qfline  = Qfname + ' = ' + str(Qfval) + '\n'
            thefilevars.write(Qfline)

            Qtline  = Qtname + ' = ' + str(Qtval) + '\n'
            thefilevars.write(Qtline)

            cftline  = cftname + ' = ' + str(cftval) + '\n'
            thefilevars.write(cftline)

            sftline  = sftname + ' = ' + str(sftval) + '\n'
            thefilevars.write(sftline)                        

        
    for genid in gens.keys():
        for k in range(T):
            gen     = gens[genid] 
            nodeID  = gen.nodeID
            GenPval = GenPvalues[genid,k]
            GenQval = GenQvalues[genid,k]
            GPname  = 'GP_' + str(genid) + '_' + str(nodeID) + '_' + str(k)
            GQname  = 'GQ_' + str(genid) + '_' + str(nodeID) + '_' + str(k)

            GPline  = GPname + ' = ' + str(GenPval) + '\n'
            thefilevars.write(GPline)

            GQline  = GQname + ' = ' + str(GenQval) + '\n'
            thefilevars.write(GQline)

        
    log.joint(' done writing AC QCQP allvars solution to .sol file\n\n')
        
    thefilevars.close()


def topfg(log,all_data): #check this

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
    

def readcuts(log,all_data):

    branches     = all_data['branches']
    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']
    
    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif ('_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename'] 
          or '_n_1_1' in all_data['casename']):
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    elif '_line' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 5]
    #elif '_line5' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = '../cuts/cuts_' + original_casename + '.txt' ###cuts | newcuts
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    firstline = lines[1].split()
    jabr      = 1
    linenum   = 2

    jabrs  = {}
    i2s    = {}
    limits = {}

    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            numi2cuts = int(thisline[3])
            all_data['addcuts_numi2cuts'] = numi2cuts
            log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
            linenum += 1
            jabr = 0
            i2   = 1
            continue

        elif thisline[0] == '#limit-envelope' and i2:
            numlimitcuts = int(thisline[3])
            all_data['addcuts_numlimitcuts'] = numlimitcuts
            log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) 
                      + '\n')
            linenum += 1
            i2       = 0
            continue

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            cutid     = int(thisline[7])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid,f,t,branch.r) not in jabrs:
                jabrs[(branchid,f,t,branch.r)] = 1
            else:
                jabrs[(branchid,f,t,branch.r)] += 1
            
            linenum += 1

        elif (jabr == 0) and i2:

            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid,f,t,branch.r) not in i2s:
                i2s[(branchid,f,t,branch.r)] = 1
            else:
                i2s[(branchid,f,t,branch.r)] += 1
            
            linenum += 1


        elif (jabr == 0) and (i2 == 0):
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            if thisline[14] == 'Pft':
                from_or_to = 'f'
            elif thisline[14] == 'Ptf':
                from_or_to = 't'
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            branch = branches[branchid]

            if (branchid,f,t,branch.r) not in limits:
                limits[(branchid,f,t,branch.r)] = 1
            else:
                limits[(branchid,f,t,branch.r)] += 1
                
            linenum += 1


    numsel = 1000
            
    mostjabrs  = dict(sorted(jabrs.items(), key = lambda x: x[1], reverse = True)[:numsel])

    mosti2     = dict(sorted(i2s.items(), key = lambda x: x[1], reverse = True)[:numsel])

    mostlimits = dict(sorted(limits.items(), key = lambda x: x[1], reverse = True)[:numsel])

    log.joint(' mostjabrs ' + str(mostjabrs) + '\n')
    #log.joint(' mosti2s ' + str(mostjabrs) + '\n')
    #log.joint(' mostlimits ' + str(mostjabrs) + '\n')

    log.joint(' number of cones with Jabr-cuts ' + str(len(jabrs)) + '\n')

    # num   = list(mostjabrs.keys())[0]
    # count = 0
    # numdic = {}
    # for numcuts in mostjabrs.values():
    #     if m

def getloads(log,all_data):

    casename = all_data['casename']
    T        = all_data['T']

    if all_data['nperturb']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_n1.txt'
    elif all_data['uniform']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_u' + str(all_data['uniform_drift']) + '.txt'
    elif all_data['uniform2']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_u2.txt'

    try:
        thefile = open(filename, "r")
        lines = thefile.readlines()
        lenlines = len(lines) - 1 # END                                                  
        thefile.close()
    except:
        log.stateandquit("Cannot open file " + datafilename)
        sys.exit("Check file with loads")

    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    Pd      = {}
    linenum = 0

    log.joint(' reading file with multi-period demands\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()
        buscount        = int(thisline[1])
        k               = int(thisline[5])
        load            = float(thisline[7]) 
        Pd[buscount,k]  = load
        linenum        += 1

    log.joint(' done reading multi-period demands\n')        
    return  Pd


def getloads_check(log,all_data):

    casename = all_data['casename']
    T        = all_data['T'] + 3

    if all_data['nperturb']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_n1.txt'
    elif all_data['uniform']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_u' + str(all_data['uniform_drift']) + '.txt'
    elif all_data['uniform2']:
        filename = '../../cutplane/newmtploads2/' + casename + '_mtploads_' + str(T) + '_u2.txt'

    try:
        thefile = open(filename, "r")
        lines = thefile.readlines()
        lenlines = len(lines) - 1 # END                                                  
        thefile.close()
    except:
        log.stateandquit("Cannot open file " + datafilename)
        sys.exit("Check file with loads")

    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    Pd      = {}
    linenum = 0

    log.joint(' reading file with multi-period demands\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()
        buscount        = int(thisline[1])
        k               = int(thisline[5])
        if k == 0:
            load            = float(thisline[7]) 
            Pd[buscount,0]  = load
        linenum        += 1

    log.joint(' Pd buscount ' + str(buscount) + ' is ' + str(Pd[buscount,0]) + '\n')
    log.joint(' done reading multi-period demands\n')        
    return  Pd

def getrampr(log,all_data):

    casename = all_data['casename']
    T        = all_data['T']
    filename = '../../cutplane/ramprates/' + casename + '_rampr_' + str(T) + '.txt'
    try:
        thefile = open(filename, "r")
        lines = thefile.readlines()
        lenlines = len(lines) - 1 # END                                                      
        thefile.close()
    except:
        log.stateandquit("Cannot open file " + datafilename)
        sys.exit("Check file with ramp rates")

    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    rampru   = {}
    ramprd   = {}
    linenum = 0

    log.joint(' reading file with multi-period ramp rates\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()
        gencount        = int(thisline[1])
        k               = int(thisline[5])
        rpru            = float(thisline[7])
        rprd            = float(thisline[9])
        rampru[gencount,k] = rpru
        ramprd[gencount,k] = rprd
        linenum        += 1

    log.joint(' done reading multi-period ramp rates\n')
    return  rampru, ramprd


def getsol_ac_mtp(log,all_data):

    T            = all_data['T']
    casename     = all_data['casename']
    casetype     = all_data['casetype']
    branches     = all_data['branches']
    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']
    gens         = all_data['gens']

    Vinit_mtp         = {}
    thetadiffinit_mtp = {}

    filename = 'ACsol_' + all_data['casename'] + '.txt'

    try:
        thefile   = open(filename, "r")
        lines     = thefile.readlines()
        lenlines  = len(lines)
        thefile.close()
    except:
        log.stateandquit('Cannot open file ' + filename + '\n')
        sys.exit("failure")
                
    linenum = 0
    log.joint(' reading file ' + filename + '\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()
        if thisline[0] == 'bus':
            buscount             = int(thisline[1])
            for h in range(T):
                Vinit_mtp[buscount,h] = float(thisline[3])
        elif thisline[0] == 'branch':
            branchcount  = int(thisline[1])
            cval         = float(thisline[15])
            sval         = float(thisline[17])
            thetadiff    = math.atan2(sval,cval)
            for h in range(T):
                thetadiffinit_mtp[branchcount,h] = thetadiff
                
        linenum += 1

    log.joint(' done loading single time-period solution\n')

    return Vinit_mtp, thetadiffinit_mtp


def initial_solution_mtp(log,all_data):

    T = all_data['T']
    
    if all_data['initial_solution']:
        log.joint(' loading initial solution...\n')
        Vinit_mtp, thetadiffinit_mtp = getsol_ac_mtp(log,all_data)
    else:
        log.joint(' flat-start (default)\n')
        Vinit_mtp         = {}
        thetadiffinit_mtp = {}
        for h in range(T):            
            for buscount in all_data['buses'].keys():
                Vinit_mtp[buscount,h] = 1  

            for branchcount in all_data['branches'].keys():
                thetadiffinit_mtp[branchcount,h] = 0

    return Vinit_mtp, thetadiffinit_mtp
