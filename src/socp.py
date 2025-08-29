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
import socket
import psutil

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

    elif all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=5 maxtime_real=1000 convex=1';") #bar_conic_enable=1  honorbnds=1
        elif all_data['multistart']:
            ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';")
        elif all_data['knitropresolveoff']:
            ampl.eval("option knitro_options 'presolve=0';")       
        else:
            # Only add bar_murule when using algorithm 1 (Interior-Direct) or 6 (Augmented Lagrangian)
            algorithm = all_data.get("knitro_algorithm", 1)  # default to 1 if not set
            if algorithm in (1, 6):
                ampl.eval(f"option knitro_options 'algorithm={algorithm} bar_murule=1 blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';")
            else:
                ampl.eval(f"option knitro_options 'algorithm={algorithm} blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';")


    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']

    IDtoCountmap = all_data['IDtoCountmap']
        
    # Setting up buses
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

        bus_gens[buscount] = bus.genidsbycount 
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
            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)
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
    
    # Setting up the generators
    gens      = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} 

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
        
    # Solve
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()

    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")


    # Get solution

    total_costvar           = ampl.get_objective("total_cost")
    v2var                   = ampl.get_variable("v")
    cvar                    = ampl.get_variable("c")
    svar                    = ampl.get_variable("s")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['v2values']    = v2var.get_values().to_dict()
    all_data['cvalues']     = cvar.get_values().to_dict()
    all_data['svalues']     = svar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    
    PLoss   = {}
    QLoss   = {}
    for branch in all_data['Pfvalues'].keys():
        PLoss[branch] = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch] 
        QLoss[branch] = all_data['Qfvalues'][branch] + all_data['Qtvalues'][branch] 

    timesofar = time.time() - all_data['T0']        
    
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" solver " + all_data['solver'] + "\n")
    log.joint(" objective " + str(all_data['objvalue']) + '\n')
    log.joint(" active power generation " + str(sum(all_data['GenPvalues'])) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(PLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(all_data['GenQvalues'])) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(QLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')    

    log.joint(" writing casename, modfile, obj and runtime to summary_socp.log\n")

    summary_socp = open("summary_socp.log","a+")

    summary_socp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile']
                       + ' solver ' + all_data['solver'] +  ' obj ' + str(all_data['objvalue'])
                       + ' runtime ' + str(timesofar) + '\n')

    summary_socp.close()


    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")

    if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)


def writesol(log,all_data):

    branches     = all_data['branches']
    buses        = all_data['buses']
    gens         = all_data['gens']
    IDtoCountmap = all_data['IDtoCountmap']
    tolerance    = 1e-05

    objvalue     = all_data['objvalue']
    v2values     = all_data['v2values']
    cvalues      = all_data['cvalues']
    svalues      = all_data['svalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']

    ampl = AMPL()

    if all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        i2fvalues = all_data['i2fvalues']
        
    if all_data['modfile'] == 'jabr.mod' or all_data['modfile'] == 'jabr_mosek.mod':
        if all_data['solver'] == 'knitroampl':
            filename      = 'JABRsol_knitro_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_knitro_' + all_data['casename'] + '.txt'
        elif all_data['solver'] == 'gurobi_ampl':
            filename      = 'JABRsol_gurobi_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_gurobi_' + all_data['casename'] + '.txt'
        elif all_data['solver'] == 'mosek':
            filename      = 'JABRsol_mosek_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_mosek_' + all_data['casename'] + '.txt'
    elif all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        if all_data['solver'] == 'knitroampl':
            filename      = 'I2sol_knitro_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/I2sol_knitro_' + all_data['casename'] + '.txt'                               
        elif all_data['solver'] == 'gurobi_ampl':
            filename      = 'I2sol_gurobi_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/I2sol_gurobi_' + all_data['casename'] + '.txt'                        
        elif all_data['solver'] == 'mosek':
            filename      = 'I2sol_mosek_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/I2sol_mosek_' + all_data['casename'] + '.txt'            
    elif all_data['modfile'] == 'acopfqcqp.mod':
        if all_data['solver'] == 'knitro':
            filename      = 'ACQCQPsol_knitro_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_knitro_' + all_data['casename'] + '.txt'                               
        elif all_data['solver'] == 'gurobi':
            filename      = 'ACQCQPsol_gurobi_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_gurobi_' + all_data['casename'] + '.txt'                        
        elif all_data['solver'] == 'mosek':
            filename      = 'ACQCQPsol_mosek_' + all_data['casename'] + '.txt'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_mosek_' + all_data['casename'] + '.txt'

                
    thefile       = open(filename,'w+')

    log.joint(' writing solution to ' + filename + '\n')


    machinename    = socket.gethostname()
    now            = time.time()
    full_text      = ampl.get_option("version")
    AMPL_version   = full_text.splitlines()[0]
    if all_data['solver'] == 'knitro':
        solver_version = 'Artelys Knitro 13.2.0'
    elif all_data['solver'] == 'gurobi':
        solver_version = 'Gurobi 12.0.2'
    elif all_data['solver'] == 'mosek':
        solver_version = 'MOSEK 10.0.43'
    opsystem       = f"{platform.system()} {platform.release()} ({platform.platform()})"
    processor      = platform.processor() or platform.machine()
    try:
        ram_gb = f"{round(psutil.virtual_memory().total / (1024**3))} GB"
        cores = f"{psutil.cpu_count(logical=False)} physical cores, {psutil.cpu_count(logical=True)} logical processors"
    except ImportError:
        ram_gb = "Unknown RAM"
        cores = f"{os.cpu_count()} cores"
    ram            = ram_gb


    if all_data['modfile'] == 'jabr.mod' or all_data['modfile'] == 'jabr_mosek.mod':
        thefile.write('/JABRsolution : ' + all_data['casename'] + '\n')
    elif all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        thefile.write('/I2solution : ' + all_data['casename'] + '\n')
    elif all_data['modfile'] == 'acopfqcqp.mod':
        thefile.write('/ACQCQPsolution : ' + all_data['casename'] + '\n')
    
    thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefile.write('/MachineName : ' + machinename + '\n')
    thefile.write('/Processor : ' + processor + '\n')
    thefile.write('/OS : ' + opsystem + '\n')
    thefile.write('/Cores : ' + cores + '\n')
    thefile.write('/RAM : ' + ram + '\n')
    thefile.write(AMPL_version + '\n')
    thefile.write('/Solver : ' + solver_version + '\n')
    thefile.write('objvalue ' + str(all_data['objvalue']) + '\n')
    
    thefile.write('voltages:\n')

    for buscount in buses.keys():
        v2val      = v2values[buscount]
        vval       = (v2val)**(0.5)
        f          = buses[buscount].nodeID
        line = 'bus ' + str(buscount) + ' M ' + str(vval) + '\n'
        thefile.write(line)

        
    thefile.write('power flows and cs variables:\n')
    
    for branchid in branches.keys():

        branch     = branches[branchid]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        Pfval      = Pfvalues[branchid]
        Ptval      = Ptvalues[branchid]
        Qfval      = Qfvalues[branchid]
        Qtval      = Qtvalues[branchid]
        cftval     = cvalues[branchid]
        sftval     = svalues[branchid]

        if all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
            i2fval = i2fvalues[branchid]
            line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' i2ft ' + str(i2fval) + '\n'
        else:
            line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + '\n'

        thefile.write(line)

    thefile.write('generation:\n')
        
    for genid in gens.keys():
        gen     = gens[genid] 
        nodeID  = gen.nodeID
        line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPvalues[genid]) + ' GQ ' + str(GenQvalues[genid]) + '\n'
        thefile.write(line_gen)
        
        
    thefile.close()

    log.joint(' done writing SOCP solution to .txt file\n\n')


def writesol_allvars(log,all_data):

    branches      = all_data['branches']
    buses         = all_data['buses']
    gens          = all_data['gens']
    IDtoCountmap  = all_data['IDtoCountmap']
    tolerance     = 1e-05

    objvalue     = all_data['objvalue']
    v2values     = all_data['v2values']
    cvalues      = all_data['cvalues']
    svalues      = all_data['svalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']

    ampl = AMPL()

    if all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        i2fvalues = all_data['i2fvalues']


    datenow       = '_01_28_24'
        
    if all_data['modfile'] == 'jabr.mod' or all_data['modfile'] == 'jabr_mosek.mod':
        if all_data['solver'] == 'knitroampl':
            filename      = 'JABRsol_knitro_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_knitro_' + all_data['casename'] + '.sol'
        elif all_data['solver'] == 'gurobi_ampl':
            filename      = 'JABRsol_gurobi_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_gurobi_' + all_data['casename'] + '.sol'
        elif all_data['solver'] == 'mosek':
            filename      = 'JABRsol_mosek_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/JABRsol_mosek_' + all_data['casename'] + '.sol'
    elif all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        if all_data['solver'] == 'knitroampl':
            filename      = 'I2sol_knitro_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/I2sol_knitro_' + all_data['casename'] + '.sol'                               
        elif all_data['solver'] == 'gurobi_ampl':
            filename      = 'I2sol_gurobi_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/I2sol_gurobi_' + all_data['casename'] + '.sol'                        
        elif all_data['solver'] == 'mosek':
            filename      = 'I2sol_mosek_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/I2sol_mosek_' + all_data['casename'] + '.sol'            
    elif all_data['modfile'] == 'acopfqcqp.mod':
        if all_data['solver'] == 'knitro':
            filename      = 'ACQCQPsol_knitro_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_knitro_' + all_data['casename'] + '.txt'                               
        elif all_data['solver'] == 'gurobi':
            filename      = 'ACQCQPsol_gurobi_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_gurobi_' + all_data['casename'] + '.txt'                        
        elif all_data['solver'] == 'mosek':
            filename      = 'ACQCQPsol_mosek_' + all_data['casename'] + '.sol'
            #filename      = 'SOCPsols' + datenow + '/ACQCQPsol_mosek_' + all_data['casename'] + '.txt'
    
    thefilevars   = open(filename,'w+')


    log.joint(' writing solution to ' + filename + '\n')
    
    machinename    = socket.gethostname()
    now            = time.time()
    full_text      = ampl.get_option("version")
    AMPL_version   = full_text.splitlines()[0]
    if all_data['solver'] == 'knitro':
        solver_version = 'Artelys Knitro 13.2.0'
    elif all_data['solver'] == 'gurobi':
        solver_version = 'Gurobi 12.0.2'
    elif all_data['solver'] == 'mosek':
        solver_version = 'MOSEK 10.0.43'
    opsystem       = f"{platform.system()} {platform.release()} ({platform.platform()})"
    processor      = platform.processor() or platform.machine()
    try:
        ram_gb = f"{round(psutil.virtual_memory().total / (1024**3))} GB"
        cores = f"{psutil.cpu_count(logical=False)} physical cores, {psutil.cpu_count(logical=True)} logical processors"
    except ImportError:
        ram_gb = "Unknown RAM"
        cores = f"{os.cpu_count()} cores"
    ram            = ram_gb

    if all_data['modfile'] == 'jabr.mod' or all_data['modfile'] == 'jabr_mosek.mod':
        thefilevars.write('/JABRsolution : ' + all_data['casename'] + '\n')
    elif all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
        thefilevars.write('/I2solution : ' + all_data['casename'] + '\n')
    elif all_data['modfile'] == 'acopfqcqp.mod':
        thefilevars.write('/ACQCQPsolution : ' + all_data['casename'] + '\n')
        
    thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefilevars.write('/MachineName : ' + machinename + '\n')
    thefilevars.write('/Processor : ' + processor + '\n')
    thefilevars.write('/OS : ' + opsystem + '\n')
    thefilevars.write('/Cores : ' + cores + '\n')
    thefilevars.write('/RAM : ' + ram + '\n')
    thefilevars.write(AMPL_version + '\n')
    thefilevars.write('/Solver : ' + solver_version + '\n')
    thefilevars.write('/Objvalue ' + str(all_data['objvalue']) + '\n')
    
    for buscount in buses.keys():
        bus        = buses[buscount]
        f          = bus.nodeID
        v2value    = v2values[buscount]        

        v2name     = 'c_' + str(f) + '_' + str(f)
        v2line     = v2name + ' = ' + str(v2value) + '\n'
        thefilevars.write(v2line)

        IPvalue = - bus.Pd
        IQvalue = - bus.Qd
        IPname  = 'IP_' + str(bus.nodeID)
        IQname  = 'IQ_' + str(bus.nodeID)

        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                IPvalue += GenPvalues[gencounter]
                IQvalue += GenQvalues[gencounter]

        IPline = IPname + ' = ' + str(IPvalue) + '\n'
        thefilevars.write(IPline)

        IQline = IQname + ' = ' + str(IQvalue) + '\n'
        thefilevars.write(IQline)
        
        
    for branchid in branches.keys():

        branch     = branches[branchid]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        Pfval      = Pfvalues[branchid]
        Ptval      = Ptvalues[branchid]
        Qfval      = Qfvalues[branchid]
        Qtval      = Qtvalues[branchid]
        cftval     = cvalues[branchid]
        sftval     = svalues[branchid]

        Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t)
        Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f)
        Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t)
        Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f)
        cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t)
        sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t)        

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

        if all_data['modfile'] == 'i2.mod' or all_data['modfile'] == 'i2_mosek.mod':
            i2fval  = i2fvalues[branchid]
            iftname = 'i2_' + str(branchid) + '_' + str(f) + '_' + str(t)
            iftline = iftname + ' = ' + str(i2fval) + '\n'

            thefilevars.write(iftline)
        
        
    for genid in gens.keys():
        gen     = gens[genid] 
        nodeID  = gen.nodeID
        GenPval = GenPvalues[genid]
        GenQval = GenQvalues[genid]
        GPname  = "GP_" + str(genid) + "_" + str(nodeID)
        GQname  = "GQ_" + str(genid) + "_" + str(nodeID)

        GPline  = GPname + ' = ' + str(GenPval) + '\n'
        thefilevars.write(GPline)

        GQline  = GQname + ' = ' + str(GenQval) + '\n'
        thefilevars.write(GQline)

        
    log.joint(' done writing SOCP allvars solution to .sol file\n\n')
        
    thefilevars.close()

#update later    
def gosocp2(log,all_data):

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
        ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 TimeLimit=1000 writemodel=i2.lp';")
        #ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 iisfind=1 writemodel=jabr.lp resultfile=jabr.ilp';")

    elif all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=5 maxtime_real=1000 convex=1';") #bar_conic_enable=1
        elif all_data['multistart']:
            ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';")
        elif all_data['knitropresolveoff']:
            ampl.eval("option knitro_options 'presolve=0';")       
        else:
            ampl.eval("option knitro_options 'blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';") 


    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']
	#mp_angle = all_data['mp_angle']
        print(mp_vm)

    IDtoCountmap = all_data['IDtoCountmap']
        
    # Setting up buses
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
    CSmax     = {}
    c_ubound  = {}
    c_lbound  = {}
    s_ubound  = {}
    s_lbound  = {}
    Cinit     = {}
    Sinit     = {}

    
    i2max       = {}
    g           = {}
    b           = {}
    bshunt      = {}
    ratio       = {}
    phase_angle = {}
    
    counter = 0 
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

        y                        = branch.y
        i2max[branchcount]       = branch.limit**2 / Vmin[branch.id_f] * Vmin[branch.id_f]
        g[branchcount]           = y.real
        b[branchcount]           = y.imag
        bshunt[branchcount]      = branch.bc
        ratio[branchcount]       = branch.ratio
        phase_angle[branchcount] = branch.angle_rad

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

            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)
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
        if maxanglerad <= math.pi/2:
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
    
    ampl.get_parameter('i2max').setValues(i2max)    
    ampl.get_parameter('g').setValues(g) 
    ampl.get_parameter('b').setValues(b)
    ampl.get_parameter('bshunt').setValues(bshunt)    
    ampl.get_parameter('ratio').setValues(ratio)
    ampl.get_parameter('phase_angle').setValues(phase_angle)


    # Setting up generators
    gens = {}
    Pmax = {}
    Pmin = {}
    Qmax = {}
    Qmin = {}
    fixedcost = {}
    lincost = {}
    quadcost = {} 

    for gen in all_data['gens'].values():
        gencount = gen.count 
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
        
    # Solve
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()

    
    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")


    # Get solution

    total_costvar           = ampl.get_objective("total_cost")
    v2var                   = ampl.get_variable("v")
    i2fvar                  = ampl.get_variable("i2f")
    cvar                    = ampl.get_variable("c")
    svar                    = ampl.get_variable("s")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['v2values']    = v2var.get_values().to_dict()
    all_data['i2fvalues']   = i2fvar.get_values().to_dict()
    all_data['cvalues']     = cvar.get_values().to_dict()
    all_data['svalues']     = svar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    PLoss   = {}
    QLoss   = {}
    for branch in all_data['Pfvalues'].keys():
        PLoss[branch] = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch] 
        QLoss[branch] = all_data['Qfvalues'][branch] + all_data['Qtvalues'][branch] 

    timesofar = time.time() - all_data['T0']        
    
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" solver " + all_data['solver'] + "\n")
    log.joint(" objective " + str(all_data['objvalue']) + '\n')
    log.joint(" active power generation " + str(sum(all_data['GenPvalues'])) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(PLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(all_data['GenQvalues'])) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(QLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')    

    log.joint(" writing casename, modfile, obj and runtime to summary_socp.log\n")

    summary_socp = open("summary_socp.log","a+")

    summary_socp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile'] +  ' obj ' + str(all_data['objvalue']) + ' runtime ' + str(timesofar) + '\n')

    summary_socp.close()

    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")

    if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)


def gosocp_mosek(log,all_data):

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
    ampl.eval("option mosek_options 'tech:optionnativeread=/opt/newampl/mosekopt outlev=1 sol:chk:feastol=1e-08 chk:mode=2';")
    #cvt:socp=0     

    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']

    IDtoCountmap = all_data['IDtoCountmap']
        
    # Setting up buses
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

        bus_gens[buscount] = bus.genidsbycount 
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
            
            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)   
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
        if maxanglerad <= math.pi/2:
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
    
    # Setting up generators
    gens      = {}
    lingens   = {}
    quadgens  = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} 
    alfa      = {}
    beta      = {}
    gamma     = {}
    sigma     = 0

    for gen in all_data['gens'].values():
        gencount = gen.count 
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


        fixedcost[gencount] = gen.costvector[2] * gen.status
        lincost[gencount] = gen.costvector[1]
        quadcost[gencount] = gen.costvector[0]
        
        if gen.costvector[0]:
            quadgens[gencount] = gen.nodeID
            alfa[gencount]     = math.sqrt(gen.costvector[0])            
            beta[gencount]     = - lincost[gencount] / (2 * alfa[gencount])
            gamma[gencount]    = fixedcost[gencount] - beta[gencount]**2
            sigma             += gamma[gencount]
            #log.joint(' alfa ' + str(alfa[gencount]) + ' beta ' + str(beta[gencount]) + ' gamma ' + str(gamma[gencount]) + '\n')
        else:
            lingens[gencount] = gen.nodeID
      
    ampl.getSet('gens').setValues(list(gens))
    ampl.getSet('quadgens').setValues(list(quadgens))
    ampl.getSet('lingens').setValues(list(lingens))
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)
    ampl.get_parameter('alfa').setValues(alfa)
    ampl.get_parameter('beta').setValues(beta)
    ampl.get_parameter('gamma').setValues(gamma)
    ampl.get_parameter('sigma').set(sigma)

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
        
    # Solve
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()

    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")


    # Get solution

    total_costvar           = ampl.get_objective("total_cost")
    v2var                   = ampl.get_variable("v")
    cvar                    = ampl.get_variable("c")
    svar                    = ampl.get_variable("s")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['v2values']    = v2var.get_values().to_dict()
    all_data['cvalues']     = cvar.get_values().to_dict()
    all_data['svalues']     = svar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    
    PLoss   = {}
    QLoss   = {}
    for branch in all_data['Pfvalues'].keys():
        PLoss[branch] = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch] 
        QLoss[branch] = all_data['Qfvalues'][branch] + all_data['Qtvalues'][branch] 

    timesofar = time.time() - all_data['T0']        

    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" solver " + all_data['solver'] + "\n")
    log.joint(" objective " + str(all_data['objvalue']) + '\n')
    log.joint(" active power generation " + str(sum(all_data['GenPvalues'])) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(PLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(all_data['GenQvalues'])) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(QLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')    

    log.joint(" writing casename, modfile, obj and runtime to summary_socp.log\n")

    summary_socp = open("summary_socp.log","a+")

    summary_socp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile']
                       + ' solver ' + all_data['solver'] +  ' obj ' + str(all_data['objvalue'])
                       + ' runtime ' + str(timesofar) + '\n')

    summary_socp.close()


    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")

    if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)


def gosocp2_mosek(log,all_data):

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
    ampl.eval("option mosek_options 'tech:optionnativeread=/opt/newampl/mosekopt outlev=1 sol:chk:feastol=1e-08 chk:mode=2';")    
    #cvt:socp=0  writemodel=i2.ptf

    if all_data['fix_point']:
        getsol_knitro(log,all_data)
        tolerance = 1e-05
        log.joint(' fixing tolerance to ' + str(tolerance) + '\n')
        mp_vm    = all_data['mp_vm']
	#mp_angle = all_data['mp_angle']

    IDtoCountmap = all_data['IDtoCountmap']
        
    # Setting up buses
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

        bus_gens[buscount] = bus.genidsbycount 

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
    CSmax     = {}
    c_ubound  = {}
    c_lbound  = {}
    s_ubound  = {}
    s_lbound  = {}
    Cinit     = {}
    Sinit     = {}

    
    i2max       = {}
    g           = {}
    b           = {}
    bshunt      = {}
    ratio       = {}
    phase_angle = {}
    
    counter = 0 
    for branch in all_data['branches'].values():
        branchcount = branch.count
        if branch.status != 1: #the reader does not consider branches whose status is 0, hence this should never be true
            log.joint(' branch ' + str(branchcount) + ' off, we skip it\n')
            #breakexit('check')
            #continue
        
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

            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)
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
        if maxanglerad <= math.pi/2:
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
    
    ampl.get_parameter('i2max').setValues(i2max)    
    ampl.get_parameter('g').setValues(g) 
    ampl.get_parameter('b').setValues(b)
    ampl.get_parameter('bshunt').setValues(bshunt)    
    ampl.get_parameter('ratio').setValues(ratio)
    ampl.get_parameter('phase_angle').setValues(phase_angle)

    # Setting up generators
    gens      = {}
    lingens   = {}
    quadgens  = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} 
    alfa      = {}
    beta      = {}
    gamma     = {}
    sigma     = 0

    for gen in all_data['gens'].values():
        gencount = gen.count 
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


        fixedcost[gencount] = gen.costvector[2] * gen.status
        lincost[gencount] = gen.costvector[1]
        quadcost[gencount] = gen.costvector[0]
        
        if gen.costvector[0]:
            quadgens[gencount] = gen.nodeID
            alfa[gencount]     = math.sqrt(gen.costvector[0])            
            beta[gencount]     = - lincost[gencount] / (2 * alfa[gencount])
            gamma[gencount]    = fixedcost[gencount] - beta[gencount]**2
            sigma             += gamma[gencount]
            #log.joint(' alfa ' + str(alfa[gencount]) + ' beta ' + str(beta[gencount]) + ' gamma ' + str(gamma[gencount]) + '\n')
        else:
            lingens[gencount] = gen.nodeID

                
    ampl.getSet('gens').setValues(list(gens))
    ampl.getSet('quadgens').setValues(list(quadgens))
    ampl.getSet('lingens').setValues(list(lingens))
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)
    ampl.get_parameter('alfa').setValues(alfa)
    ampl.get_parameter('beta').setValues(beta)
    ampl.get_parameter('gamma').setValues(gamma)
    ampl.get_parameter('sigma').set(sigma)


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
        
    # Solver
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()


    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")


    # Get solution

    total_costvar           = ampl.get_objective("total_cost")
    v2var                   = ampl.get_variable("v")
    i2fvar                  = ampl.get_variable("i2f")
    cvar                    = ampl.get_variable("c")
    svar                    = ampl.get_variable("s")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['v2values']    = v2var.get_values().to_dict()
    all_data['i2fvalues']   = i2fvar.get_values().to_dict()
    all_data['cvalues']     = cvar.get_values().to_dict()
    all_data['svalues']     = svar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    PLoss   = {}
    QLoss   = {}
    for branch in all_data['Pfvalues'].keys():
        PLoss[branch] = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch] 
        QLoss[branch] = all_data['Qfvalues'][branch] + all_data['Qtvalues'][branch] 

    timesofar = time.time() - all_data['T0']        
    
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" solver " + all_data['solver'] + "\n")
    log.joint(" objective " + str(all_data['objvalue']) + '\n')
    log.joint(" active power generation " + str(sum(all_data['GenPvalues'])) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(PLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(all_data['GenQvalues'])) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(QLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')    

    log.joint(" writing casename, modfile, obj and runtime to summary_socp.log\n")

    summary_socp = open("summary_socp.log","a+")

    summary_socp.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile'] +  ' obj ' + str(all_data['objvalue']) + ' runtime ' + str(timesofar) + '\n')

    summary_socp.close()

    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")


    if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)
    
    
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
    
