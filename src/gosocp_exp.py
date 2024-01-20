from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math
import random


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
        ampl.eval("option gurobi_options 'method=2 barhomogeneous=1 numericfocus=1 barconvtol=1e-6 outlev=1 iisfind=1 writemodel=jabr.lp resultfile=jabr.ilp';")

    if all_data['solver'] == 'knitroampl':
        if all_data['mytol']:
            ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';")
        elif all_data['multistart']:
            ampl.eval("option knitro_options 'ms_enable=1 ms_numthreads=10 ms_maxsolves=5 ms_terminate =1';")
        elif all_data['knitropresolveoff']:
            ampl.eval("option knitro_options 'presolve=0';")       
        else:
            ampl.eval("option knitro_options 'blasoptionlib=1 numthreads=20 linsolver=7 maxtime_real=1000';") 

            
    #bar_conic_enable=1 --> too slow for some problems ...
    #bar_refinement=1
            
    #if all_data['solver'] == 'knitroampl' and all_data['conic']:
    #    ampl.eval("option knitro_options 'bar_conic_enable=1';")


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
    branches = {}
    U        = {}
    Gtt      = {}
    Btt      = {}
    Gff      = {}
    Bff      = {}
    Gtf      = {}
    Btf      = {}
    Gft      = {}
    Bft      = {}
    bus_f    = {}
    bus_t    = {}
    CSmax    = {}
    Cinit    = {}
    Sinit    = {}
    slackjabr = {}
    
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

        slackjabr[branchcount] = 0

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

    ampl.get_parameter('slackjabr').setValues(slackjabr)
    
    #gens
    gens = {}
    Pmax = {}
    Pmin = {}
    Qmax = {}
    Qmin = {}
    fixedcost = {}
    lincost = {}
    quadcost = {} #we assume at most quadratic cost (for the moment)

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
    

    i = 1
    all_data['losstest_dic']        = {}
    all_data['losstest_loss']       = {}
    all_data['losstest_branchloss'] = {}
    all_data['losstest_obj']        = {}
    all_data['losstest_Pg']         = {}
    numbranches = all_data['numbranches']


    Yes = True

    
    while i <= 100:
        
        randombranch = random.randrange(1,numbranches+1)
        branch       = all_data['branches'][randombranch]
        
        while branch.status == 0:
            randombranch = random.randrange(1,numbranches+1)
            branch       = all_data['branches'][randombranch]

        f            = branch.f
        t            = branch.t
        Vmaxf        = Vmax[IDtoCountmap[f]]
        Vmaxt        = Vmax[IDtoCountmap[t]]
        slack        = 2 * ( (Vmaxf**2) * (Vmaxt**2) )
        #( Vmaxf * Vmaxf + Vmaxt * Vmaxt )
        letbranch    = "let slackjabr[" + str(randombranch) + "] := " + str(slack) + ";" 

        log.joint(' setting slack of random branch ...\n')
        log.joint(letbranch)

        if Yes:
            ampl.eval(letbranch)


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

        letbranch    = "let slackjabr[" + str(randombranch) + "] := " + str(0) + ";"
        log.joint('setting slack of random branch ...\n')

        log.joint(letbranch + '\n')
        ampl.eval(letbranch)
        
        
        log.joint(" ------------------------\n")


        #GET SOLUTION
        total_cost = ampl.get_objective("total_cost")
        Pf         = ampl.get_variable("Pf")
        Pt         = ampl.get_variable("Pt")
        Pg         = ampl.get_variable("Pg")
        dic_Pg     = Pg.get_values().to_dict()
        dic_Pf     = Pf.get_values().to_dict()
        dic_Pt     = Pt.get_values().to_dict()
        dicPLoss   = {}

        for branch in dic_Pf.keys():
            dicPLoss[branch] = dic_Pf[branch] + dic_Pt[branch]
            log.joint(' loss branch ' + str(branch) + ' = '
                      + str(dicPLoss[branch]) + '\n')

        objval      = total_cost.get().value()
        timesofar   = time.time() - all_data['T0']
        rbranchloss = dic_Pf[randombranch] + dic_Pt[randombranch]
        
        all_data['losstest_dic'][i] = (randombranch,
                                       objval,
                                       sum(dicPLoss.values()),
                                       rbranchloss)

        all_data['losstest_loss'][i]       = sum(dicPLoss.values())
        all_data['losstest_branchloss'][i] = rbranchloss
        all_data['losstest_obj'][i]        = objval
        all_data['losstest_Pg'][i]         = sum(dic_Pg.values())
        
        avg = sum(all_data['losstest_loss'].values()) / i        
        avg_branch = sum(all_data['losstest_branchloss'].values()) / i
        avg_obj    = sum(all_data['losstest_obj'].values()) / i
        avg_Pg     = sum(all_data['losstest_Pg'].values()) / i

        maxloss        = max(all_data['losstest_loss'].values())
        maxloss_branch = max(all_data['losstest_branchloss'].values())

        minloss        = min(all_data['losstest_loss'].values())
        minloss_branch = min(all_data['losstest_branchloss'].values())
        minPg          = min(all_data['losstest_Pg'].values())


        log.joint(' iteration ' + str(i) + '\n')
        log.joint(" case " + all_data['casefilename'] + "\n")
        log.joint(" modfile " + all_data['modfile'] + "\n")
        log.joint(" objective " + str(objval) + '\n')
        log.joint(" total active power generation "
                  + str(sum(dic_Pg.values())) + '\n')
        log.joint(" total active power demand " + str(sum(Pd.values())) + '\n')
        log.joint(" total active power loss "
                  + str(sum(dicPLoss.values())) + '\n')

        log.joint(' current Pg  avg ' + str(avg_Pg) + '\n' )
        log.joint(' current obj  avg ' + str(avg_obj) + '\n' )
        log.joint(' current loss avg ' + str(avg) + '\n' )
        log.joint(' current loss avg (branch)' + str(avg_branch) + '\n' )
        log.joint(' current max loss ' + str(maxloss) + '\n')
        log.joint(' current max loss (at the branch)' + str(maxloss_branch)
                  + '\n')
        log.joint(' current min loss ' + str(minloss) + '\n')
        log.joint(' current min loss (at the branch)' + str(minloss_branch)
                  + '\n')
        log.joint(' current min Pg ' + str(minPg) + '\n')

        log.joint(" time so far " + str(timesofar) + '\n')

        i += 1

        #log.joint('losstest loss ' + str(all_data['losstest_loss']) + '\n')
        #log.joint('losstest loss (branch)'
        #          + str(all_data['losstest_branchloss']) + '\n')
        
        #breakexit('c')
        log.joint(" ------------------------\n")
    
