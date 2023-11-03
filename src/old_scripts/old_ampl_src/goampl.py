from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math


def goampl(log,all_data):
    ampl = AMPL()

    all_data['ampl_object'] = ampl
    
    log.joint(' ampl object created \n')
    
    modfile = all_data['modfile']
    solver = all_data['solver']

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

    if all_data['solver'] == 'gurobi_ampl':
        ampl.eval("options option gurobi_options 'method=2';")
        ampl.eval("options option gurobi_options 'bar=1';")

    log.joint(" ampl object created, modfile read, and solver chosen\n")
    
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

    number_buses = 0

    for bus in all_data['buses'].values():
        number_buses += 1
        buscount = bus.count
        buses[buscount] = bus.nodeID 
        Pd[buscount] = bus.Pd
        Qd[buscount] = bus.Qd
        Vmax[buscount] = bus.Vmax
        Vmin[buscount] = bus.Vmin
        branches_f[buscount] = []
        branches_t[buscount] = []
        bus_gens[buscount] = bus.genidsbycount #load this as a sparse ds
        if bus.Gs != 0:
            bus_Gs[buscount] = bus.Gs
        if bus.Bs != 0:
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

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set = ampl.getSet('bus_gens')
    
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
        Pmax[gencount] = gen.Pmax
        Pmin[gencount] = gen.Pmin
        Qmax[gencount] = gen.Qmax
        Qmin[gencount] = gen.Qmin
        #print(" cost vector",gen.costvector) #check small case5.m
        #print(" cost degree",gen.costdegree)
        fixedcost[gencount] = gen.costvector[2]
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

        amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb,_var[j].ub);'     #shows full model                                 

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

    breakexit("now the solution")

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
    dic_v = v.get_values().to_dict()
    dic_c = c.get_values().to_dict()
    dic_s = s.get_values().to_dict()
    dicPg = Pg.get_values().to_dict()
    dicPf = Pf.get_values().to_dict()
    dicPt = Pt.get_values().to_dict()
    dicQf = Qf.get_values().to_dict()
    dicQt = Qt.get_values().to_dict()
    
    dicPLoss = {}
    dicQLoss = {}
    
    for branch in dicPf.keys():
        dicPLoss[branch] = dicPf[branch] + dicPt[branch] 
        dicQLoss[branch] = dicQf[branch] + dicQt[branch] 

    log.joint(" objective " + str(total_cost.get().value()) + '\n')
    log.joint(" total active power generation " + str(sum(dicPg.values())) + '\n')
    log.joint(" total active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" total active power loss " + str(sum(dicPLoss.values())) + '\n')
    log.joint(" total net reactive power loss " + str(sum(dicQLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(time.time()-all_data['timestart']) + '\n')
    
    log.joint(" ------------------------\n")
        
    
