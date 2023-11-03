from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
from cut_4 import cutheuristic_4
from supernodes_4 import supernodes_4
from cutjabr_4 import cutjabr_4
from gurobipy import * 
#import amplpy_gurobi as ampls_gurobi

#import igraph as ig
#import matplotlib.pyplot as plt


def goampl_4(log,all_data):

    ampl = AMPL()
    log.joint(' ampl object created \n')
    
    modfile = all_data['modfile']
    solver = all_data['solver']

    log.joint(' reading modfile ...\n')
    
    t0 = time.time()
    ampl.read(modfile)
    t1 = time.time()
    log.joint(' modfile read in time = {0}\n',t1-t0)
    print(" modfile read in time = {0}\n".format(t1-t0))
    
    ampl.eval("option display_precision 0;")
    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)    
    #ampl.setOption('presolve',0)
    log.joint(" solver set to = {:s} \n",solver)
    print(" solver set to = {0} \n".format(solver))

    log.joint(" ampl object created, modfile read, and solver chosen\n")
    
    ##### cuts ######
    

    #buses
    buses = {}
    Pd = {}
    Qd = {}
    Vmax = {}
    Vmin = {}    
    branches_f = {}
    branches_t = {}
    bus_gens = {}
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

        for branchid in bus.frombranchids.values():
            branches_f[buscount].append(branchid) 
            
        for branchid in bus.tobranchids.values():
            branches_t[buscount].append(branchid)      

    
    ampl.getSet('buses').setValues(list(buses))
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
    Rft = {}
    bus_f = {}
    bus_t = {}
    CSmax = {}
    zero_r = []

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
        Rft[branchcount] = branch.r

        #if branchcount == 1:

        CSmax[branchcount] = Vmax[branch.id_f] * Vmax[branch.id_t]
        bus_f[branchcount] = branch.id_f
        bus_t[branchcount] = branch.id_t
        
        if branch.r <= 0:
            zero_r.append(branchcount)
        
    #print("branch 8792",U[8792]) #case6468rte.m
    #print("buses",buses)
    #print("branches",branches)
    #breakexit("check branch")

    ampl.getSet('branches').setValues(list(branches))
    ampl.get_parameter('Gtt').setValues(Gtt)
    ampl.get_parameter('Btt').setValues(Btt)
    ampl.get_parameter('Gff').setValues(Gff)
    ampl.get_parameter('Bff').setValues(Bff)
    ampl.get_parameter('Gtf').setValues(Gtf)
    ampl.get_parameter('Btf').setValues(Btf)
    ampl.get_parameter('Gft').setValues(Gft)
    ampl.get_parameter('Bft').setValues(Bft) 
    ampl.get_parameter('Rft').setValues(Rft)
    ampl.get_parameter('U').setValues(U)
    ampl.get_parameter('bus_f').setValues(bus_f)
    ampl.get_parameter('bus_t').setValues(bus_t)
    ampl.get_parameter('CSmax').setValues(CSmax)


    #gens
    PMAX = [1.947467335079817,0.3676754551741309,0.273956066366522,0,0.09268611915958354]
    #PMAX = [-1,0,0,-1,-1]
    PMIN = [1.947467335079817,0.3676754551741309,0.273956066366522,0,0.09268611915958354]
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
        print("gencount =",gencount)
        gens[gencount] = gen.nodeID
        Pmax[gencount] = PMAX[gencount - 1]
        Pmin[gencount] = PMIN[gencount - 1]
        Qmax[gencount] = gen.Qmax
        Qmin[gencount] = gen.Qmin
        #print(" cost vector",gen.costvector) #check small case5.m
        #print(" cost degree",gen.costdegree)
        fixedcost[gencount] = gen.costvector[2]
        lincost[gencount] = gen.costvector[1]
        quadcost[gencount] = gen.costvector[0]
    
    print("Pmax",Pmax)
    ampl.getSet('gens').setValues(list(gens))
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)

    #CUTS (initialized as an empty list in main)
    ampl.getSet('cuts').setValues(all_data['cuts'])
    ampl.get_parameter('num_jabr_cuts').set(all_data['num_jabr_cuts'])
    ampl.getSet('jabr_cuts').setValues(all_data['jabr_cuts'])
    ampl.get_parameter('t_jabr').setValues(all_data['t_jabr'])
    ampl.get_parameter('c_jabr').setValues(all_data['c_jabr'])
    ampl.get_parameter('s_jabr').setValues(all_data['s_jabr'])
    ampl.get_parameter('vk_vm_jabr').setValues(all_data['vk_vm_jabr'])

    log.joint(" sets and parameters loaded\n")
    
    breakexit("done loading data")

    log.joint(" saving processed data to all_data\n")

    all_data['goampl_num_buses'] = len(buses)
    all_data['goampl_buses'] = buses
    all_data['goampl_branches'] = branches
    all_data['goampl_bus_f'] = bus_f
    all_data['goampl_bus_t'] = bus_t
    all_data['goampl_branches_f'] = branches_f
    all_data['goampl_branches_t'] = branches_t
    all_data['goampl_Pd'] = Pd
    all_data['goampl_gens'] = gens
    all_data['goampl_zero_r'] = zero_r

    #m = ampl.exportGurobiModel(gurobiDriver='gurobi')
    
    #m.write("ampltogurobi.lp")
    
    ampl.eval("option gurobi_options 'iisfind=1';")

    #EXPAND MODEL
    #ampl.eval("expand;")

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

    max_iterations = 50
    iteration = 0

    
    while iteration <= max_iterations:
       
        iteration += 1
        
        ampl.eval("expand;")

        #SOLVE
        log.joint(" solving model ...\n")
        t0 = time.time()
        ampl.solve()
        t1 = time.time()
        
        #GET SOLUTION
        total_cost = ampl.get_objective("total_cost")
        v = ampl.get_variable("v")
        c = ampl.get_variable("c")
        s = ampl.get_variable("s")
        Pf = ampl.get_variable("Pf")
        Pt = ampl.get_variable("Pt")
        Pg = ampl.get_variable("Pg")
        Loss = ampl.get_variable("PLoss")
        I = ampl.get_variable("I")
        dic_v = v.get_values().to_dict()
        dic_c = c.get_values().to_dict()
        dic_s = s.get_values().to_dict()
        dicPf = Pf.get_values().to_dict()
        dicPt = Pt.get_values().to_dict()
        dicPg = Pg.get_values().to_dict()
        dicLoss = Loss.get_values().to_dict()
        for loss in dicLoss:
            if dicLoss[loss] <= 1e-06:
                dicLoss[loss] = 0
        
        #saving data
        all_data['dic_v'] = dic_v
        all_data['dic_c'] = dic_c
        all_data['dic_s'] = dic_s
        all_data['dicPf'] = dicPf
        all_data['dicPt'] = dicPt
        all_data['dicPg'] = dicPg
        all_data['PLoss'] = dicLoss 

        #print(" voltages\n",v.get_values())
        print(" active power generation\n",Pg.get_values())
        print(" active power_f\n",Pf.get_values())
        print(" active power_t\n",Pt.get_values())
        print(" v\n",v.get_values())
        #print("type of v",type(v))
        #print(" v.lb\n",v.lb())
        #print(" v.ub\n",v.ub())
        print(" c\n",c.get_values())
        print(" s\n",s.get_values())
        print(" active power loss\n",Loss.get_values())
        print(" squared-current\n",I.get_values())
        
        print("non-zero active power loss\n")
        for branch in dicLoss.keys():
            if dicLoss[branch] > 0:
                print("branch id = {0}, branch = {1}, loss = {2}\n".format(branch,branches[branch],dicLoss[branch]))
        print(" ------------------------ ")
        
        print(" objective is:",total_cost.get().value())
        print(" total active power generation",sum(dicPg.values()))
        print(" total active power demand",sum(Pd.values()))
        print(" total active power loss",sum(dicLoss.values()))
        print(" solver runtime:",t1-t0)

        breakexit("solved")

        print(" ------------------------\n")
        
        #supernodes_4(log,all_data)   

        print(" ------------------------\n")
        
        #cutheuristic_4(log,all_data)

        print(" ------------------------\n")
        
        cutjabr_4(log,all_data)
        
        print(" ------------------------\n")

        
        log.joint(" updating AMPL model with new cuts ...\n")
        #ampl.getSet('cuts').setValues(all_data['cuts'])
        ampl.get_parameter('num_jabr_cuts').set(all_data['num_jabr_cuts'])
        ampl.getSet('jabr_cuts').setValues(all_data['jabr_cuts'])
        ampl.get_parameter('t_jabr').setValues(all_data['t_jabr'])
        ampl.get_parameter('c_jabr').setValues(all_data['c_jabr'])
        ampl.get_parameter('s_jabr').setValues(all_data['s_jabr'])
        ampl.get_parameter('vk_vm_jabr').setValues(all_data['vk_vm_jabr'])
        
        #constrs = ampl.get_constraint('j_cuts')
        #print("jabr_cuts constraints:")
        #for constr in constrs:
        #    print("contraint_name =",constr[1].name())
            
        print("iteration =",iteration)
        breakexit("iteration done")

        print(" ------------------------\n")
        
    
