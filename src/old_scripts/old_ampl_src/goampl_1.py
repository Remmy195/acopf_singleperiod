from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
from cut_1 import cutheuristic_1

import igraph as ig
import matplotlib.pyplot as plt


def goampl_1(log,all_data):

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
    ampl.setOption('presolve',0)
    log.joint(" solver set to = {:s} \n",solver)
    print(" solver set to = {0} \n".format(solver))

    log.joint(" ampl object created, modfile read, and solver chosen\n")
    
    ##### cuts ######
    ampl.getSet("cuts").setValues(all_data['cuts'])

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
    reverse_branches = {}
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


    for branch in all_data['branches'].values():
        branchcount = branch.count
        branches[branchcount] = (branch.id_f,branch.id_t)
        reverse_branches[(branch.id_f,branch.id_t)] = branchcount #THIS IS NOT CORRECT since there exist parallel lines, planning on using it with igraph
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
        CSmax[branchcount] = Vmax[branch.id_f] * Vmax[branch.id_t]
        bus_f[branchcount] = branch.id_f
        bus_t[branchcount] = branch.id_t

    #print("branch 8792",U[8792]) #case6468rte.m
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
    
    breakexit("done loading data")

    log.joint(" saving processed data to all_data\n")

    all_data['goampl_buses'] = buses
    all_data['goampl_gens'] = gens
    all_data['goampl_branches'] = branches
    all_data['goampl_bus_gens'] = bus_gens
    all_data['goampl_bus_f'] = bus_f
    all_data['goampl_bus_t'] = bus_t
    all_data['goampl_Pd'] = Pd
    all_data['goampl_U'] = U    
    all_data['goampl_branches_f'] = branches_f
    all_data['goampl_branches_t'] = branches_t

    #EXPAND MODEL
    ampl.eval("expand;")
    
    #SOLVE
    t0 = time.time()
    ampl.solve()
    t1 = time.time()
    breakexit("AMPL model expand")

    #GET SOLUTION

    total_cost = ampl.get_objective("total_cost")
    v = ampl.get_variable("v")
    current = ampl.get_variable("I")
    Pf = ampl.get_variable("Pf")
    Pt = ampl.get_variable("Pt")
    Pg = ampl.get_variable("Pg")
    Loss = ampl.get_variable("PLoss")
    dicPf = Pf.get_values().to_dict()
    dicPt = Pt.get_values().to_dict()
    dicPg = Pg.get_values().to_dict()
    dicLoss = Loss.get_values().to_dict()
    for loss in dicLoss:
        if dicLoss[loss] <= 1e-06:
            dicLoss[loss] = 0
    #print("voltages \n",v.get_values())
    print(" voltages\n",v.get_values())
    print(" sq-current\n",current.get_values())
    print(" active power generation\n",Pg.get_values())
    print(" active power_f\n",Pf.get_values())
    print(" active power_t\n",Pt.get_values())
    print(" active power loss\n",Loss.get_values())
    #print("type=",type(Loss.get_values()))
    print(" ------------------------ ")
    print(" objective is:",total_cost.get().value())
    print(" total active power generation",100*sum(dicPg.values()))
    print(" total active power demand",100*sum(Pd.values()))
    print(" total active power loss",100*sum(dicLoss.values()))
    print(" solver runtime:",t1-t0)

    breakexit("solved")
    ########### defining source and target for cut heuristic ########
    
    total_bus_gen = {}
    for gen in gens:
        if gens[gen] in total_bus_gen.keys():
            total_bus_gen[gens[gen]] += dicPg[gen]
        else:
            total_bus_gen[gens[gen]] = dicPg[gen]

    netgen = {}
    for bus in Pd.keys():
        if bus in total_bus_gen.keys():
            netgen[bus] = total_bus_gen[bus] - Pd[bus]
        else:
            netgen[bus] = - Pd[bus]

    print("netgen dictionary",netgen)
    
    #N = 1
    #set_s = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[:N])
    #set_t = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[-N:])

    #print("set_s",set_s)
    #print("set_t",set_t)
    #print(set_s.keys())
    breakexit("check dictionary")

    s = max(netgen, key = netgen.get)
    t = min(netgen, key = netgen.get)

    print("source = {0}, target = {1}".format(s,t))
    breakexit(" source and target computed")

    ###################### saving data #########################

    all_data['goampl_Pg'] = dicPg   
    all_data['goampl_PLoss'] = dicLoss
    all_data['goampl_source'] = s
    all_data['goampl_target'] = t
    all_data['goampl_delta_s'] = branches_f[s]
    
    breakexit("data saved")
    
    print(" ------------------------\n")

    #print(" setting up the undirected graph ")
    
    #edges = []
    #for edge in branches.values():
    #    edges.append((edge[0]-1,edge[1]-1))
    #edges = list(branches.values())
    
    #print("the edges are:",edges)

    #for edge in edges:
    #    print(edge - 1)
    
    #epsilon = 0.01

    #capacities = {}
    #for i in dicLoss.keys():
        #capacities[i] = (dicLoss[i] + epsilon ) * (1/max(abs(dicPf[i]),abs(dicPt[i])))
    #    capacities[i] = dicLoss[i]
    #print(capacities)
    #g = ig.Graph(number_buses,edges,directed=False)
    #g.es["capacity"] = list(capacities.values())

    #cut = g.mincut(capacity = list(capacities.values()))

    #print("min cut",cut.value)
    
    #L1 = []
    #L2 = []
    #for bus in cut.partition[0]:
    #    L1.append(bus+1)

    #for bus in cut.partition[1]:
    #    L2.append(bus+1)
    #print("buses partition given by the cut:",L1,L2)

    #ct = []
    #for edgeid in cut.cut:
    #    e = g.es[edgeid]
    #    s = e.source + 1
    #    t = e.target + 1
    #    ct.append((s,t))

    #print("cut:",ct)
    #print(" printing branches, branch id's, and flow:")
    #for branch in ct:
    #    b = reverse_branches[branch]
    #    print("branch,branchid,power_f,power_t= ",branch,b,dicPf[b],dicPt[b])

    #fig, ax = plt.subplots()
    #ig.plot(
    #g,
    #target=ax,
    #layout="kk",
    #vertex_label=range(1,g.vcount()+1), #changed the labels
    #vertex_color="lightblue",
    #vertex_size = 0.3
    #)
    #plt.show()
