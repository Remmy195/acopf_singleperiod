from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time

def goampl(log,all_data):

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
    

    ampl.setOption('solver',solver)    
    log.joint(" solver set to = {:s} \n",solver)
    print(" solver set to = {0} \n".format(solver))

    log.joint(" ampl object created, modfile read, and solver chosen\n")

    print("dict of branches",all_data['branches'])
    print("values of branches",all_data['branches'].values())
    breakexit("check branches")
    
    #buses
    buses = {}
    Pd = {}
    Qd = {}
    Vmax = {}
    Vmin = {}
    
    branches_f = {}
    branches_t = {}
    
    allbranches = all_data['branches']
    for bus in all_data['buses'].values():
        buses[bus.nodeID] = bus.count
        Pd[bus.nodeID] = bus.Pd
        Qd[bus.nodeID] = bus.Qd
        Vmax[bus.nodeID] = bus.Vmax
        Vmin[bus.nodeID] = bus.Vmin
        #branches_f[bus] = [] #list or dictionary?
        #branches_t[bus] = []
        branches_f[bus] = {}
        branches_t[bus] = {}
        
        #print(bus.frombranchids)
        #breakexit("check")
        #print("\n branches_f \n")
        for branchid in bus.frombranchids.values():
            branches_f[bus] = branchid 
            branches_f[bus].append((allbranches[branchid].f,allbranches[branchid].t))
            #allbranches[branchid].show(log)
        
        #print("bus id",bus.nodeID)
        #print("branches_f ",branches_f[bus])
        #print(" ---- \n")
        
        #print(" branches_t \n")
        for branchid in bus.tobranchids.values():
            branches_t[bus].append((allbranches[branchid].f,allbranches[branchid].t))
            #allbranches[branchid].show(log)
        
        #print("bus id",bus.nodeID)
        #print("branches_t ",branches_t[bus])
        #print("\n")
        #breakexit(" END BUS ")

    print("Vmax",Vmax)
    print("Vmin",Vmin)
    
    ampl.getSet('buses').setValues(list(buses))
    ampl.get_parameter('Pd').set_values(Pd)
    ampl.get_parameter('Qd').set_values(Qd)
    ampl.get_parameter('Vmax').set_values(Vmax)
    ampl.get_parameter('Vmin').set_values(Vmin)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    
    for bus in all_data['buses'].values():
        print("nodeID",bus.nodeID)
        branches_f_set[bus.nodeID].setValues(branches_f[bus])
        branches_t_set[bus.nodeID].setValues(branches_t[bus])
        
        #ampl.getSet(branches_set[bus]).setValues(branches_f[bus])

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

    for branch in all_data['branches'].values():
        branches[(branch.id_f,branch.id_t)] = branch.count
        U[(branch.id_f,branch.id_t)] = branch.limit
        Gtt[(branch.id_f,branch.id_t)] = branch.Gtt
        Btt[(branch.id_f,branch.id_t)] = branch.Btt
        Gff[(branch.id_f,branch.id_t)] = branch.Gff
        Bff[(branch.id_f,branch.id_t)] = branch.Bff
        Gtf[(branch.id_f,branch.id_t)] = branch.Gtf
        Btf[(branch.id_f,branch.id_t)] = branch.Btf
        Gft[(branch.id_f,branch.id_t)] = branch.Gft
        Bft[(branch.id_f,branch.id_t)] = branch.Bft

    print("Gff",Gff)
    print("Gtt",Gtt)
    print("Bff",Bff)
    print("Btf",Btf)
    print("Bft",Bft)
    print("Gft",Gft)
    print("Gtf",Gtf)
    print("U",U)

    ampl.getSet('branches').setValues(list(branches))
    ampl.get_parameter('Gtt').setValues(Gtt)
    ampl.get_parameter('Btt').setValues(Btt)
    ampl.get_parameter('Gff').setValues(Gff)
    ampl.get_parameter('Bff').setValues(Bff)
    ampl.get_parameter('Gtf').setValues(Gtf)
    ampl.get_parameter('Btf').setValues(Btf)
    ampl.get_parameter('Gft').setValues(Gft)
    ampl.get_parameter('Bft').setValues(Btf)
    ampl.get_parameter('U').setValues(U)

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
        gens[gen.nodeID] = gen.count #check this identifier... COUNT  
        Pmax[gen.nodeID] = gen.Pmax
        Pmin[gen.nodeID] = gen.Pmin
        Qmax[gen.nodeID] = gen.Qmax
        Qmin[gen.nodeID] = gen.Qmin
        #print(" cost vector",gen.costvector)
        #print(" cost degree",gen.costdegree)
        fixedcost[gen.nodeID] = gen.costvector[2]
        lincost[gen.nodeID] = gen.costvector[1]
        quadcost[gen.nodeID] = gen.costvector[0]

    print("Pmax",Pmax)
    print("Pmin",Pmin)
    print("Qmax",Qmax)
    print("Qmin",Qmin)
    
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

    #LOADING OBJECTIVE
    
    #SOLVE

    ampl.eval("expand;")
    
    ampl.solve()

    breakexit("solved?")
    #GET SOLUTION

    total_cost = ampl.get_objective("total_cost")

    print("Objective is:",total_cost.get().value())
    v = ampl.get_variable("v")
    Pf = ampl.get_variable("Pf")
    Pt = ampl.get_variable("Pt")
    print(v.get_values())
    print(Pf.get_values())
    print(Pt.get_values())
    #print(v.get_values())
