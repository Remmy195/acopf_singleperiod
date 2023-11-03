from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time



def cutheuristic_2(log,all_data):
    ampl = AMPL()

    modfile = all_data['modfile_cut']
    solver = all_data['solver_cut']

    log.joint(" resetting ampl ... \n")

    ampl.eval("reset;")

    log.joint(" reading modfile ...\n")

    t0 = time.time()
    ampl.read(modfile)
    t1 = time.time()
    log.joint(' modfile read in time = {0}\n',t1-t0)
    print(" modfile read in time = {0}\n".format(t1-t0))

    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)
    log.joint(" solver set to = {:s} \n",solver)
    print(" solver set to = {0} \n".format(solver))

    log.joint(" ampl object created, modfile read, and solver chosen\n")

    log.joint(" loading data from all_data\n")

    buses = all_data['goampl_buses']
    gens = all_data['goampl_gens'] 
    branches = all_data['goampl_branches']
    bus_gens = all_data['goampl_bus_gens']
    bus_f = all_data['goampl_bus_f'] 
    bus_t = all_data['goampl_bus_t'] 
    branches_f = all_data['goampl_branches_f']
    branches_t = all_data['goampl_branches_t'] 
    Pd = all_data['goampl_Pd'] 
    Pg = all_data['goampl_Pg']
    PLoss = all_data['goampl_PLoss']
    s = all_data['goampl_source']
    t = all_data['goampl_target']
    #delta_s = all_data['goampl_delta_s']

    ######### addding new source and new target ########

    total_capacity = sum(PLoss.values())
    extracapacity = 1 #if losses are exactly 0

    num_buses_original = len(buses)
    num_branches_original = len(branches)

    print("number of buses, number of branches",len(buses),len(branches))
    print("original branches",branches)
    print("s",s)
    print("t",t)
    
    print("old bus_f",bus_f)
    print("old bus_t",bus_t)

    super_source = num_buses_original + 1
    super_target = num_buses_original + 2

    buses[super_source] = super_source
    buses[super_target] = super_target
    
    branches_f[super_source] = []
    branches_t[super_source] = [] #will be empty
    branches_t[super_target] = []
    branches_f[super_target] = [] #will be empty

    newbranch = 0
    for i in s:
        newbranch += 1
        branches[num_branches_original + newbranch] = (super_source,i)
        bus_f[num_branches_original + newbranch] = super_source
        bus_t[num_branches_original + newbranch] = i
        branches_f[super_source].append(num_branches_original + newbranch)
        branches_t[i].append(num_branches_original + newbranch)
        print("added branch",(super_source,i))
        PLoss[num_branches_original + newbranch] = total_capacity + extracapacity
        

    for i in t:
        newbranch += 1
        branches[num_branches_original + newbranch] = (i,super_target)
        bus_f[num_branches_original + newbranch] = i
        bus_t[num_branches_original + newbranch] = super_target
        branches_f[i].append(num_branches_original + newbranch)
        branches_t[super_target].append(num_branches_original + newbranch)
        print("added branch",(i,super_target))
        PLoss[num_branches_original + newbranch] = total_capacity + extracapacity

    print("number of buses, number of branches",len(buses),len(branches))
    breakexit("check")
    
    print("new branches",branches)
    print("branches_f source",branches_f[super_source])
    print("branches_t target",branches_t[super_target])
    print("capacities",PLoss)
    print("new bus_f",bus_f)
    print("new bus_t",bus_t)
    breakexit("check")


    ####################################################

    #redefining delta_s, s and t
    delta_s = branches_f[super_source]
    s = super_source
    t = super_target

    ######### loading processed data to AMPL ###########
    
    ampl.getSet('buses').setValues(list(buses))
    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    #bus_gens_set = ampl.getSet('bus_gens')

    for buscount in buses.values():
        branches_f_set[buscount].setValues(branches_f[buscount])
        branches_t_set[buscount].setValues(branches_t[buscount])
        #bus_gens_set[bus.count].setValues(bus_gens[bus.count])

    ampl.getSet('branches').setValues(list(branches))
    #ampl.getSet('gens').setValues(list(gens))
    ampl.getSet('delta_s').setValues(delta_s)
    ampl.getSet('st').setValues([s,t])
    ampl.get_parameter('bus_f').setValues(bus_f)
    ampl.get_parameter('bus_t').setValues(bus_t)
    #ampl.get_parameter('Pg').set_values(Pg)
    #ampl.get_parameter('Pd').set_values(Pd)
    ampl.get_parameter('PLoss').set_values(PLoss)
    
    ampl.eval("expand;")

    breakexit("AMPL model expanded")

    #SOLVE                                                                                                           
    t0 = time.time()
    ampl.solve()
    t1 = time.time()
    
    #GET SOLUTION                                                                                                    
    max_flow = ampl.get_objective("flow")
    print("objective is:",max_flow.get().value())
    print("solver runtime:",t1-t0)
    constrs = ampl.get_constraint("Cap_f")
    for constr in constrs:
        if constr[1].dual() > 1e-6:
            constrname = constr[1].name()
            constrID = constrname.split('[')[1].split(']')[0]
            if int(constrID) not in all_data['cuts']:
                all_data['cuts'].append(int(constrID))
            print("constrname",constrname)
            print("dual",constr[1].dual())

    print("new cuts associated to branches =",all_data['cuts'])

    breakexit("AMPL model solved and set of cuts updated")

    ampl.eval("display Pf, Pt,  _conname, _con.lb , _con.ub , _con.dual;") 
