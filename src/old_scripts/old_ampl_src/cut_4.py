from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time

def cutheuristic_4(log,all_data):

    log.joint(" ***********************\n")

    log.joint(" min-cut heuristic\n")

    log.joint(" ***********************\n")

    ampl = AMPL()

    modfile = all_data['modfile_cut']
    solver = all_data['solver_cut']

    #log.joint(" resetting ampl ... \n")

    ampl.eval("reset;")

    log.joint("\n reading modfile for max-flow ...\n")

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

    buses = all_data['goampl_buses'].copy()
    branches = all_data['goampl_branches'].copy()
    bus_f = all_data['goampl_bus_f'].copy() 
    bus_t = all_data['goampl_bus_t'].copy()
    branches_f = all_data['goampl_branches_f'].copy()
    branches_t = all_data['goampl_branches_t'].copy()
    zero_r = all_data['goampl_zero_r'].copy()
    PLoss = all_data['PLoss'].copy()
    s = all_data['source'].copy()
    t = all_data['target'].copy()

    ######### addding new source and new target ########

    total_loss = sum(PLoss.values())
    print("total loss = ",total_loss)
    extracapacity = 1
    print("capacities for auxiliary edges and 0-resistance branches = ",total_loss + extracapacity)
    breakexit("capacities")

    num_buses_original = len(buses)
    num_branches_original = len(branches)

    print("number of buses = {0}, number of branches = {1}".format(len(buses),len(branches)))
    #print("original branches",branches)
    print("super source",s)
    print("super target",t)
    
    #print("old bus_f",bus_f)
    #print("old bus_t",bus_t)

    super_source = num_buses_original + 1
    super_target = num_buses_original + 2

    buses[super_source] = 0
    buses[super_target] = -1

    #buses[super_source] = super_source 0
    #buses[super_target] = super_target 0
    
    branches_f[super_source] = []
    branches_t[super_source] = [] #will be empty
    branches_t[super_target] = []
    branches_f[super_target] = [] #will be empty

    #newbranch = 0 #old 0
    #for i in s:
    #    newbranch += 1
    #    branches[num_branches_original + newbranch] = (super_source,i)
    #    bus_f[num_branches_original + newbranch] = super_source
    #    bus_t[num_branches_original + newbranch] = i
    #    branches_f[super_source].append(num_branches_original + newbranch)
    #    if (num_branches_original + newbranch) not in branches_t[i]:    
    #        branches_t[i].append(num_branches_original + newbranch)
    #    print("added branch",(super_source,i))
    #    PLoss[num_branches_original + newbranch] = total_capacity + extracapacity
        

    #print("BBBBB",branches_t[1])

    #for i in t:
    #    newbranch += 1
    #    branches[num_branches_original + newbranch] = (i,super_target)
    #    bus_f[num_branches_original + newbranch] = i
    #    bus_t[num_branches_original + newbranch] = super_target
    #    if (num_branches_original + newbranch) not in branches_f[i]:    
    #        branches_f[i].append(num_branches_original + newbranch)
    #    branches_t[super_target].append(num_branches_original + newbranch)
    #    print("added branch",(i,super_target))
    #    PLoss[num_branches_original + newbranch] = total_capacity + extracapacity

    #assigning large capacity to branches with 0 or negative resistance (Polish cases)
    for branchid in PLoss:
        if branchid in zero_r:
            PLoss[branchid] = total_loss + extracapacity


    newbranch = 0
    for i in s:
        newbranch += 1
        branches[num_branches_original + newbranch] = (0,i)
        bus_f[num_branches_original + newbranch] = 0
        bus_t[num_branches_original + newbranch] = i
        branches_f[super_source].append(num_branches_original + newbranch)
        if (num_branches_original + newbranch) not in branches_t[i]:    
            branches_t[i].append(num_branches_original + newbranch)
        print("added branch",(0,i))
        PLoss[num_branches_original + newbranch] = total_loss + extracapacity
        

    for i in t:
        newbranch += 1
        branches[num_branches_original + newbranch] = (i,-1)
        bus_f[num_branches_original + newbranch] = i
        bus_t[num_branches_original + newbranch] = -1
        if (num_branches_original + newbranch) not in branches_f[i]:    
            branches_f[i].append(num_branches_original + newbranch)
        branches_t[super_target].append(num_branches_original + newbranch)
        print("added branch",(i,-1))
        PLoss[num_branches_original + newbranch] = total_loss + extracapacity

    print("new number of buses = {0}, new number of branches = {1}".format(len(buses),len(branches)))
    #breakexit("check")

    #print("buses",buses)
    #print("new branches",branches)
    #print("branches_f source",branches_f[super_source])
    #print("branches_t target",branches_t[super_target])
    #print("branches_f target",branches_f[super_target])
    #print("branches_t source",branches_t[super_source])
    
    #print("capacities",PLoss)
    #print("new bus_f",bus_f)
    #print("new bus_t",bus_t)
    breakexit("check")

    ####################################################

    #redefining delta_s, s and t
    delta_s = branches_f[super_source]
    delta_t = branches_t[super_target]
    s = super_source
    t = super_target
    
    print("delta_s",delta_s)
    print("delta_t",delta_t)
    print("buses",buses)
    print("branches",branches)
    print("branches_f[super_source]",branches_f[super_source])
    print("branches_t[super_source]",branches_t[super_target])
    print("PLoss",PLoss)
    #s = 0
    #t = -1

    ######### loading processed data to AMPL ###########
    
    ampl.getSet('buses').setValues(list(buses))
    ampl.getSet('branches').setValues(list(branches))
    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')


    for bus in all_data['buses'].values():
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
    
    branches_f_set[super_source].setValues(branches_f[super_source])
    branches_t_set[super_source].setValues(branches_t[super_source]) #before branches_f everywhere
    branches_f_set[super_target].setValues(branches_f[super_target])
    branches_t_set[super_target].setValues(branches_t[super_target])
    

    #for buscount in buses.values(): 0
    #    branches_f_set[buscount].setValues(branches_f[buscount]) 0
    #    branches_t_set[buscount].setValues(branches_t[buscount]) 0


    ampl.getSet('delta_s').setValues(delta_s)
    ampl.getSet('delta_t').setValues(delta_t)
    ampl.getSet('st').setValues([s,t])
    ampl.get_parameter('bus_f').setValues(bus_f)
    ampl.get_parameter('bus_t').setValues(bus_t)
    ampl.get_parameter('PLoss').set_values(PLoss)
    
    #ampl.eval("expand;")

    breakexit("AMPL model expanded")

    #SOLVE                                                                                                           
    t0 = time.time()
    ampl.solve()
    t1 = time.time()
    breakexit("check Pf[24]")
    #GET SOLUTION                                                                                                    
    max_flow = ampl.get_objective("flow")
    print("objective is:",max_flow.get().value())
    print("solver runtime:",t1-t0)
    print("super source",all_data['source'])
    print("super target",all_data['target'])
    constrs_f = ampl.get_constraint("Cap_f")
    for constr in constrs_f:
        if constr[1].dual() > 1e-6:
            added = 0
            constrname = constr[1].name()
            constrID = int(constrname.split('[')[1].split(']')[0])
            if constrID not in all_data['cuts'] and (constrID not in delta_s and constrID not in delta_t):
                all_data['new_cut'] = 1
                all_data['cuts'].append(int(constrID))
                added = 1
            if (constrID in delta_s) or (constrID in delta_t):
                breakexit("restart super-source and super-target sets!")
                all_data['size_supernodes'] = 0
                branches_f = all_data['goampl_branches_f'].copy()
                branches_t = all_data['goampl_branches_t'].copy()
                print("branches_f",branches_f)
                return True
            print("constrname = {0}, (id_f,id_t) = {1}, dual = {2}, new = {3}".format(constrname,branches[constrID],constr[1].dual(),added))

    constrs_t = ampl.get_constraint("Cap_t")
    for constr in constrs_t:
        if constr[1].dual() > 1e-6:
            added = 0
            constrname = constr[1].name()
            constrID = int(constrname.split('[')[1].split(']')[0])
            if constrID not in all_data['cuts'] and (constrID not in delta_s and constrID not in delta_t):
                all_data['new_cut'] = 0
                all_data['cuts'].append(int(constrID))
                added = 1
            if (constrID in delta_s) or (constrID in delta_t):
                breakexit("restart super-source and super-target sets!")
                all_data['size_supernodes'] = 0
                branches_f = all_data['goampl_branches_f'].copy()
                branches_t = all_data['goampl_branches_t'].copy()
                return True 
            print("constrname = {0}, (id_f,id_t) = {1}, dual = {2}, new = {3}".format(constrname,branches[constrID],constr[1].dual(),added))
            
    print("current set of cuts (by branchids) =",all_data['cuts'])
    print("current number of cuts (min-cut heuristic) =",len(all_data['cuts']))

    #ampl.eval("display _conname, _con.dual;") 

    branches_f = all_data['goampl_branches_f'].copy()
    branches_t = all_data['goampl_branches_t'].copy()
    
    breakexit("min-cut heuristic ok")
