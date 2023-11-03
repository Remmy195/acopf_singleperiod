import numpy as np
from myutils import breakexit
from log import danoLogger
import copy
import math

def supernodes_5(log,all_data):

    log.joint(" computing the super source and super target ...\n")

    gens = copy.deepcopy(all_data['goampl_gens'])
    Pd = copy.deepcopy(all_data['goampl_Pd'])
    dicPf = copy.deepcopy(all_data['dicPf'])
    dicPt = copy.deepcopy(all_data['dicPt'])
    dicPg = copy.deepcopy(all_data['dicPg'])

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

    N = all_data['size_supernodes']
    print("new cuts?",all_data['new_mincut'])

    if all_data['new_mincut'] == 0 and N < math.floor(all_data['goampl_num_buses']/2) - 1:
        N += 1
        all_data['size_supernodes'] = N 
        print("N (size of supernodes) =",N)
        
    all_data['new_mincut'] = 0

    set_s = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[:N])
    set_t = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[-N:])

    all_data['source'] = list(set_s.keys()) 
    all_data['target'] = list(set_t.keys()) 

    print("super source",list(set_s.keys()))
    print("super target",list(set_t.keys()))
    #breakexit("supernodes function ok")
