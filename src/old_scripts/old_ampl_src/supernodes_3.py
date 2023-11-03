import numpy as np
from myutils import breakexit
from log import danoLogger

def supernodes_3(log,all_data):

    log.joint(" computing the super source and super target ...\n")

    gens = all_data['goampl_gens'].copy()
    Pd = all_data['goampl_Pd'].copy()
    dicPf = all_data['dicPf'].copy()
    dicPt = all_data['dicPt'].copy()
    dicPg = all_data['dicPg'].copy()

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

    #print("netgen dictionary",netgen)

    N = all_data['size_supernodes']
    print("new cuts?",all_data['new_cut'])

    breakexit("check")

    if all_data['new_cut'] == 0 and N < all_data['goampl_num_buses']:
        N += 1
        all_data['size_supernodes'] = N 
        print("N =",N)
    all_data['new_cut'] = 0


    set_s = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[:N])
    set_t = dict(sorted(netgen.items(), key = lambda x: x[1], reverse = True)[-N:])

    all_data['source'] = list(set_s.keys()) #before s                                       
    all_data['target'] = list(set_t.keys()) #before t 

    print("super source",list(set_s.keys()))
    print("super target",list(set_t.keys()))
    breakexit("supernodes function ok")
