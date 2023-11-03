import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def cutjabr_4(log,all_data):
    
    log.joint(" ***********************\n")

    log.joint(" jabr-cuts heuristic\n")
    
    log.joint(" ***********************\n")
    
    log.joint("\n loading data from all_data\n")

    buses = all_data['goampl_buses'].copy()
    branches = all_data['goampl_branches'].copy()
    bus_f = all_data['goampl_bus_f'].copy() 
    bus_t = all_data['goampl_bus_t'].copy()
    
    v = all_data['dic_v']
    c = all_data['dic_c']
    s = all_data['dic_s']
    
    violation = {}
    violated = 0
    for i in branches:
        violation[i] = c[i] + s[i] - v[bus_f[i]] * v[bus_t[i]] 
        if violation[i] > 0:
            violated += 1
    
    print("number of violated jabr's =",violated)
    print("violations =",violation)

    log.joint(" checking for most violated jabrs ... \n")

    N = math.ceil(violated/3)
    
    most_violated = dict(sorted(violation.items(), key = lambda x: x[1], reverse = True)[:N])

    print(" we pick the {0}% most violated jabrs\n".format(100 * N/violated))
    print(" the {0} most violated cuts = {1}".format(N,most_violated))
    
    breakexit("check")

    log.joint(" computing new jabr cuts \n")

    num_cuts = all_data['num_jabr_cuts']
    new_cuts = []
    new_t_jabr = []
    new_c_jabr = []
    new_s_jabr = []
    new_vk_vm_jabr = []

    counter = 0
    for i in most_violated.keys():
        counter += 1
        new_cuts.append((num_cuts + counter,i))
        new_t_jabr.append(math.sqrt( (2 * c[i])**2 + (2 * s[i])**2 + (v[bus_f[i]] - v[bus_t[i]])**2) )
        new_c_jabr.append(2 * c[i])
        new_s_jabr.append(2 * s[i])
        new_vk_vm_jabr.append(v[bus_f[i]] - v[bus_t[i]])

    log.joint(" saving data to all_data ... \n")

    #print("before\n")
    #print("num",all_data['num_jabr_cuts'])
    #print("jabr_cuts",all_data['jabr_cuts'])    
    #print("t_jabr",all_data['t_jabr'])    

    all_data['num_jabr_cuts'] += counter
    all_data['jabr_cuts'].extend(new_cuts)
    all_data['t_jabr'].extend(new_t_jabr)
    all_data['c_jabr'].extend(new_c_jabr)
    all_data['s_jabr'].extend(new_s_jabr)
    all_data['vk_vm_jabr'].extend(new_vk_vm_jabr)

    #print("\n after")
    print("current number of jabr cuts",all_data['num_jabr_cuts'])
    print("current jabr cuts",all_data['jabr_cuts'])
    
    breakexit("cutjabr heuristic ok")
    

    
 
    

