import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def cutjabr_3(log,all_data):


    log.joint(" loading data from all_data\n")

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
    print("dictionary of violations =",violation)

    log.joint(" checking for most violated jabrs ... \n")

    N =  math.ceil(violated/3)
    
    most_violated = dict(sorted(violation.items(), key = lambda x: x[1], reverse = True)[:N])

    print("the N = {0} most violated cuts = {1}".format(N,most_violated))
    
    breakexit("check")

    log.joint(" saving data to all_data ... \n")

    print("most violated jabr constraints =",list(most_violated.keys()))

    all_data['jabr_cuts'] = list(most_violated.keys())
    all_data['t_jabr'] = [math.sqrt( (2 * c[i])**2 + (2 * s[i])**2 + (v[bus_f[i]] - v[bus_t[i]])**2 ) for i in most_violated.keys()]
    all_data['c_jabr'] = [2 * c[i] for i in most_violated.keys()]
    all_data['s_jabr'] = [2 * s[i] for i in most_violated.keys()]
    all_data['vk_vm_jabr'] =  [v[bus_f[i]] - v[bus_t[i]] for i in most_violated.keys()]

    print("t_jabr",all_data['t_jabr'])

    
 
    

