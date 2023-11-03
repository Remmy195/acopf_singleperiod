import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def cutjabr_6(log,all_data):
    
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

    num_cuts = all_data['num_jabr_cuts']
    new_cuts = []
    new_c_jabr = []
    new_s_jabr = []
    new_vf_jabr = []
    new_vt_jabr = []
    
    violation = {}
    violated = 0
    
    log.joint(" computing new jabr cuts... \n")
    for i in branches:
        violation[i] = c[i]**2 + s[i]**2 - v[bus_f[i]] * v[bus_t[i]] 
        if violation[i] > 1e-05:
            violated += 1
            print("branch = {0}, violation = {1}, cut id = {2}".format(i,violation[i],num_cuts + violated))
            print("variables: cft = {0}, sft = {1}, vf = {2}, vt = {3}".format(c[i],s[i],v[bus_f[i]],v[bus_t[i]]))
            cutnorm = math.sqrt( (2 * c[i])**2 + (2 * s[i])**2 + (v[bus_f[i]] - v[bus_t[i]])**2 )
            new_cuts.append((num_cuts + violated,i))
            coeff_cft = 4 * c[i]
            coeff_sft = 4 * s[i]
            coeff_vf = v[bus_f[i]] - v[bus_t[i]] - cutnorm
            coeff_vt = - (v[bus_f[i]] - v[bus_t[i]]) - cutnorm
            new_c_jabr.append(coeff_cft)
            new_s_jabr.append(coeff_sft)
            new_vf_jabr.append(coeff_vf)
            new_vt_jabr.append(coeff_vt)
            print("LHS coeff: cft = {0} , sft = {1}, vf = {2}, vt = {3}".format(coeff_cft,coeff_sft,coeff_vf,coeff_vt))
            print("cutnorm = {0}\n".format(cutnorm))
            print("-----\n")

    all_data['violated'] = violated
    
    if all_data['violated'] == 0:
        print(" no more cuts to add")
        return None

    most_violated_branch = max(violation, key=violation.get) 
        
    print("number of violated jabr's =",violated)
    print("max error = {0} at branch = {1}".format(violation[most_violated_branch],most_violated_branch))
    log.joint(" saving jabr cuts to all_data ... \n")

    all_data['num_jabr_cuts'] += violated
    all_data['jabr_cuts'].extend(new_cuts)
    all_data['c_jabr'].extend(new_c_jabr)
    all_data['s_jabr'].extend(new_s_jabr)
    all_data['vf_jabr'].extend(new_vf_jabr)
    all_data['vt_jabr'].extend(new_vt_jabr)
    #breakexit("cutjabr heuristic ok")
    

    
 
    

