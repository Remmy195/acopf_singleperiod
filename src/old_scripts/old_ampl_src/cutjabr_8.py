import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def cutjabr_8(log,all_data):
    
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
    rnd = all_data['round']
    jabr_cuts = all_data['jabr_cuts']
    new_c_jabr = []
    new_s_jabr = []
    new_vf_jabr = []
    new_vt_jabr = []
    
    threshold = all_data['threshold']
    violated = {}
    violated_count = 0
    log.joint(" computing violations of jabr inequalities ... \n")
    for i in branches:
        violation = c[i]**2 + s[i]**2 - v[bus_f[i]] * v[bus_t[i]] 
        if violation > threshold:
            violated_count += 1
            violated[i] = violation
            
    log.joint(" computing jabr-cuts ...\n")

    
    num_selected =  math.ceil(violated_count * all_data['most_violated_fraction'] )
    most_violated = dict(sorted(violated.items(), key = lambda x: x[1], reverse = True)[:num_selected])
    most_violated_count = 0

    for i in most_violated.keys():
        if (most_violated_count == 0):
            most_violated_branch = i
            max_error = most_violated[i]
        most_violated_count += 1
        print("--------\n")
        print("branch = {0}, violation = {1}, cut id = {2}".format(i,most_violated[i],num_cuts + most_violated_count))
        print("variables current solution: cft = {0}, sft = {1}, vf = {2}, vt = {3}".format(c[i],s[i],v[bus_f[i]],v[bus_t[i]]))
        cutnorm = math.sqrt( (2 * c[i])**2 + (2 * s[i])**2 + (v[bus_f[i]] - v[bus_t[i]])**2 )
        coeff_cft = 4 * c[i]
        coeff_sft = 4 * s[i]
        coeff_vf = v[bus_f[i]] - v[bus_t[i]] - cutnorm
        coeff_vt = - (v[bus_f[i]] - v[bus_t[i]]) - cutnorm
        jabr_cuts[(num_cuts + most_violated_count,i)] = (rnd,most_violated[i],coeff_cft,coeff_sft,coeff_vf,coeff_vt)
        new_c_jabr.append(coeff_cft)
        new_s_jabr.append(coeff_sft)
        new_vf_jabr.append(coeff_vf)
        new_vt_jabr.append(coeff_vt)
        print("LHS coeff: cft = {0} , sft = {1}, vf = {2}, vt = {3}".format(coeff_cft,coeff_sft,coeff_vf,coeff_vt))
        print("cutnorm = {0}\n".format(cutnorm))

    #if all_data['drop_jabr_cuts'] == 0 and num_cuts:
    #    all_data['drop_jabr_cuts'] = 1
        
    all_data['most_violated'] = most_violated
    all_data['violated'] = most_violated_count
    all_data['num_jabr_cuts_rnd'][rnd] = most_violated_count
    
    if all_data['violated'] == 0:
        print(" largest violations below threshold {0}\n".format(all_data['threshold']))
        print(" no more cuts to add for current threshold\n")
        return None
    
    print("number of violated jabr's = {0} number of jabr-cuts added = {1}".format(violated_count,most_violated_count))
    print("max error = {0} at branch = {1}".format(max_error,most_violated_branch))
    log.joint(" saving jabr cuts to all_data ... \n")

    all_data['num_jabr_cuts'] += most_violated_count

    all_data['c_jabr'].extend(new_c_jabr)
    all_data['s_jabr'].extend(new_s_jabr)
    all_data['vf_jabr'].extend(new_vf_jabr)
    all_data['vt_jabr'].extend(new_vt_jabr)
    #breakexit("cutjabr heuristic ok")
    

    
 
    

