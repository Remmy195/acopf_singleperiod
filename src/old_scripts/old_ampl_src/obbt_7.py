import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math

def obbt_7(log,all_data):
    
    log.joint(" ***********************\n")

    log.joint(" bound-tightening  heuristic\n")
    
    log.joint(" ***********************\n")
    
    log.joint("\n loading data from all_data\n")

    buses = all_data['goampl_buses'].copy()
    branches = all_data['goampl_branches'].copy()
    bus_f = all_data['goampl_bus_f'].copy() 
    bus_t = all_data['goampl_bus_t'].copy()
    most_violated = all_data['most_violated']

    #only for %of most_violated
    
    Vmax = all_data['goampl_Vmax']
    Vmin = all_data['goampl_Vmin']
    
    log.joint(" calling ampl object ... ")
    
    ampl = all_data['ampl_object']

    #solver_obbt = all_data['solver_obbt']
    print("most_violated",most_violated)
    for branch in most_violated.keys():
        bus = bus_f[branch]
        print("from bus of most violated branch =",bus)
        #first we minimize
        ampl.get_parameter('obbt_bus').set(bus)
        ampl.get_parameter('obbt_min').set(1)
        #ampl.eval("expand;")
        #breakexit("check")
        t0 = time.time()
        ampl.solve()
        t1 = time.time()
        print(" time to solve = {0}".format(t1-t0))
        v = ampl.get_variable("v")
        vdic = v.get_values().to_dict()
        print("new Vmin[{0}] = {1}".format(bus,vdic[bus]))
        Vmin[bus] = math.sqrt( vdic[bus] )

        ampl.get_parameter('obbt_min').set(-1)
        #ampl.eval("expand;")
        #breakexit("check")
        t0 = time.time()
        ampl.solve()
        t1 = time.time()
        print(" time to solve = {0}".format(t1-t0))
        v = ampl.get_variable("v")
        vdic = v.get_values().to_dict()
        print("new Vmax[{0}] = {1}".format(bus,vdic[bus]))
        Vmax[bus] = math.sqrt( vdic[bus] )
        

    all_data['goampl_Vmax'] = Vmax
    all_data['goampl_Vmin'] = Vmin
        
    ampl.get_parameter('Vmax').setValues(Vmax)
    ampl.get_parameter('Vmin').setValues(Vmin)

    print(" setting obbt_min = 0\n")

    all_data['obbt_min'] = 0
    ampl.get_parameter('obbt_min').set(all_data['obbt_min'])
    print(" done with obbt round\n")
    #breakexit("check")

    
    

    
 
    

