###############################################################################
##                                                                           ##      
## This code was written and is being maintained by Matias Villagra,         ##     
## PhD Student in Operations Research @ Columbia, supervised by              ## 
## Daniel Bienstock.                                                         ##      
##                                                                           ##    
## Please report any bugs or issues (for sure there will be) to              ##     
##                        mjv2153@columbia.edu                               ##      
##                                                                           ##  
## Oct 2023                                                                  ## 
###############################################################################

import sys
import os
import numpy as np
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from ac import *
from socp import *
from graphics import *

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    casefilename      = 'NONE'
    modfile           = 'NONE'
    solver            = 'knitro'
    fix               = 0
    fix_point         = 0
    initial_point     = 0
    multistart        = 0
    mytol             = 0
    writesol          = 0
    knitro_algorithm = 1         # default to Interior-Direct
    knitro_bar_murule = None
    knitropresolveoff = 0

    linenum       = 0
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename  = thisline[1]

            elif thisline[0] == 'modfile':
                modfile       = thisline[1]
                
            elif thisline[0] == 'lpfilename':
                lpfilename    = thisline[1]

            elif thisline[0] == 'solver':
                solver        = thisline[1]

            elif thisline[0] == 'initial_point':
                initial_point = 1

            elif thisline[0] == 'fix':
                fix           = 1

            elif thisline[0] == 'fix_point':
                fix_point     = 1                

            elif thisline[0] == 'multistart':
                multistart    = 1

            elif thisline[0] == 'mytol':
                mytol    = 1

            elif thisline[0] == 'writesol':
                writesol    = 1                                

            elif thisline[0] == 'knitropresolveoff':
                knitropresolveoff = 1

            elif thisline[0] == 'knitro_algorithm':
                knitro_algorithm = int(thisline[1])

            elif thisline[0] == 'knitro_bar_murule':
                knitro_bar_murule = int(thisline[1])

            elif thisline[0] == 'END':
                break

            else:
                sys.exit("illegal input " + thisline[0] + "bye")

        linenum += 1

    all_data                      = {}
    all_data['casefilename']      = casefilename
    all_data['casename']          = casefilename.split('data/')[1].split('.m')[0]
    all_data['modfile']           = modfile.split('../modfiles/')[1]
    all_data['solver']            = solver
    all_data['initial_point']     = initial_point
    all_data['fix_point']         = fix_point    
    all_data['fix']               = fix
    all_data['multistart']        = multistart
    all_data['mytol']             = mytol
    all_data['writesol']          = writesol
    all_data['knitropresolveoff'] = knitropresolveoff    
    all_data['knitro_algorithm']  = knitro_algorithm
    all_data['knitro_bar_murule'] = knitro_bar_murule
    
    return all_data

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: main.py file.config [logfile]\n')
        exit(0)

    T0        = time.time()
    mylogfile = "main.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = danoLogger(mylogfile)
    stateversion(log)

    all_data       = read_config(log,sys.argv[1])
    all_data['T0'] = T0
    
    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    if all_data['modfile'] == 'ac.mod' or all_data['modfile'] == 'ac2.mod':
        s = all_data['solver'].lower()
        if not (s.startswith('knitro') or s.startswith('baron') or s.startswith('ipopt') or s.startswith('gurobi')):
            log.joint('Wrong solver/modfile, please check config file\n')
            exit(1)
        print('Objective value: ')
        goac(log,all_data)
        

    elif all_data['modfile'] == 'acopfqcqp.mod' or all_data['modfile'] == 'jabr.mod':
        if all_data['solver'] != 'mosek':
            gosocp(log, all_data)
        else:
            log.joint('Wrong solver/modfile, please check config file\n')
            exit(1)

    elif all_data['modfile'] == 'jabr_mosek.mod':
        if all_data['solver'] == 'mosek':
            gosocp_mosek(log,all_data)
        else:
            log.joint('Wrong solver/modfile, please check config file\n')
            exit(1)            
            
    elif all_data['modfile'] == 'i2.mod':
        if all_data['solver'] != 'mosek':
            gosocp2(log,all_data)
        else:
            log.joint('Wrong solver/modfile, please check config file\n')
            exit(1)

    elif all_data['modfile'] == 'i2_mosek.mod':
        if all_data['solver'] == 'mosek':
            gosocp2_mosek(log,all_data)
        else:
            log.joint('Wrong solver/modfile, please check config file\n')
            exit(1)
            
    else:
        log.joint('Wrong solver/modfile, please check config file\n')
        exit(1)

    solution = build_solution_dict_from_ampl(log,all_data)
    solution_plot(all_data, solution, None)



    log.closelog()

