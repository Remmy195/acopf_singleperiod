###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                        mjv2153@columbia.edu                               ##
##                                                                           ##
## April 2024                                                                ##
###############################################################################

import sys
import os
import numpy as np
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from ac_mtp2 import *
from ac_mtp import *
from socp_mtp import *

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        exit(1)

    casefilename       = '../data/none'
    modfile            = '../modfiles/none'
    solver             = 'knitroampl'
    fix_solution       = 0
    fix_tolerance      = 1e-5
    initial_solution   = 0
    initial_solution_0 = 0
    writesol           = 0

    # Solver parameters
    feastol            = "1e-6"
    opttol             = "1e-6"
    linsolver          = "5"
    max_time           = "1800"
    wstart             = '0'
    bar_initmu         = '1e-1'
    multistart         = 0
    
    T                  = 2
    expand             = 0

    AMPL_presolve      = 0

    nperturb           = 0
    uniform            = 1
    uniform_drift      = 0.02    
    uniform2           = 0
    ac                 = 0
    socp               = 0
    heuristic1         = 0
    heuristic2         = 0

    
    linenum = 0
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename  = thisline[1]

            elif thisline[0] == 'modfile':
                modfile       = thisline[1]
                
            elif thisline[0] == 'ac':
                ac   = 1
                socp = 0

            elif thisline[0] == 'socp':
                ac   = 0
                socp = 1

            elif thisline[0] == 'heuristic1':
                heuristic1   = 1
                heuristic2   = 0

            elif thisline[0] == 'heuristic2':
                heuristic1   = 0
                heuristic2   = 1                                
            
            elif thisline[0] == 'solver':
                solver        = thisline[1]

            elif thisline[0] == 'initial_solution':
                initial_solution = 1

            elif thisline[0] == 'initial_solution_0':
                initial_solution_0 = 1                

            elif thisline[0] == 'fix_solution':
                fix_solution     = 1

            elif thisline[0] == 'fix_tolerance':
                fix_tolerance    = 1                                

            elif thisline[0] == 'multistart':
                multistart    = 1

            elif thisline[0] == 'writesol':
                writesol    = 1

            elif thisline[0] == 'T':
                T           = int(thisline[1])

            elif thisline[0] == 'expand':
                expand           = 1

            elif thisline[0] == 'AMPL_presolve':
                AMPL_presolve = 1

            elif thisline[0] == 'feastol':
                feastol      = thisline[1]

            elif thisline[0] == 'opttol':
                opttol       = thisline[1]

            elif thisline[0] == 'linsolver':
                linsolver    = thisline[1]

            elif thisline[0] == 'max_time':
                max_time     = thisline[1]

            elif thisline[0] == 'wstart':
                wstart       = '1'

            elif thisline[0] == 'bar_initmu':
                bar_initmu   = thisline[1]                

            elif thisline[0] == 'nperturb':
                nperturb           = 1
                uniform            = 0
                uniform2           = 0

            elif thisline[0] == 'uniform':
                uniform           = 1
                uniform2          = 0
                nperturb          = 0

            elif thisline[0] == 'uniform_drift':
                uniform_drift     = float(thisline[1])
                
            elif thisline[0] == 'uniform2':
                uniform2           = 1
                uniform            = 0
                nperturb           = 0                                
                
            elif thisline[0] == 'END':
                break

            else:
                exit("illegal input " + thisline[0] + "bye")

        linenum += 1


    all_data                       = {}
    all_data['casefilename']       = casefilename
    all_data['casename']           = casefilename.split('data/')[1].split('.m')[0]
    all_data['modfile']            = modfile.split('../modfiles/')[1]
    all_data['solver']             = solver
    all_data['AMPL_presolve']      = AMPL_presolve
    all_data['initial_solution']   = initial_solution
    all_data['initial_solution_0'] = initial_solution_0    
    all_data['fix_solution']       = fix_solution
    all_data['fix_tolerance']      = fix_tolerance
    all_data['bar_initmu']         = bar_initmu
    all_data['multistart']         = multistart
    all_data['writesol']           = writesol
    all_data['T']                  = T
    all_data['expand']             = expand
    all_data['uniform2']           = uniform2
    all_data['uniform']            = uniform
    all_data['uniform_drift']      = uniform_drift    
    all_data['nperturb']           = nperturb
    all_data['feastol']            = feastol
    all_data['opttol']             = opttol
    all_data['linsolver']          = linsolver
    all_data['max_time']           = max_time
    all_data['wstart']             = wstart
    all_data['ac']                 = ac
    all_data['socp']               = socp
    all_data['heuristic1']         = heuristic1
    all_data['heuristic2']         = heuristic2    
    
    
    casetype = ''

    if all_data['nperturb']:
        casetype = 'gaussian'
        log.joint(' case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform']:
        casetype = 'uniform_'+str(uniform_drift)
        log.joint(' case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform2']:
        casetype = 'uniform2'
        log.joint(' case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')

    all_data['casetype'] = casetype
    
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

    if all_data['ac']:
        all_data['solver'] = 'knitroampl'
        log.joint('\n==================================================='
                  + '===========================\n')
        log.joint('===================================================='
                  + '==========================\n')        
        log.joint('\n                  Initializing Multi-Time Period ACOPF\n')
        if all_data['heuristic1']:
            all_data['modfile'] == 'ac_mtp.mod'
            log.joint('\n**Running Heuristic 1: We give to Knitro the '
                      + 'full Multi-Time Period Formulation\n')
            retcode = goac_mtp(log,all_data)
        elif all_data['heuristic2']:
            all_data['modfile'] == 'ac.mod'
            log.joint('\n\nRunning Heuristic 2: We give to Knitro '
                      + 'one period at a time,\n')
            log.joint('and impose the ramping constraints using '
                      + 'generation from t-1\n')            
            retcode = goac_mtp2(log,all_data)            
               
        # Set as default option when running MTP ACOPF
        else:
            log.joint('\n**Running Heuristic 3: First we run H1. If '
                      + 'it fails, relax tolerances and run H2.\n')
            log.joint(' feastol (abs) = ' + str(all_data['feastol'])
                      + ' opttol (abs) = ' + str(all_data['opttol'])
                      + ' max time (secs) ' + str(all_data['max_time'])
                      + '\n')
            all_data['modfile'] = 'ac_mtp.mod'
            retcode  = goac_mtp(log,all_data)
            if retcode != 0:
                log.joint('\n\nRunning Heuristic 2: We give to Knitro '
                          + 'one period at a time,\n')
                log.joint('and impose the ramping constraints using '
                          + 'generation from t-1\n')
                log.joint(' We relax feasibility and optimality tolerances\n')
                all_data['modfile']  = 'ac.mod'
                all_data['max_time'] = "1200"
                all_data['feastol']  = "1e-4"
                all_data['opttol']   = "1e-4"
                log.joint(' feastol (abs) = ' + str(all_data['feastol'])
                          + ' opttol (abs) = ' + str(all_data['opttol'])
                          + ' max time (secs) per period ' + str(all_data['max_time'])
                          + '\n')
                retcode = goac_mtp2(log,all_data)

    elif all_data['socp']:
        log.joint('\n==================================================='
                  + '===========================\n')
        log.joint('==================================================='
                  + '===========================\n')        
        log.joint('\n      Initializing the JABR SOCP relaxation for '
                  + 'Multi-Time Period ACOPF\n\n')
        if all_data['solver'] == 'mosek':
            all_data['modfile'] = 'jabr_mosek_mtp.mod'
            retcode = gosocp_mosek_mtp(log,all_data)
        else:
            all_data['modfile'] = 'jabr_mtp.mod'
            retcode = gosocp_mtp(log,all_data)
    else:
        log.joint(' modfile chosen ' + all_data['modfile'] + '\n')
        log.joint(' solver chosen ' + all_data['solver'] + '\n')
        log.joint('main_mtp: ac or socp? Please check config file.\n')
        exit(1)
        
    log.closelog()

