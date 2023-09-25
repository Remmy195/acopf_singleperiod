#!/usr/bin/python                                                                                                                         
import sys
import os
import numpy as np
import subprocess
import time
import reader
import pglibreader
from myutils import breakexit
from versioner import *
from log import danoLogger
from gosocp import *
from gosocp2 import gosocp2

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
    lpfilename        = 'lp.lp'
    solver            = 'ipopt'
    fix_point         = 0
    multistart        = 0
    conic             = 0
    knitropresolveoff = 0
    mytol             = 0
    
    linenum = 0

    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]

            elif thisline[0] == 'modfile':
                modfile = thisline[1]
                
            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]

            elif thisline[0] == 'solver':
                solver = thisline[1]

            elif thisline[0] == 'fix_point':
                fix_point  = 1

            elif thisline[0] == 'multistart':
                multistart = 1

            elif thisline[0] == 'conic':
                conic = 1

            elif thisline[0] == 'knitropresolveoff':
                knitropresolveoff = 1

            elif thisline[0] == 'mytol':
                mytol = 1
                
            elif thisline[0] == 'END':
                break

            else:
                sys.exit("illegal input " + thisline[0] + "bye")

        linenum += 1

    all_data                      = {}
    all_data['casefilename']      = casefilename
    all_data['casename']          = casefilename.split('data/')[1].split('.m')[0]
    all_data['modfile']           = modfile
    all_data['lpfilename']        = lpfilename
    all_data['solver']            = solver
    all_data['fix_point']         = fix_point
    all_data['multistart']        = multistart
    all_data['conic']             = conic
    all_data['knitropresolveoff'] = knitropresolveoff 
    all_data['mytol']             = mytol
    
    return all_data

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: main.py file.config [logfile]\n')
        exit(0)
        
    T0 = time.time()
    mylogfile = "main.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = danoLogger(mylogfile)
    stateversion(log)

    all_data       = read_config(log,sys.argv[1])
    all_data['T0'] = T0
    
    #readcode       = reader.readcase(log,all_data,all_data['casefilename'])
    readcode       = pglibreader.readcase(log,all_data,all_data['casefilename'])

    if all_data['modfile'] == 'jabr.mod':
        gosocp(log,all_data)
    elif all_data['modfile'] == 'jabr_i2.mod':
        gosocp2(log,all_data)    
    else:
        log.joint(' check config file\n')
        exit(0)
        
    log.closelog()
