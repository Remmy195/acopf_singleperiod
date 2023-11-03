#!/usr/bin/python                                                                                                                         
import sys
import os
import numpy as np
import subprocess
import time
import reader
from myconstants import setconstants
from myutils import breakexit
from versioner import *
from log import danoLogger
from golp import golp

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    casefilename = 'NONE'
    lpfilename = 'NOE'
    linenum = 0
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]
                
            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]
                
            elif thisline[0] == 'END':
                break

            else:
                sys.exit("illegal input " + thisline[0] + "bye")

        linenum += 1

    all_data = {}
    all_data['casefilename'] = casefilename
    all_data['lpfilename'] = lpfilename
    all_data['cuts'] = []
    
    return all_data

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: main.py file.config [logfile]\n')
        exit(0)
    t0 = time.time()

    mylogfile = "main.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = danoLogger(mylogfile)
    stateversion(log)

    all_data = read_config(log,sys.argv[1])
    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    setconstants(log,all_data)

    golp(log,all_data)
    
    log.closelog()
