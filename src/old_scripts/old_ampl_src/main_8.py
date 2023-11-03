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
from goampl_8 import goampl_8


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
    modfile = 'NONE'
    modfile_cut = 'NONE'
    lpfilename = 'lp.lp'
    solver = 'ipopt'
    solver_cut = 'gurobi_ampl'
    linenum = 0
    threshold = 1e-05
    most_violated_fraction = 1
    obbt = 0
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]

            elif thisline[0] == 'modfile':
                modfile = thisline[1]

            elif thisline[0] == 'modfile_cut':
                modfile_cut = thisline[1]
                
            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]

            elif thisline[0] == 'solver':
                solver = thisline[1]

            elif thisline[0] == 'solver_cut':
                solver_cut = thisline[1]

            elif thisline[0] == 'threshold':
                threshold = float(thisline[1])

            elif thisline[0] == 'most_violated_fraction':
                most_violated_fraction = float(thisline[1])

            elif thisline[0] == 'obbt':
                obbt = 1
                
            elif thisline[0] == 'END':
                break

            else:
                sys.exit("illegal input " + thisline[0] + "bye")

        linenum += 1

    all_data = {}
    all_data['timestart'] = time.time()
    all_data['casefilename'] = casefilename
    all_data['modfile'] = modfile
    all_data['modfile_cut'] = modfile_cut
    all_data['lpfilename'] = lpfilename
    all_data['solver'] = solver
    all_data['solver_cut'] = solver_cut
    all_data['obbt'] = obbt
    
    all_data['mincut_cuts'] = [] ####
    all_data['num_mincut_cuts'] = 0
    all_data['new_mincut'] = 0 ####
    all_data['size_supernodes'] = 0 ##
    
    all_data['most_violated_fraction'] = most_violated_fraction
    all_data['threshold'] = threshold
    
    all_data['num_jabr_cuts'] = 0
    #all_data['drop_jabr_cuts'] = 1
    all_data['jabr_cuts'] = {}
    all_data['c_jabr'] = []
    all_data['s_jabr'] = []
    all_data['vf_jabr'] = []
    all_data['vt_jabr'] = []
    all_data['num_jabr_cuts_rnd'] = {}

    
    all_data['obbt_min'] = 0
    
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

    goampl_8(log,all_data)    

    log.closelog()
