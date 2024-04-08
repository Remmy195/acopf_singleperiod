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
import time
from versioner import *
from log import danoLogger
from myutils import breakexit

if __name__ == '__main__':
    if len(sys.argv) != 1:
        print ('Usage: acopf_exp.py\n')
        exit(0)

    T0        = time.time()
    mylogfile = "exp_" + str(time.strftime('%m-%d-%Y_%H:%M:%S_%Z', time.localtime(T0))) + ".log"
    ac        = 1
    socp      = 1
    
    log = danoLogger(mylogfile)
    stateversion(log)

    if ac:

        log.joint("========================================\n")
        log.joint("========================================\n")
        log.joint("                ACOPF                   \n")
        log.joint("========================================\n")
        log.joint("========================================\n")

        t_ac = time.time()
        
        ACsols = 'ACsols_' + str(time.strftime('%m-%d-%Y_%H:%M:%S_%Z', time.localtime(t_ac)))
        os.system("mkdir " + ACsols) 

        directory = '../data/'
        files = os.listdir(directory)
        index = 0
        oldfilename = "case0.m"
        os.system("sed -i s/solver0/knitroampl/ ac_mtp.conf")
        os.system("sed -i s/nlp/ac/ ac_mtp.conf")        
        while index < len(files):
            filename = files[index]
            if filename.endswith('.m'):
                log.joint('filename ' + str(filename) + '\n')
                os.system("sed -i s/" + oldfilename + "/" + filename + "/ ac_mtp.conf")
                log.joint('running program with ' + filename + '\n')
                os.system("echo " + filename) 
                os.system("python3 ../src/main_mtp.py ac_mtp.conf " + ACsols)
                oldfilename = filename
            index += 1

        os.system("sed -i s/" + oldfilename + "/case0.m/ ac_mtp.conf")
        os.system("sed -i s/knitroampl/solver0/ ac_mtp.conf")
        os.system("sed -i s/ac/nlp/ ac_mtp.conf")                
                  
        log.joint('done with ACOPF experiments\n')
        log.joint('total time = ' + str(time.time() - t_ac) + '\n')

    if socp:
        log.joint("========================================\n")
        log.joint("========================================\n")
        log.joint("                SOCPS                   \n")
        log.joint("========================================\n")
        log.joint("========================================\n")

        t_socp = time.time()
        
        os.system("sed -i s/nlp/socp/ ac_mtp.conf")
        solvers = ["knitroampl","gurobi_ampl","mosek"]
        directory = '../data/'
        files = os.listdir(directory)
        oldsolver   = "solver0"
        for solver in solvers:
            log.joint("========================================\n")            
            log.joint('\n\nRunning solver ' + solver + '...\n')
            os.system("sed -i s/" + oldsolver + "/" + solver + "/ ac_mtp.conf")
            SOCPsols = 'SOCPsols_' + solver + "_" + str(time.strftime('%m-%d-%Y_%H:%M:%S_%Z', time.localtime(t_ac)))
            os.system("mkdir " + SOCPsols)
            index = 0
            oldfilename = "case0.m"            
            while index < len(files):
                filename = files[index]
                if filename.endswith('.m'):
                    log.joint('filename ' + str(filename) + '\n')
                    os.system("sed -i s/" + oldfilename + "/" + filename + "/ ac_mtp.conf")
                    log.joint('running program with ' + filename + '\n')
                    if index == 10:
                        breakexit('c')
                    os.system("echo " + filename) 
                    os.system("python3 ../src/main_mtp.py ac_mtp.conf " + SOCPsols) # I could pass the exp.log file here as well
                    oldfilename = filename
                index += 1

            t1    = time.time()
            oldsolver  = solver
            os.system("sed -i s/" + oldfilename + "/case0.m/ ac_mtp.conf")

            log.joint('done with SOCP experiments with solver ' + solver + '\n')
            log.joint('total time = ' + str(time.time() - t_socp) + '\n')


        os.system("sed -i s/" + oldsolver + "/solver0/ ac_mtp.conf")
        os.system("sed -i s/socp/nlp/ ac_mtp.conf")                
        log.joint('done with SOCP experiments\n')
        log.joint('total time = ' + str(time.time() - t_socp) + '\n')

    T1 = time.time()
    log.joint('done with experiments\n')
    log.joint('total time = ' + str(T1 - T0) + '\n')
        
    log.closelog()

