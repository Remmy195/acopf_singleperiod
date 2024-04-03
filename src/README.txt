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

Quick guide on how to run the scripts:

- To run Multi-Time Period ACOPF:
1) Edit config file 'ac_mtp.conf' in directory 'runs':
1.1) Provide the casefilename of the instance;
1.2) If the string "ac" is given, the script runs heuristic3 by default (see
   below);
1.3) If "heuristic1" is given, the script runs heuristic1,
   if "heuristic2" is given, then it runs heuristic2.
2) In 'runs', type 'python ../src/main_mtp.py ac_mpt.conf'
Some heuristics for finding AC primal bounds using Knitro:

heuristic1: We give to Knitro the full Multi-Time Period Formulation
heuristic2: We give to Knitro one period at a time and impose
	    the ramping constraints using generation from the previous
	    time period
heuristic3: First we run heuristic1. If it fails, relax tolerances and
	    run heuristic2	    

- To run the Jabr relaxation of Multi-Time Period ACOPF:
1) Edit config file 'ac_mtp.conf' in directory 'runs':
1.1) Provide the casefilename of the instance;
1.2) Provide the string 'socp' to run the Jabr SOCP
1.3) Provide the string "solver " + "<solvername>"
2) In 'runs', type 'python ../src/main_mtp.py ac_mpt.conf'

Currently, the data provided corresponds to 4 time periods (c.f. directories
/data/ramprates and /data/mtploads) and a small number of instances from
PGLIB, the Pegase project, and some ACTIVSg big cases.
