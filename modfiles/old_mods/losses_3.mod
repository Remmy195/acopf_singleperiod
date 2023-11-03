#################

#A conic relaxation (Jabr) of ACOPF 
  
#################

param nI >= 0, integer; #number of buses
param nJ >= 0, integer; #number of branches
param nG >= 0, integer; #number of generators
                                                              
#SETS                                
set buses;
set gens;
set branches;
set branches_f {i in buses};
set branches_t {i in buses};
set bus_gens {i in buses};
set cuts;
set jabr_cuts;
#set jabr_index
     
#PARAMETERS
param fixedcost {i in gens} >= 0;
param lincost {i in gens}; #linear cost each generator (for the moment, then add option for piecewise linear)    
param quadcost {i in gens} >= 0; #quadratic cost coefficient each generator

param Gtt {i in branches}; 
param Btt {i in branches};
param Gff {i in branches}; 
param Bff {i in branches}; 
param Gtf {i in branches}; 
param Btf {i in branches}; 
param Gft {i in branches}; 
param Bft {i in branches}; 
param Rft {i in branches};

param bus_f {i in branches};
param bus_t {i in branches};                                                                                  
param Vmax {i in buses} >= 0; #max voltage (non squared)
param Vmin {i in buses} >= 0; #min voltage (non squared)
param CSmax {i in branches} >= 0;
param Pmax {i in gens}; #max active power gen
param Pmin {i in gens}; #min active power gen
param Qmax {i in gens}; #max reactive power gen
param Qmin {i in gens}; #min reactive power gen                                                    
param Pd {i in buses}; #active power demand
param Qd {i in buses}; #reactive power demand
param U {i in branches} >= 0; #branch power flow limit

param t_jabr {i in jabr_cuts};
param c_jabr {i in jabr_cuts};
param s_jabr {i in jabr_cuts};
param vk_vm_jabr {i in jabr_cuts};

    
#VARIABLES                                
var c {i in branches} >= - CSmax[i], <= CSmax[i];
var s {i in branches} >= - CSmax[i], <= CSmax[i]; 
var v {i in buses} >= Vmin[i] ^2, <= Vmax[i] ^2; #squared-voltage magnitudes (1:= 1?)
var Pf {i in branches} >= - U[i], <= U[i]; #limit branches
var Pt {i in branches} >= - U[i], <= U[i]; #limit branches
var Qf {i in branches} >= - U[i], <= U[i]; #limit branches
var Qt {i in branches} >= - U[i], <= U[i]; #limit branches                     
var PLoss {i in branches} = Pf[i] + Pt[i];
var I {i in branches} >= 0;
                                  
var Pg {i in gens} >= Pmin[i], <= Pmax[i]; #active power generator
var Qg {i in gens} >= Qmin[i], <= Qmax[i]; #reactive power generator
                                                                                                  
#OBJECTIVE                                                  
minimize total_cost: sum {i in gens} (fixedcost[i] + lincost[i] * Pg[i] + quadcost[i] * Pg[i] ^2);
           
#CONSTRAINTS

#def
subject to Pf_def {i in branches}:
                   Pf[i] = Gff[i] * v[bus_f[i]] + Gft[i] * c[i] + Bft[i] * s[i];

subject to Pt_def {i in branches}: 
                   Pt[i] = Gtt[i] * v[bus_t[i]] + Gtf[i] * c[i] - Btf[i] * s[i];

subject to Qf_def {i in branches}:
                   Qf[i] = - Bff[i] * v[bus_f[i]] - Bft[i] * c[i] + Gft[i] * s[i];

subject to Qt_def {i in branches}:
                   Qt[i] = - Btt[i] * v[bus_t[i]] - Btf[i] * c[i] - Gtf[i] * s[i];

#def_I

subject to I_def {i in branches: Rft[i] > 0}:
                   I[i] = (1/Rft[i]) * ( Gff[i] * v[bus_f[i]] + Gtt[i] * v[bus_t[i]] + c[i] * ( Gft[i] + Gtf[i] ) + s[i] * ( Bft[i] - Btf[i] ) ); 

#balance

subject to Pbalance {i in buses}:
         (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k] ) else 0) - Pd[i];

subject to Qbalance {i in buses}:
         (sum {j in branches_f[i]} Qf[j] ) + (sum {j in branches_t[i]} Qt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k] ) else 0) - Qd[i];

#Loss inequalities                                                                                                       
subject to Loss_ineq {i in branches}: Pf[i] + Pt[i] >= 0;

#subject to cs_linear {i in branches}: 2 * c[i] <= v[bus_f[i]] + v[bus_t[i]];

subject to linear_cuts_f {i in cuts}: 2 * Pf[i] <= v[bus_f[i]] + I[i];

subject to linear_cuts_t {i in cuts}: 2 * Pt[i] <= v[bus_t[i]] + I[i];

#subject to conic_cuts {i in cuts}: Pf[i] ^2 + Qf[i] ^2 <= v[bus_f[i]] * I[i];

#subject to jabr {i in branches}: c[i] ^2 + s[i] ^2 <= v[bus_f[i]] * v[bus_t[i]];

subject to j_cuts {i in jabr_cuts}: 
    t_jabr[i] * ( v[bus_f[i]] + v[bus_t[i]] ) >= c_jabr[i] * 2 * c[i] + s_jabr[i] * 2 * s[i] + vk_vm_jabr[i] * (v[bus_f[i]] - v[bus_t[i]]);

#subject to NewLoss1_2: PLoss[1] + v[bus_f[1]] >= 2 * Rft[1] ^0.5 * Pf[1];
#subject to NewLossP {i in branches}: PLoss[i] + v[bus_f[i]] >= 2 * Rft[i] ^0.5 * Pf[i];
#subject to NewLossQ {i in branches}: PLoss[i] + v[bus_f[i]] >= 2 * Rft[i] ^0.5 * Qf[i];


#limits
#subject to limits_f {i in branches}: Pf[i] ^2 + Qf[i] ^2 <= U[i] ^2;

#subject to limits_t {i in branches}: Pt[i] ^2 + Qt[i] ^2 <= U[i] ^2;
             

#reference bus
