

#################

#Cut heuristic
  
#################
                                                              
#SETS                                
set buses;
set branches;
set branches_f {i in buses};
set branches_t {i in buses};
set delta_s;
set delta_t;
set st;

#set gens;
#set bus_gens {i in buses};


#PARAMETERS

param bus_f {i in branches};
param bus_t {i in branches};                                                                                  
param PLoss {i in branches};
#param Pd {i in buses}; #active power demand
#param Pg {i in gens}; #active power generation
             
#VARIABLES                                
var Pf {i in branches}; 
var Pt {i in branches}; 

#OBJECTIVE                                                  
maximize flow: sum {i in delta_s} Pf[i];
           
#CONSTRAINTS

#balance

subject to Pbalance {i in buses: i not in st}:
         (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = 0;

#subject to Pbalance {i in buses: i not in st}:
#         (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k] ) else 0) - Pd[i];

#netflow

subject to Net {i in branches}: Pf[i] + Pt[i] = 0;

#capacities

subject to Cap_f {i in branches}: Pf[i] <= PLoss[i];

subject to Cap_t {i in branches}: Pt[i] <= PLoss[i];
            
#delta_t

subject to non_neg_delta_t {i in delta_t}: Pf[i] >= 0;
 
#losses for super_source

#subject to loss_super_source {j in delta_s}: Pf[j] + Pt[j] = 0;

