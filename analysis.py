from pytriqs.operators.operators import Operator, c, c_dag, n, dagger
from pytriqs.archive import *
from pytriqs.applications.impurity_solvers.cthyb import SortedSpaces

f = HDFArchive('dca_beta20_np975_histogram.out.h5','r')
h_diag =  f['h_diag-0']

fc = HDFArchive('check.h5','w')
fc['h_diag-0'] = h_diag 

A = PostProcess(f['density_matrix-0'],f['h_loc-0'])

CHEC = HDFArchive('densmat.h5','w')
CHEC['density_matrix'] = f['density_matrix-0']


c_dag_0_up = ( c_dag("00-up",0) + c_dag("10-up",0) + c_dag("01-up",0) + c_dag("11-up",0))/2
c_dag_1_up = ( c_dag("00-up",0) - c_dag("10-up",0) + c_dag("01-up",0) - c_dag("11-up",0))/2
c_dag_2_up = ( c_dag("00-up",0) - c_dag("10-up",0) - c_dag("01-up",0) + c_dag("11-up",0))/2
c_dag_3_up = ( c_dag("00-up",0) + c_dag("10-up",0) - c_dag("01-up",0) - c_dag("11-up",0))/2

c_dag_0_down = ( c_dag("00-down",0) + c_dag("10-down",0) + c_dag("01-down",0) + c_dag("11-down",0))/2
c_dag_1_down = ( c_dag("00-down",0) - c_dag("10-down",0) + c_dag("01-down",0) - c_dag("11-down",0))/2
c_dag_2_down = ( c_dag("00-down",0) - c_dag("10-down",0) - c_dag("01-down",0) + c_dag("11-down",0))/2
c_dag_3_down = ( c_dag("00-down",0) + c_dag("10-down",0) - c_dag("01-down",0) - c_dag("11-down",0))/2

ops_0hole = [
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)/2,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)/2
]

ops_1hole = [
c_dag_0_up*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up),
c_dag_0_down*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up),
(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)*c_dag_3_up,
(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)*c_dag_3_down,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_2_up,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_2_down,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_1_up,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_1_down,
c_dag_0_up*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up),
c_dag_0_down*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up),
(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)*c_dag_1_up,
(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)*c_dag_1_down,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_2_up,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_2_down,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_3_up,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_3_down
]

# n in space
n0up = c_dag_0_up * dagger(c_dag_0_up)
n1up = c_dag_1_up * dagger(c_dag_1_up)
n2up = c_dag_2_up * dagger(c_dag_2_up)
n3up = c_dag_3_up * dagger(c_dag_3_up)
print  "n0up = ",  A.average(n0up) 
print  "n1up = ",  A.average(n1up) 
print  "n2up = ",  A.average(n2up) 
print  "n3up = ",  A.average(n3up) 

n0down = c_dag_0_down * dagger(c_dag_0_down)
n1down = c_dag_1_down * dagger(c_dag_1_down)
n2down = c_dag_2_down * dagger(c_dag_2_down)
n3down = c_dag_3_down * dagger(c_dag_3_down)
print  "n0down = ",  A.average(n0down) 
print  "n1down = ",  A.average(n1down) 
print  "n2down = ",  A.average(n2down) 
print  "n3down = ",  A.average(n3down) 

# n in k space (to compare with gf)
n_00 =  n("00-up", 0)
n_10 =  n("10-up", 0)
n_01 =  n("01-up", 0)
n_11 =  n("11-up", 0)
print  "<n00> = ",  A.average(n_00) 
print  "<n10> = ",  A.average(n_10) 
print  "<n01> = ",  A.average(n_01) 
print  "<n11> = ",  A.average(n_11) 

# dimer states are NOT orthogonals.
print  "dimer prod" ,  A.dot_product_from_creation_op(ops_0hole[0], ops_0hole[0]) 
print  "dimer prod" ,  A.dot_product_from_creation_op(ops_0hole[0], ops_0hole[1]) 
print  "dimer prod" ,  A.dot_product_from_creation_op(ops_0hole[1], ops_0hole[1]) 
print  "dimer prod" ,  A.dot_product_from_creation_op(ops_1hole[1], ops_0hole[1]) 

for op in ops_1hole:
   print  "<op0|op>" ,  A.dot_product_from_creation_op(op, ops_1hole[2]) 

for n, op in enumerate(ops_0hole): 
    print n, A.average_projector_from_creation_op(op)

print "\n"

for n, op in enumerate(ops_1hole): 
    print n, A.average_projector_from_creation_op(op)/2
