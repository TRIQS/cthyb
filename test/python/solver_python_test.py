import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators.operators2 import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
import itertools
from collections import OrderedDict

# H_loc parameters
beta = 10.0
num_orbitals = 2
mu = 1.0
U = 2.0
J = 0.2
n_n_only = False

# Poles of delta
epsilon = 2.3

# Hybridization matrices
V = 1.0 * np.eye(num_orbitals) + 0.1 * (np.ones(num_orbitals) - np.eye(num_orbitals))

# Block structure of GF
gf_struct = OrderedDict()
gf_struct['up'] = range(0,num_orbitals)
gf_struct['down'] = range(0,num_orbitals)

# Construct solver    
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=10001)

# Hamiltonian
H = Operator() 

for o in range(0,num_orbitals):
    H += U*n("up",o)*n("down",o)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*n("up",o1)*n("down",o2)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o2>=o1: continue;
    H += (U-3*J)*n("up",o1)*n("up",o2)
    H += (U-3*J)*n("down",o1)*n("down",o2)
      
if not n_n_only: # spin flips and pair hopping
    for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
        if o1==o2: continue
        H += -J*c_dag("up",o1)*c_dag("down",o1)*c("up",o2)*c("down",o2)
        H += -J*c_dag("up",o1)*c_dag("down",o2)*c("up",o2)*c("down",o1)

# Quantum numbers
qn = [Operator(),Operator()]
for o in range(0,num_orbitals):
    qn[0] += n("up",o)
    qn[1] += n("down",o)

# Set hybridization function
delta_w = GfImFreq(indices = range(0,num_orbitals), beta=beta)
delta_w <<= inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon)
delta_w.from_L_G_R(V, delta_w, V)

S.G0_iw["up"] <<= inverse(iOmega_n + mu - delta_w)
S.G0_iw["down"] <<= inverse(iOmega_n + mu - delta_w)

# Parameters
p = SolverCore.solve_parameters()
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50
p["n_cycles"] = 500

S.solve(h_loc=H, params=p, quantum_numbers=qn, use_quantum_numbers=True)
  
if mpi.rank==0:
    Results = HDFArchive("solver_python_test.output.h5",'w')
    Results["G_up"] = S.G_tau["up"]
    Results["G_down"] = S.G_tau["down"]
