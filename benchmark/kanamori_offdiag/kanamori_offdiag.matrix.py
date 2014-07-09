#!/bin/env pytriqs

from datetime import datetime
import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb_matrix import *
from pytriqs.gf.local import *
from itertools import product

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

# Remove cthyb-specific parameters 
cthyb_only_params = [
    'beta',
    'verbosity']
p = {k:v for k, v in p.items() if not k in cthyb_only_params}

print_master("Welcome to Kanamori benchmark, matrix version.")

gf_struct = {}
for spin in spin_names:
    for o in range(num_orbitals):
        bn, i = mkind(spin,o)
        gf_struct.setdefault(bn,[]).append(i)

# Hamiltonian
H = Operator()

for o in range(num_orbitals):
    H += U *N(*mkind("up",o))*N(*mkind("dn",o))
    
for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*N(*mkind("up",o1))*N(*mkind("dn",o2))

for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if o2>=o1: continue
    H += (U-3*J)*N(*mkind("up",o1))*N(*mkind("up",o2))
    H += (U-3*J)*N(*mkind("dn",o1))*N(*mkind("dn",o2))

for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if(o1==o2): continue
    H += -J*Cdag(*mkind("up",o1))*Cdag(*mkind("dn",o1))*C(*mkind("up",o2))*C(*mkind("dn",o2))
    H += -J*Cdag(*mkind("up",o1))*Cdag(*mkind("dn",o2))*C(*mkind("up",o2))*C(*mkind("dn",o1))

QN={}
if use_qn:
    QN["N_up"] = sum([N(*mkind("up",o)) for o in range(num_orbitals)],Operator())
    QN["N_dn"] = sum([N(*mkind("dn",o)) for o in range(num_orbitals)],Operator())
    for o in range(num_orbitals):
        dn = N(*mkind("up",o)) - N(*mkind("dn",o))
        QN["P_"+str(o)] = dn*dn
 
print_master("Constructing the solver...")

# Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct.items())

print_master("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = range(num_orbitals), beta=beta, n_points=n_iw)
delta_w_part = delta_w.copy()
for e, v in zip(epsilon,V):
    delta_w_part <<= inverse(iOmega_n - e)
    delta_w_part.from_L_G_R(np.transpose(v),delta_w_part,v)
    delta_w += delta_w_part

for spin in spin_names:
    S.G0[spin] <<= inverse(iOmega_n + mu - delta_w)

print_master("Running the simulation...")

# Solve the problem
p['H_local'] = H
p['quantum_numbers'] = QN
p['legendre_accumulation'] = False
p['time_accumulation'] = True

start_time = datetime.now()
S.solve(**p)
print_master("Simulation lasted: " + str((datetime.now() - start_time).total_seconds()) + " seconds")

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(matrix_results_file_name,'w')
    for b in gf_struct: Results[b] = S.G_tau[b]

