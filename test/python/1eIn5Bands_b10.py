import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_krylov import *
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from collections import OrderedDict
import itertools
import time
import os

fileName=os.path.splitext(os.path.basename(__file__))[0]

beta = 10.0
# H_loc parameters
num_orbitals = 5
U = 5.0
J = 0.1
half_bandwidth = 1.0
mu = 40.0
# Half-filling chemical potential
# mu = (U/2.0)*((2.0 * num_orbitals) - 1.0) - (5.0*J/2.0)*(num_orbitals - 1.0)

# Parameters
p = Parameters()
p["beta"] = beta
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 500
p["n_warmup_cycles"] = 10000
p["n_cycles"] = 100000
p["n_tau_delta"] = 100000
p["n_tau_g"] = 100000
p["krylov_bs_use_cutoff"] = True
p["krylov_bs_prob_cutoff"] = -1.0
p["krylov_small_matrix_size"] = 100
p["make_path_histograms"] = True
p["use_old_trace"] = False
p["use_truncation"] = True
p["measure_gt"] = True
p["trace_estimator"] = "None"

# Block structure of GF
gf_struct = OrderedDict()
for o in range(0,num_orbitals):
  gf_struct['up-%s'%o] = [0,]
  gf_struct['down-%s'%o] = [0,]

# Hamiltonian -- must include both quadratic and quartic terms
H = Operator()
for o in range(0,num_orbitals):
    H += -mu*(N("up-%s"%o,0) + N("down-%s"%o,0))

for o in range(0,num_orbitals):
    H += U*N("up-%s"%o,0)*N("down-%s"%o,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*N("up-%s"%o1,0)*N("down-%s"%o2,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o2>=o1: continue;
    H += (U-3*J)*N("up-%s"%o1,0)*N("up-%s"%o2,0)
    H += (U-3*J)*N("down-%s"%o1,0)*N("down-%s"%o2,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += -J*C_dag("up-%s"%o1,0)*C_dag("down-%s"%o1,0)*C("up-%s"%o2,0)*C("down-%s"%o2,0)
    H += -J*C_dag("up-%s"%o1,0)*C_dag("down-%s"%o2,0)*C("up-%s"%o2,0)*C("down-%s"%o1,0)

# Quantum numbers
qn = []
for o in range(0,num_orbitals+2): qn.append(Operator())
for o in range(0,num_orbitals):
    qn[0] += N("up-%s"%o,0) + N("down-%s"%o,0) # Ntot
    qn[1] += (N("up-%s"%o,0) - N("down-%s"%o,0)) # Sz
    qn[2+o] += (N("up-%s"%o,0) - N("down-%s"%o,0))*(N("up-%s"%o,0) - N("down-%s"%o,0)) # Seniority number = Sz^2

# Construct the solver
S = Solver(parameters=p, H_local=H, quantum_numbers=qn, gf_struct=gf_struct)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)

for o in range(0,num_orbitals):
    S.Delta_tau["up-%s"%o] <<= InverseFourier(delta_w)
    S.Delta_tau["down-%s"%o] <<= InverseFourier(delta_w)

n_loops=1
# Now do the DMFT loop
for IterNum in range(n_loops):

  starttime = time.clock()

  g_tau = S.G_tau.copy()
  g_w = GfImFreq(indices = [0], beta=beta)

  if IterNum > 0:
    # Compute S.Delta_tau with the self-consistency condition while imposing
    # paramagnetism
    for o in range(0,num_orbitals):
      g_w <<= g_w + 0.5 * (1.0/num_orbitals) * Fourier(S.G_tau["up-%s"%o] + S.G_tau["down-%s"%o])

# first three moments (w, const, 1/w), number of moments, range of freq to fit
    g_w.fit_tail([[[0.0,0.0,1.0]]],5,100,500)

    for name, d0block in S.Delta_tau:
      d0block <<= InverseFourier( (half_bandwidth/2.0)**2 * g_w )

  S.solve(parameters = p)

  endtime = time.clock()
  print( "Time for loop %s: "%IterNum, endtime - starttime, " seconds" )

  if mpi.rank==0:
      Results = HDFArchive(fileName+".h5",'a')
      Results["G-%s"%IterNum] = S.G_tau
      Results["Delta-%s"%IterNum] = S.Delta_tau
      Results["Time-%s"%IterNum] = endtime - starttime

# Save script and parameters in HDFArchive
import sys, pytriqs.version as version
if mpi.rank==0:
  Results = HDFArchive(fileName+".h5",'a')
  try:
    Results.create_group("log")
  except:
    pass
  log = Results["log"]
  log["code_version"] = version.release
  log["script"] = open(sys.argv[0]).read() # read myself !

