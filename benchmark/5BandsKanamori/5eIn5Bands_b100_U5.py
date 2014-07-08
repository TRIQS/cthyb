import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators.operators2 import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from collections import OrderedDict
import itertools
import time
import os

fileName=os.path.splitext(os.path.basename(__file__))[0]

beta = 100.0
# H_loc parameters
num_orbitals = 5
U = 5.0
J = 0.1
half_bandwidth = 1.0
#mu = 35.0
# Half-filling chemical potential
mu = (U/2.0)*((2.0 * num_orbitals) - 1.0) - (5.0*J/2.0)*(num_orbitals - 1.0)

# Parameters
p = SolverCore.solve_parameters()
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 10000
p["n_cycles"] = 100000
p["make_histograms"] = True
p["measure_g_tau"] = False
p["use_trace_estimator"] = False

# Block structure of GF
gf_struct = OrderedDict()
for o in range(0,num_orbitals):
  gf_struct['up-%s'%o] = [0,]
  gf_struct['down-%s'%o] = [0,]

# Hamiltonian -- must include both quadratic and quartic terms
H = Operator()
for o in range(0,num_orbitals):
    H += -mu*(n("up-%s"%o,0) + n("down-%s"%o,0))

for o in range(0,num_orbitals):
    H += U*n("up-%s"%o,0)*n("down-%s"%o,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*n("up-%s"%o1,0)*n("down-%s"%o2,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o2>=o1: continue;
    H += (U-3*J)*n("up-%s"%o1,0)*n("up-%s"%o2,0)
    H += (U-3*J)*n("down-%s"%o1,0)*n("down-%s"%o2,0)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += -J*c_dag("up-%s"%o1,0)*c_dag("down-%s"%o1,0)*c("up-%s"%o2,0)*c("down-%s"%o2,0)
    H += -J*c_dag("up-%s"%o1,0)*c_dag("down-%s"%o2,0)*c("up-%s"%o2,0)*c("down-%s"%o1,0)

# Quantum numbers
qn = []
for o in range(0,num_orbitals+2): qn.append(Operator())
for o in range(0,num_orbitals):
    qn[0] += n("up-%s"%o,0) + n("down-%s"%o,0) # Ntot
    qn[1] += (n("up-%s"%o,0) - n("down-%s"%o,0)) # Sz
    qn[2+o] += (n("up-%s"%o,0) - n("down-%s"%o,0))*(n("up-%s"%o,0) - n("down-%s"%o,0)) # Seniority number = Sz^2

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=10001, n_iw=1025)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta, n_points=1025)
delta_w <<= (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)

for o in range(0,num_orbitals):
    S.G0_iw["up-%s"%o] <<= inverse(iOmega_n - delta_w)
    S.G0_iw["down-%s"%o] <<= inverse(iOmega_n - delta_w)

n_loops=3
# Now do the DMFT loop
for IterNum in range(n_loops):

  starttime = time.clock()

  g_tau = S.G_tau.copy()
  g_w = GfImFreq(indices = [0], beta=beta, n_points=1025)

  if IterNum > 0:
    # Compute S.Delta_tau with the self-consistency condition while imposing
    # paramagnetism
    for o in range(0,num_orbitals):
      g_w <<= g_w + 0.5 * (1.0/num_orbitals) * Fourier(S.G_tau["up-%s"%o] + S.G_tau["down-%s"%o])

# first three moments (w, const, 1/w), number of moments, range of freq to fit
    g_w.fit_tail([[[0.0,0.0,1.0]]],5,100,500)

    for name, g0block in S.G0_iw:
      g0block <<= inverse(iOmega_n - (half_bandwidth/2.0)**2 * g_w )

  S.solve(h_loc=H, params=p)

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
  Results.create_group("log")
  log = Results["log"]
  log["code_version"] = version.release
  log["script"] = open(sys.argv[0]).read() # read myself !
