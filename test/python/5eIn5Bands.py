import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators.operators2 import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
import itertools
from collections import OrderedDict
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
p["verbosity"] = 2
p["length_cycle"] = 10
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 5000

# Block structure of GF
gf_struct = OrderedDict()
for o in range(0,num_orbitals):
  gf_struct['up-%s'%o] = [0,]
  gf_struct['down-%s'%o] = [0,]

# Construct solver    
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau_delta=1000, n_tau_g=1000)

# Hamiltonian -- must include both quadratic and quartic terms
H = c_dag("up-0",0) - c_dag("up-0",0)
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

# first three moments (w, const, 1/w), number of moments, range of freq to fit -- last 20 percent of (2n+1)*pi/beta
    g_w.fit_tail([[[0.0,0.0,1.0]]],5,1600*3.14/beta,2000*3.14/beta)

    for name, d0block in S.Delta_tau:
      d0block <<= InverseFourier( (half_bandwidth/2.0)**2 * g_w )

  S.solve(h_loc=H, params=p)

  endtime = time.clock()
  print( "Time for loop %s: "%IterNum, endtime - starttime, " seconds" )

  if mpi.rank==0:
      Results = HDFArchive(fileName+".h5",'a')
      Results["G-%s"%IterNum] = S.G_tau
      Results["Delta-%s"%IterNum] = S.Delta_tau
      Results["Time-%s"%IterNum] = endtime - starttime
      del Results

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
  del Results
