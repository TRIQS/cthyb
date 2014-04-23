import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters import Parameters
#from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.applications.impurity_solvers.cthyb_segment import *
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
mu = 35.0
# Half-filling chemical potential
#mu = (U/2.0)*((2.0 * num_orbitals) - 1.0) - (5.0*J/2.0)*(num_orbitals - 1.0)

length_cycle = 50
n_warmup_cycles = 1000
n_cycles = 5000
n_tau_delta = 5000
n_tau_g = 5000

names = ['up-0','down-0','up-1','down-1','up-2','down-2']

Uintra=[[0,U],[U,0]] # intra orbital
Uinter=[[U-3*J,U-2*J],[U-2*J,U-3*J]] # inter orbital, inc Hund's
Umat = np.zeros((10,10))
for i in range(0,num_orbitals):
  Umat[2*i:2*i+2,2*i:2*i+2] = Uintra # diagonal
  for j in range(0,5):
    if i!=j:
      Umat[2*j:2*j+2,2*i:2*i+2] = Uinter # off diagonal

# Construct the solver
S = Solver(beta = beta, block_names = names, n_tau_delta = n_tau_delta, n_tau_g = n_tau_g)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)

for name, g0block in S.G0:
  g0block <<= inverse( iOmega_n + mu - delta_w )
#  g0block <<= inverse( iOmega_n + mu - (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth) )

n_loops=1
# Now do the DMFT loop
for IterNum in range(n_loops):

  starttime = time.clock()

  if IterNum > 0:
    # Compute S.Delta_tau with the self-consistency condition while imposing
    # paramagnetism
    g_w = GfImFreq(indices = [0], beta=beta)
    for o in range(0,num_orbitals):
      g_w <<= g_w + 0.5 * (1.0/num_orbitals) * Fourier(S.G_tau["up-%s"%o] + S.G_tau["down-%s"%o])

# first three moments (w, const, 1/w), number of moments, range of freq to fit
    g_w.fit_tail([[[0.0,0.0,1.0]]],5,100,500)

    for name, g0block in S.G0:
      g0block <<= inverse( iOmega_n + mu - (half_bandwidth/2.0)**2 * g_w )

  S.G0 = mpi.bcast(S.G0)

  S.solve(U = Umat, n_cycles = n_cycles, length_cycle = length_cycle, n_warmup_cycles = n_warmup_cycles)

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
  try :
      Results.create_group("log")
  except :
      pass 
  log = Results["log"]
  log["code_version"] = version.release
  log["script"] = open(sys.argv[0]).read() # read myself !
