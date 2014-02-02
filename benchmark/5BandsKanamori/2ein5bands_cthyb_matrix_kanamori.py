import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_matrix import *
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
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
# mu = (U/2.0)*((2.0 * num_orbitals) - 1.0) - (5.0*J/2.0)*(num_orbitals - 1.0)

# Block structure of GF
gf_struct=[]
for o in range(0,num_orbitals):
    gf_struct.append(('up-%s'%o, [0]))
    gf_struct.append(('down-%s'%o, [0]))
    
# Hamiltonian -- quartic terms only (quadratic terms are implicit in G0)
H = Operator()
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
    H += -J*Cdag("up-%s"%o1,0)*Cdag("down-%s"%o1,0)*C("up-%s"%o2,0)*C("down-%s"%o2,0)
    H += -J*Cdag("up-%s"%o1,0)*Cdag("down-%s"%o2,0)*C("up-%s"%o2,0)*C("down-%s"%o1,0)
    
# Quantum numbers
qn = []
for o in range(0,num_orbitals+2): qn.append(Operator())
for o in range(0,num_orbitals):
    qn[0] += N("up-%s"%o,0) + N("down-%s"%o,0) # Ntot
    qn[1] += (N("up-%s"%o,0) - N("down-%s"%o,0)) # Sz
    qn[2+o] += (N("up-%s"%o,0) - N("down-%s"%o,0))*(N("up-%s"%o,0) - N("down-%s"%o,0)) # Seniority number = Sz^2
qndict={}
for i, op in enumerate(qn): qndict["%s"%i] = op
    
# Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct)

# Initalize the Green's function to a semi circular
S.G <<= SemiCircular(half_bandwidth)

n_loops=1
# Now do the DMFT loop
for IterNum in range(n_loops):

  starttime = time.clock()

  g = S.G["up-1"].copy()
  g <<= 0.0

  # Compute S.G0 with the self-consistency condition while imposing
  # paramagnetism
  for o in range(0,num_orbitals):
    g <<= g + S.G["up-%s"%o] + S.G["down-%s"%o]
  g <<= 0.5 * (1.0/num_orbitals) * g
  for name, g0block in S.G0:
    g0block <<= inverse( iOmega_n + mu - (half_bandwidth/2.0)**2 * g )

  S.solve(H_local = H,
          quantum_numbers = qndict,
          n_cycles = 5000,
          length_cycle = 50,
          n_warmup_cycles = 1000,
          n_legendre = 200,
          random_name = "mt19937",
          use_segment_picture = False)

  endtime = time.clock()
  print( "Time for loop %s: "%IterNum, endtime - starttime, " seconds" )

  if mpi.rank==0:
      Results = HDFArchive(fileName+".h5",'a')
      Results["G-%s"%IterNum] = S.G
      Results["Delta-%s"%IterNum] = S.Delta_tau
      Results["GLegendre-%s"%IterNum] = S.G_legendre
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
