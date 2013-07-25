@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.utility import *
from pytriqs.utility.mpi import rank, is_master_node
import numpy as np

# Set up a few parameters
beta = 10.0
n_orb = 2
U = 2.0
J = 0.2
V = 1.0*np.identity(n_orb) + 0.1*(np.ones((n_orb,n_orb))-np.identity(n_orb))
mu = 1.0
epsilon = 2.3

# Construct a CTQMC solver
from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver
struct = [ ('up',range(n_orb)), ('down',range(n_orb)) ]
S = Solver(beta = beta, gf_struct = struct)

# Compute Delta and G0
Delta = S.G0.copy()
for n,d in Delta:
  d <<= inverse(iOmega_n + epsilon) + inverse(iOmega_n - epsilon)
  d.from_L_G_R(V, d, V)
  S.G0[n] <<= inverse(iOmega_n + mu - d)

# Kanamori Hamiltonian
from pytriqs.applications.impurity_solvers.operators import *
H = Operator()

for o in range(n_orb):
  H += U * N("up",o)*N("down",o)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o1==o2): continue
    H += (U-2*J)*N("up",o1)*N("down",o2)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o2>=o1): continue
    H += (U-3*J)*N("up",o1)*N("up",o2)
    H += (U-3*J)*N("down",o1)*N("down",o2)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o1==o2): continue
    H += -J*Cdag("up",o1)*Cdag("down",o1)*C("up",o2)*C("down",o2)
    H += -J*Cdag("up",o1)*Cdag("down",o2)*C("up",o2)*C("down",o1)

# Quantum numbers
N_tot = reduce(lambda x,y: x+y, [N("up",o)+N("down",o) for o in range(n_orb)])
Sz_tot = reduce(lambda x,y: x+y, [N("up",o)-N("down",o) for o in range(n_orb)])

# Solve
S.solve(H_local = H,
        quantum_numbers = { 'Ntot' : N_tot, 'Sz_tot' : Sz_tot },
        n_cycles  = 10000000,
        length_cycle = 50,
        n_warmup_cycles = 50000,
        n_time_slices_delta = 1000,
        n_time_slices_gtau = 1000,
        legendre_accumulation = False,
        time_accumulation = True,
        random_name = "",
	random_seed = 123*rank + 567,
        use_segment_picture = False)

# Calculation is done. Now save a few things
if is_master_node():
	Results = HDFArchive("kanamori_offdiag.hyb.h5",'w')
	Results["G"] = S.G
	Results["G_tau"] = S.G_tau

