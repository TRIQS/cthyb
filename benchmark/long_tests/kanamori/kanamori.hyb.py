@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.utility import *
from pytriqs.utility.mpi import rank, is_master_node
import numpy

# Set up a few parameters
beta = 10.0
n_orb = 2
U = 2.0
J = 0.2
V = 1.0
mu = 1.0
epsilon = 2.3

# Construct a CTQMC solver
from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver
struct = [ ('up-%i'%o,[0]) for o in range(n_orb) ]
struct += [ ('down-%i'%o,[0]) for o in range(n_orb) ]
S = Solver(beta = beta, gf_struct = struct)

# Compute G0
for o in range(n_orb):
  S.G0['up-%i'%o] <<= inverse( iOmega_n + mu - V*V*inverse(iOmega_n - epsilon) - V*V*inverse(iOmega_n + epsilon) )
  S.G0['down-%i'%o] <<= inverse( iOmega_n + mu - V*V*inverse(iOmega_n - epsilon) - V*V*inverse(iOmega_n + epsilon) )

# Kanamori Hamiltonian
from pytriqs.applications.impurity_solvers.operators import *
H = Operator()

for o in range(n_orb):
  H += U * N("up-%i"%o,0)*N("down-%i"%o,0)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o1==o2): continue
    H += (U-2*J)*N("up-%i"%o1,0)*N("down-%i"%o2,0)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o2>=o1): continue
    H += (U-3*J)*N("up-%i"%o1,0)*N("up-%i"%o2,0)
    H += (U-3*J)*N("down-%i"%o1,0)*N("down-%i"%o2,0)

for o1 in range(n_orb):
  for o2 in range(n_orb):
    if (o1==o2): continue
    H += -J*Cdag("up-%i"%o1,0)*Cdag("down-%i"%o1,0)*C("up-%i"%o2,0)*C("down-%i"%o2,0)
    H += -J*Cdag("up-%i"%o1,0)*Cdag("down-%i"%o2,0)*C("up-%i"%o2,0)*C("down-%i"%o1,0)

# Quantum numbers
N_tot = reduce(lambda x,y: x+y, [N("up-%i"%o,0)+N("down-%i"%o,0) for o in range(n_orb)])
Sz_tot = reduce(lambda x,y: x+y, [N("up-%i"%o,0)-N("down-%i"%o,0) for o in range(n_orb)])

# Solve
S.solve(H_local = H,
        quantum_numbers = { 'Ntot' : N_tot, 'Sz_tot' : Sz_tot },
        n_cycles  = 3000000,
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
	Results = HDFArchive("kanamori.hyb.h5",'w')
	Results["G"] = S.G
	Results["G_tau"] = S.G_tau

