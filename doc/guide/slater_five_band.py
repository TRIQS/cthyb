from triqs.gf import *
from h5 import HDFArchive
from triqs.applications.impurity_solvers.cthyb import Solver
import triqs.operators.util as op
import triqs.utility.mpi as mpi

# General parameters
filename = 'slater_five_band.h5'             # Name of file to save data to
beta = 40.0                                  # Inverse temperature
l = 2                                        # Angular momentum
n_orbs = 2*l + 1                             # Number of orbitals
U = 4.0                                      # Screened Coulomb interaction
J = 1.0                                      # Hund's coupling
half_bandwidth = 1.0                         # Half bandwidth
mu = 25.0                                    # Chemical potential
spin_names = ['up','down']                   # Outer (non-hybridizing) blocks
orb_names = ['%s'%i for i in range(n_orbs)]  # Orbital indices
off_diag = False                             # Include orbital off-diagonal elements?
n_loop = 2                                   # Number of DMFT self-consistency loops

# Solver parameters
p = {}
p["n_warmup_cycles"] = 10000     # Number of warmup cycles to equilibrate
p["n_cycles"] = 100000           # Number of measurements
p["length_cycle"] = 500          # Number of MC steps between consecutive measures
p["move_double"] = True          # Use four-operator moves

# Block structure of Green's functions
# gf_struct = [ ('up',5), ('down',5) ]
# This can be computed using the TRIQS function as follows:
gf_struct = op.set_operator_structure(spin_names,orb_names,off_diag=off_diag) 

# Construct the 4-index U matrix U_{ijkl}
# The spherically-symmetric U matrix is parametrised by the radial integrals
# F_0, F_2, F_4, which are related to U and J. We use the functions provided
# in the TRIQS library to construct this easily:
U_mat = op.U_matrix(l=l, U_int=U, J_hund=J, basis='spherical')

# Set the interacting part of the local Hamiltonian
# Here we use the full rotationally-invariant interaction parametrised 
# by the 4-index tensor U_{ijkl}.
# The TRIQS library provides a function to build this Hamiltonian from the U tensor:
H = op.h_int_slater(spin_names,orb_names,U_mat,off_diag=off_diag)

# Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct)

# Set the hybridization function and G0_iw for the Bethe lattice
delta_iw = GfImFreq(indices=[0], beta=beta)
delta_iw << (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)
for name, g0 in S.G0_iw: g0 << inverse(iOmega_n + mu - delta_iw)

# Now start the DMFT loops
for i_loop in range(n_loop):

  # Determine the new Weiss field G0_iw self-consistently
  if i_loop > 0: 
    g_iw = GfImFreq(indices=[0], beta=beta)
    # Impose paramagnetism
    for name, g in S.G_tau: g_iw << g_iw + (0.5/n_orbs)*Fourier(g)
    # Compute S.G0_iw with the self-consistency condition
    for name, g0 in S.G0_iw: 
        g0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * g_iw )

  # Solve the impurity problem for the given interacting Hamiltonian and set of parameters
  S.solve(h_int=H, **p)

  # Save quantities of interest on the master node to an h5 archive
  if mpi.is_master_node():
      with HDFArchive(filename,'a') as Results:
          Results['G0_iw-%s'%i_loop] = S.G0_iw
          Results['G_tau-%s'%i_loop] = S.G_tau
