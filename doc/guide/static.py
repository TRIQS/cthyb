from triqs.gf import *
from triqs.operators import *
from triqs.applications.impurity_solvers.cthyb import Solver

# Parameters
D = 1.0         # Half-bandwidth of the bath
V = 0.2         # Hybridisation amplitude
U = 4.0         # Coulomb interaction
e_f = -U/2      # Local energy level
h = 0.01        # External field
beta = 50       # Inverse temperature

# Construct the impurity solver with the inverse temperature
# and the structure of the Green's functions
S = Solver(beta = beta, gf_struct = [ ('up',1), ('down',1) ])

# Initialize the non-interacting Green's function S.G0_iw
# External magnetic field introduces Zeeman energy splitting between
# different spin components
S.G0_iw['up']   << inverse(iOmega_n - e_f + h/2 - V**2 * Wilson(D))
S.G0_iw['down'] << inverse(iOmega_n - e_f - h/2 - V**2 * Wilson(D))

# Run the solver
S.solve(h_int = U * n('up',0) * n('down',0),     # Local Hamiltonian
        n_cycles  = 500000,                      # Number of QMC cycles
        length_cycle = 200,                      # Length of one cycle
        n_warmup_cycles = 10000,                 # Warmup cycles
        measure_density_matrix = True,           # Measure the reduced density matrix
        use_norm_as_weight = True)               # Required to measure the density matrix

# Extract accumulated density matrix
rho = S.density_matrix

# Object containing eigensystem of the local Hamiltonian
h_loc_diag = S.h_loc_diagonalization

from triqs.atom_diag import trace_rho_op

# Evaluate occupations
print("<N_up> =", trace_rho_op(rho, n('up',0), h_loc_diag))
print("<N_down> = ", trace_rho_op(rho, n('down',0), h_loc_diag))

# Evaluate double occupancy
print("<N_up*N_down> =", trace_rho_op(rho, n('up',0)*n('down',0), h_loc_diag))
