import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.operators import *
#from atom_diag import trace_rho_op
from triqs.atom_diag import trace_rho_op
from triqs_cthyb import *
from triqs.gf import *
import numpy as np

def is_hermitian(matrix):
    N, M = matrix.shape
    return any([matrix[i,j] == matrix[j,i].conjugate() for i in range(N) for j in range(M)])

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
h = 0.1
V = 1.0
t = 0.1
epsilon = 2.3

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 10000
p["measure_G_tau"] = True
p["move_double"] = False
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = True

H = U*n("up",0)*n("dn",0) + U*n("up",1)*n("dn",1)
H = H + 0.5*h*(n("up",0) - n("dn",0)) + 0.5*h*(n("up",1) - n("dn",1))

# Construct the solver
S = Solver(beta=beta, gf_struct=[["dn",2], ["up",2]], n_tau=n_tau, n_iw=n_iw)

# Set hybridization function
delta_w = GfImFreq(beta=beta, target_shape=(2,2))
delta_w << (V**2)*(inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon))
for bn, g in S.G0_iw: g << inverse(iOmega_n - np.matrix([[-mu,t],[t,-mu]]) - delta_w)

# Solve!
S.solve(h_int=H, **p)


if mpi.is_master_node():
    # Measure expectation values
    dm = S.density_matrix
    assert all([is_hermitian(b) for b in dm]), "density matrix is not hermitian"
    
    G_moments = S.G_moments
    assert all([is_hermitian(g) for name in G_moments for g in G_moments[name]]), "Gf moments are not hermitian!"

    Sigma_moments = S.Sigma_moments
    assert all([is_hermitian(g) for name in G_moments for g in Sigma_moments[name]]), "Sigma moments are not hermitian!"
