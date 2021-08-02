import numpy as np
from triqs.gf import *
from triqs.operators import *
from triqs.atom_diag import AtomDiag, atomic_g_tau
from triqs.utility.comparison_tests import *

# for the imaginary GFs
beta = 40
n_tau = 5001

H = np.array([[0,0.1j],[-0.1j,0]])

# construct the atomic GF by hand:
#  G_tau = inverse_fourier( (iOmega - H)^(-1) )
G_iw = GfImFreq(beta=beta,indices=[0,1])
G_iw << iOmega_n - H
G_iw.invert()
G_tau = GfImTime(beta=beta,indices=[0,1],n_points=n_tau)
G_tau.set_from_fourier(G_iw)

# construct the Hamiltonian as operator
H_op = Operator()

for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        H_op += H[i,j] * c_dag('ud',i)*c('ud',j)

# initialize the atomic diagonalization
AD = AtomDiag(H_op,[('ud',0),('ud',1)])

# atomic G(tau) from solver
G_at = atomic_g_tau(AD,beta,[['ud', 2]],n_tau)

assert_gfs_are_close(G_at['ud'],G_tau)

