import numpy as np
from pytriqs.atom_diag import AtomDiag, atomic_g_tau
from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.utility.comparison_tests import *
from pytriqs.archive import HDFArchive


# the single-particle Hamiltonian of the problem
H_matrix = [[-0.2+0.j  , 0.0+0.1j , 0.5+0.j ,  0.1+0.j ],
            [ 0.0-0.1j ,-0.2+0.j  , 0.1+0.j ,  0.5+0.j ],
            [ 0.5+0.j  , 0.1+0.j  , 0.0+0.j ,  0.0+0.j ],
            [ 0.1+0.j  , 0.5+0.j  , 0.0+0.j ,  0.0+0.j ]]
H_matrix = np.array(H_matrix)

# fops for the AtomDiag solver
fops = []
for i in range(H_matrix.shape[0]):
    fops.append(('ud',i))

H = Operator()

for i in range(H_matrix.shape[0]):
    for j in range(H_matrix.shape[1]):
        H += H_matrix[i,j] * c_dag('ud',i)*c('ud',j)

# interaction term between the 1st and 2nd orbital
H+=3*n('ud',0)*n('ud',1)

S = AtomDiag(H,fops)
G_tau = atomic_g_tau(S,beta=10,gf_struct=[['ud', range(4)]],n_tau=41)

with HDFArchive('atomdiag_ed.out.h5','w') as ar:
    ar['G_tau']=G_tau

from pytriqs.utility.h5diff import h5diff
# the reference solution was checked against an external ED result
h5diff('atomdiag_ed.ref.h5','atomdiag_ed.out.h5')


