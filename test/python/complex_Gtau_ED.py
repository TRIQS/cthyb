import numpy as np
from pytriqs.operators import *
from triqs_cthyb import *
from pytriqs.gf import *
from h5 import HDFArchive
from pytriqs.utility.h5diff import h5diff

""" This test was benchmarked against the an ED-solver that is 
independent of TRIQS and was supplied by Wei Wu (May 2016). 
Due to noise in the Monte Carlo data, it is not possible to 
automatically check whether the CTHYB solver gives the same result
as the ED code.
Thus, this test just checks against the CTHYB data it produced when
the test was created, which were checked against the ED result.
The data from the ED solver can be found in the file 
complex_Gtau_ED.bench.h5 for future reference.
"""

# the Hamiltonian of the system
# the first two orbitals are the impurity, the other two are the bath
H_mat = np.array([[-0.2  , 0.1j , 0.5 ,  0.1 ],
                  [-0.1j ,-0.3  , 0.1 ,  0.5 ],
                  [ 0.5  , 0.1  , 0.1 ,  0.0 ],
                  [ 0.1  , 0.5  , 0.0 ,  0.0 ]])
corr_dim = 2

G0_iw = GfImFreq(beta=10,indices=list(range(len(H_mat))),n_points=101)
G0_iw << inverse(iOmega_n - H_mat)

H_int = 3*n("ud",0)*n("ud",1)

p = {}
p["random_seed"] = 123 
p["length_cycle"] = 100
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 5000

S = Solver(beta=10,gf_struct=[["ud",list(range(corr_dim))]],n_tau=203,n_iw=101)
S.G0_iw << G0_iw[:corr_dim,:corr_dim]
S.solve(h_int=H_int,**p)

with HDFArchive("complex_Gtau_ED.out.h5","w") as ar:
    ar["G_tau"]=S.G_tau

h5diff("complex_Gtau_ED.ref.h5","complex_Gtau_ED.out.h5")

