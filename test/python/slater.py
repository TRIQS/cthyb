from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators.operators2 import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from U_matrix import U_matrix

beta = 100.0
# H_loc parameters
num_orbitals = 5
U = 5.0
J = 0.1
F0 = U
F2 = J*(14.0/(1.0 + 0.63)) # Sasha uses coefficient 0.63 instead of 0.625
F4 = F2*0.63
L = 2 # angular momentum
half_bandwidth = 1.0
mu = 32.5  # 3 electrons in 5 bands

mkind = lambda sn, cn: (sn+'-'+cn, 0)
spin_names = ("up","down")
cubic_names=['%s'%i for i in range(num_orbitals)]
U_matrix = U_matrix(L, radial_integrals=[F0,F2,F4], basis="cubic")

# Parameters
p = SolverCore.solve_parameters()
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50
p["n_cycles"] = 5000
p["measure_g_l"] = True

# Block structure of GF
gf_struct = {}
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    gf_struct[bn] = [i]

# Local Hamiltonian
H = Operator()

a_range=range(num_orbitals)
for s1, s2 in product(spin_names,spin_names):
    for ap1, ap2, a1, a2 in product(a_range,a_range,a_range,a_range):
        U_val = U_matrix[ap1,ap2,a1,a2]
        if abs(U_val.imag) > 1e-10:
            raise RuntimeError("Cubic harmonics are real, so should be the matrix elements of U.")

        H_term = 0.5*U_val.real*c_dag(*mkind(s1,cubic_names[ap1]))*c_dag(*mkind(s2,cubic_names[ap2]))*c(*mkind(s2,cubic_names[a1]))*c(*mkind(s1,cubic_names[a2]))
        H += H_term

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=100000)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)
for name, g0 in S.G0_iw:
    g0 <<= inverse(iOmega_n + mu - delta_w)

S.solve(h_loc=H, params=p)

if mpi.rank==0:
    Results = HDFArchive("slater.output.h5",'w')
    Results["G_tau"] = S.G_tau
    Results["G_leg"] = S.G_l
