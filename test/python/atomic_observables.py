import pytriqs.utility.mpi as mpi
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.U_matrix import U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from pytriqs.operators.util.observables import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *

beta = 100.0
# H_int parameters
L = 1 # angular momentum
F0 = 4.0
F2 = 0.5

spin_names = ("up","dn")
cubic_names = map(str,range(2*L+1))
U_mat = U_matrix(L, radial_integrals=[F0,F2], basis="spherical")

# Parameters
p = {}
p["verbosity"] = 0
p["length_cycle"] = 50
p["n_warmup_cycles"] = 0
p["n_cycles"] = 0

# Block structure of GF
gf_struct = set_operator_structure(spin_names,cubic_names,False)

# Local Hamiltonian
H = h_int_slater(spin_names,cubic_names,U_mat,False)

# Observables
N = N_op(spin_names,cubic_names,False)
S2 = S2_op(spin_names,cubic_names,False)
Sz = S_op('z',spin_names,cubic_names,False)
L2 = L2_op(spin_names,cubic_names,False)
Lz = L_op('z',spin_names,cubic_names,False)
LS = LS_op(spin_names,cubic_names,False)

# Additional splitting terms to lift the degeneracies
H += 0.22*Sz
H += 0.33*Lz

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=10000)

# Set G0(iw)
S.G0_iw << inverse(iOmega_n)

S.solve(h_int=H, **p)

obs = {'E':H,'N':N,'S2':S2,'Sz':Sz,'L2':L2,'Lz':Lz} #,'LS':LS}
res = dict ( (name, [item for v in quantum_number_eigenvalues(op,S.h_loc_diagonalization) for item in v ]) for (name,op) in obs.items())

print ("%9s "*7) % ("Energy","N","S^2","S_z","L^2","L_z","L*S")
print ("%9s "*7) % ("======","=","===","===","===","===","===")

E_sorted_res = sorted(zip(res['E'],res['N'],res['S2'],res['Sz'],res['L2'],res['Lz']),
#E_sorted_res = sorted(zip(res['E'],res['N'],res['S2'],res['Sz'],res['L2'],res['Lz'],res['LS']),
                      cmp=lambda r1,r2: cmp(r1[0],r2[0]))

filter_minus_0 = lambda x: 0 if (x<0 and abs(x)<1e-10) else x
for v in E_sorted_res: print ("%9.4f "*6) % tuple(map(filter_minus_0,v))
