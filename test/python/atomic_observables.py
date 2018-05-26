import pytriqs.utility.mpi as mpi
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.U_matrix import U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from pytriqs.operators.util.observables import *
from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.atom_diag import quantum_number_eigenvalues
#from atom_diag import * # quantum_number_eigenvalues
from cthyb import *

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

OUT = ("%9s "*7) % ("Energy","N","S^2","S_z","L^2","L_z","L*S")  + '\n'
OUT += ("%9s "*7) % ("======","=","===","===","===","===","===")  + '\n'

E_sorted_res = sorted(zip(res['E'],res['N'],res['S2'],res['Sz'],res['L2'],res['Lz']),
#E_sorted_res = sorted(zip(res['E'],res['N'],res['S2'],res['Sz'],res['L2'],res['Lz'],res['LS']),
                      cmp=lambda r1,r2: cmp(r1[0],r2[0]))

filter_minus_0 = lambda x: 0 if (x<0 and abs(x)<1e-10) else x
for v in E_sorted_res: OUT += ("%9.4f "*6) % tuple(map(filter_minus_0,v)) + '\n'

ref = """   Energy         N       S^2       S_z       L^2       L_z       L*S 
   ======         =       ===       ===       ===       ===       === 
  -0.4400    1.0000    0.7500   -0.5000    2.0000   -1.0000 
  -0.2200    1.0000    0.7500    0.5000    2.0000   -1.0000 
  -0.1100    1.0000    0.7500   -0.5000    2.0000    0.0000 
   0.0000    0.0000    0.0000    0.0000    0.0000    0.0000 
   0.1100    1.0000    0.7500    0.5000    2.0000    0.0000 
   0.2200    1.0000    0.7500   -0.5000    2.0000    1.0000 
   0.4400    1.0000    0.7500    0.5000    2.0000    1.0000 
   3.3500    2.0000    2.0000   -1.0000    2.0000   -1.0000 
   3.3600    2.0000    0.0000    0.0000    6.0000   -2.0000 
   3.5700    2.0000    2.0000    0.0000    2.0000   -1.0000 
   3.6800    2.0000    2.0000   -1.0000    2.0000    0.0000 
   3.6900    2.0000    0.0000    0.0000    6.0000   -1.0000 
   3.7900    2.0000    2.0000    1.0000    2.0000   -1.0000 
   3.9000    2.0000    2.0000    0.0000    2.0000    0.0000 
   4.0100    2.0000    2.0000   -1.0000    2.0000    1.0000 
   4.0200    2.0000    0.0000    0.0000    6.0000    0.0000 
   4.1200    2.0000    2.0000    1.0000    2.0000    0.0000 
   4.2000    2.0000    0.0000    0.0000    0.0000    0.0000 
   4.2300    2.0000    2.0000    0.0000    2.0000    1.0000 
   4.3500    2.0000    0.0000    0.0000    6.0000    1.0000 
   4.4500    2.0000    2.0000    1.0000    2.0000    1.0000 
   4.6800    2.0000    0.0000    0.0000    6.0000    2.0000 
  11.1100    3.0000    0.7500   -0.5000    6.0000   -2.0000 
  11.3300    3.0000    0.7500    0.5000    6.0000   -2.0000 
  11.3700    3.0000    3.7500   -1.5000    0.0000    0.0000 
  11.4400    3.0000    0.7500   -0.5000    6.0000   -1.0000 
  11.5600    3.0000    0.7500   -0.5000    2.0000   -1.0000 
  11.5900    3.0000    3.7500   -0.5000    0.0000    0.0000 
  11.6600    3.0000    0.7500    0.5000    6.0000   -1.0000 
  11.7700    3.0000    0.7500   -0.5000    6.0000    0.0000 
  11.7800    3.0000    0.7500    0.5000    2.0000   -1.0000 
  11.8100    3.0000    3.7500    0.5000    0.0000    0.0000 
  11.8900    3.0000    0.7500   -0.5000    2.0000    0.0000 
  11.9900    3.0000    0.7500    0.5000    6.0000    0.0000 
  12.0300    3.0000    3.7500    1.5000    0.0000    0.0000 
  12.1000    3.0000    0.7500   -0.5000    6.0000    1.0000 
  12.1100    3.0000    0.7500    0.5000    2.0000    0.0000 
  12.2200    3.0000    0.7500   -0.5000    2.0000    1.0000 
  12.3200    3.0000    0.7500    0.5000    6.0000    1.0000 
  12.4300    3.0000    0.7500   -0.5000    6.0000    2.0000 
  12.4400    3.0000    0.7500    0.5000    2.0000    1.0000 
  12.6500    3.0000    0.7500    0.5000    6.0000    2.0000 
  23.1500    4.0000    2.0000   -1.0000    2.0000   -1.0000 
  23.1600    4.0000    0.0000    0.0000    6.0000   -2.0000 
  23.3700    4.0000    2.0000    0.0000    2.0000   -1.0000 
  23.4800    4.0000    2.0000   -1.0000    2.0000    0.0000 
  23.4900    4.0000    0.0000    0.0000    6.0000   -1.0000 
  23.5900    4.0000    2.0000    1.0000    2.0000   -1.0000 
  23.7000    4.0000    2.0000    0.0000    2.0000    0.0000 
  23.8100    4.0000    2.0000   -1.0000    2.0000    1.0000 
  23.8200    4.0000    0.0000    0.0000    6.0000    0.0000 
  23.9200    4.0000    2.0000    1.0000    2.0000    0.0000 
  24.0000    4.0000    0.0000    0.0000    0.0000    0.0000 
  24.0300    4.0000    2.0000    0.0000    2.0000    1.0000 
  24.1500    4.0000    0.0000    0.0000    6.0000    1.0000 
  24.2500    4.0000    2.0000    1.0000    2.0000    1.0000 
  24.4800    4.0000    0.0000    0.0000    6.0000    2.0000 
  39.1600    5.0000    0.7500   -0.5000    2.0000   -1.0000 
  39.3800    5.0000    0.7500    0.5000    2.0000   -1.0000 
  39.4900    5.0000    0.7500   -0.5000    2.0000    0.0000 
  39.7100    5.0000    0.7500    0.5000    2.0000    0.0000 
  39.8200    5.0000    0.7500   -0.5000    2.0000    1.0000 
  40.0400    5.0000    0.7500    0.5000    2.0000    1.0000 
  59.4000    6.0000    0.0000    0.0000    0.0000    0.0000 
"""

assert OUT == ref




