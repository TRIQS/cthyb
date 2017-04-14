import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.applications.impurity_solvers.cthyb.util import estimate_nfft_buf_size
from pytriqs.gf.local import *
import numpy as np

# Input parameters
beta = 10.0
num_orb = 2
mu = 1.0
U = 2.0
J = 0.2
epsilon = [-2.3,2.3]
V = [1.0*np.eye(num_orb) + 0.1*(np.ones((num_orb,num_orb)) - np.eye(num_orb))]*2

spin_names = ("up","dn")
orb_names = range(num_orb)
n_iw = 1024

g2_n_iw = 5
g2_n_inu = 10
g2_n_l = 4
g2_blocks = set([("up","up"),("up","dn"),("dn","up")])

p = {}
p["verbosity"] = 2
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 500
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = False
p["measure_g_tau"] = False
p["measure_g_l"] = False

# Parameters for the preliminary run
p_pre = p.copy()
p_pre["measure_pert_order"] = True

p["measure_g2_inu"] = True
p["measure_g2_legendre"] = True
p["measure_g2_pp"] = True
p["measure_g2_ph"] = True
p["measure_g2_block_order"] = "AABB"
p["measure_g2_blocks"] = g2_blocks
p["measure_g2_n_iw"] = g2_n_iw
p["measure_g2_n_inu"] = g2_n_inu
p["measure_g2_n_l"] = g2_n_l

mpi.report("Welcome to measure_g2 benchmark.")

gf_struct = set_operator_structure(spin_names,orb_names,True)
mkind = get_mkind(True,None)

# Hamiltonian
H = h_int_kanamori(spin_names,orb_names,
                   np.array([[0,U-3*J],[U-3*J,0]]),
                   np.array([[U,U-2*J],[U-2*J,U]]),
                   J,True)

mpi.report("Constructing the solver...")

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=n_iw)

mpi.report("Preparing the hybridization function...")

# Set hybridization function
delta_w = GfImFreq(indices = orb_names, beta=beta, n_points=n_iw)
delta_w_part = delta_w.copy()
for e, v in zip(epsilon,V):
    delta_w_part << inverse(iOmega_n - e)
    delta_w_part.from_L_G_R(np.transpose(v),delta_w_part,v)
    delta_w += delta_w_part

S.G0_iw << inverse(iOmega_n + mu - delta_w)

# Solve the problem
mpi.report("Preliminary run...")
S.solve(h_int=H, **p_pre)

p["nfft_buf_sizes"] = estimate_nfft_buf_size(gf_struct, S.perturbation_order)

mpi.report("Running the simulation...")
S.solve(h_int=H, **p)

# Check shapes of g2 containers
ref_shape_inu = (2*g2_n_iw-1,2*g2_n_inu,2*g2_n_inu,
                 num_orb,num_orb,num_orb,num_orb)
ref_shape_l = (2*g2_n_iw-1,g2_n_l,g2_n_l,
               num_orb,num_orb,num_orb,num_orb)
for bn in g2_blocks:
    assert S.G2_iw_inu_inup_pp[bn].data.shape == ref_shape_inu
    assert S.G2_iw_inu_inup_ph[bn].data.shape == ref_shape_inu
    assert S.G2_iw_l_lp_pp[bn].data.shape == ref_shape_l
    assert S.G2_iw_l_lp_ph[bn].data.shape == ref_shape_l

# Save the results
if mpi.is_master_node():
    with HDFArchive("measure_g2.h5",'w') as ar:
        ar['beta'] = beta
        ar['U'] = U
        ar['mu'] = mu
        ar['J'] = J
        ar['epsilon'] = epsilon
        ar['V'] = V
        p['measure_g2_blocks'] = list(p['measure_g2_blocks'])
        ar['solve_params'] = p
        ar['G2_iw_inu_inup_pp'] = S.G2_iw_inu_inup_pp
        ar['G2_iw_inu_inup_ph'] = S.G2_iw_inu_inup_ph
        ar['G2_iw_l_lp_pp'] = S.G2_iw_l_lp_pp
        ar['G2_iw_l_lp_ph'] = S.G2_iw_l_lp_ph
