#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import n
from pytriqs.atom_diag import trace_rho_op
from triqs_cthyb import SolverCore
from pytriqs.gf import GfImFreq, iOmega_n, inverse
import numpy as np

spin_names = ("up","dn")
mkind = lambda sn, on: (sn,on)
gf_struct = [ ["up", [1, 2]], ["dn", [1, 2]] ]

# Input parameters
beta = 100.0
U = 2.0
mu = 1.0
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
p["n_warmup_cycles"] = 100000
p["n_cycles"] = 1000000
p["measure_G_tau"] = True
p["measure_density_matrix"] = True
p["use_norm_as_weight"] = True

move_global_prob = 0.05

H = U*n(*mkind("up",1))*n(*mkind("dn",1)) + U*n(*mkind("up",2))*n(*mkind("dn",2))

# Global moves
gm_flip_spins_1   = {'flip_spins_1' :   {mkind("up",1) : mkind("dn",1), mkind("dn",1) : mkind("up",1)}}
gm_flip_spins_all = {'flip_spins_all' : {mkind("up",1) : mkind("dn",1), mkind("dn",1) : mkind("up",1),
                                         mkind("up",2) : mkind("dn",2), mkind("dn",2) : mkind("up",2)}}
gm_swap_atoms     = {'swap_atoms' :     {mkind("up",1) : mkind("up",2), mkind("dn",1) : mkind("dn",2),
                                         mkind("up",2) : mkind("up",1), mkind("dn",2) : mkind("dn",1)}}

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

# Set hybridization function
delta_w = GfImFreq(indices = [1,2], beta=beta)
delta_w << (V**2)*(inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon))
for sn in spin_names:
    S.G0_iw[sn] << inverse(iOmega_n - np.matrix([[-mu,t],[t,-mu]]) - delta_w)

if mpi.is_master_node():
    arch = HDFArchive("move_global_beta%.0f_prob%.2f.h5" %(beta,move_global_prob) ,'w')
    arch['beta'] = beta
    arch['move_global_prob'] = move_global_prob

static_observables = {"N1_up" : n(*mkind("up",1)), "N1_dn" : n(*mkind("dn",1)),
                      "N2_up" : n(*mkind("up",2)), "N2_dn" : n(*mkind("dn",2))}

global_moves = [('none',{}),
                ('flip_spins_1',gm_flip_spins_1),
                ('flip_spins_all',gm_flip_spins_all),
                ('swap_atoms',gm_swap_atoms),
                ('swap_and_flip',dict(gm_flip_spins_all,**gm_swap_atoms))]

for gm_name, gm in global_moves:
    mpi.report("Running with global moves set '%s'" % gm_name)
    if gm_name != 'none':
        p["move_global"] = gm
        p["move_global_prob"] = move_global_prob

    # Solve the problem
    S.solve(h_int=H, **p)

    if mpi.is_master_node():
        # Save the results
        arch.create_group(gm_name)
        arch[gm_name]['G_tau'] = S.G_tau

        dm = S.density_matrix
        for on, o in static_observables.items():
            arch[gm_name][on] = trace_rho_op(dm,o,S.h_loc_diagonalization)

