#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from h5 import HDFArchive
from pytriqs.operators import n, Operator
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from triqs_cthyb import SolverCore
from pytriqs.gf import GfImFreq, iOmega_n, inverse
import numpy as np

# Input parameters
beta = 10.0
num_orbitals = 2
mu = 1.0
U = 2.0
J = 0.2
epsilon = [-2.3,2.3]
V = [1.0*np.eye(num_orbitals) + 0.1*(np.ones((num_orbitals,num_orbitals)) - np.eye(num_orbitals))]*2

spin_names = ("up","dn")
orb_names = range(num_orbitals)

for use_qn in (True, False):

    n_iw = 1025
    n_tau = 10001

    p = {}
    p["max_time"] = -1
    p["random_name"] = ""
    p["random_seed"] = 123 * mpi.rank + 567
    p["length_cycle"] = 100
    p["n_warmup_cycles"] = 50000
    p["n_cycles"] = int(1.e7 / mpi.size)
    p["use_norm_as_weight"] = True
    p["measure_density_matrix"] = False
    p["move_double"] = True

    results_file_name = "kanamori_offdiag." + ("qn." if use_qn else "") + "h5"

    mpi.report("Welcome to Kanamori (off-diagonal) benchmark.")

    gf_struct = set_operator_structure(spin_names,orb_names,True)
    mkind = get_mkind(True,None)

    # Hamiltonian
    H = h_int_kanamori(spin_names,orb_names,
                       np.array([[0,U-3*J],[U-3*J,0]]),
                       np.array([[U,U-2*J],[U-2*J,U]]),
                       J,True)

    if use_qn:
        QN = [sum([n(*mkind("up",o)) for o in orb_names],Operator()),
              sum([n(*mkind("dn",o)) for o in orb_names],Operator())]
        for o in orb_names:
            dn = n(*mkind("up",o)) - n(*mkind("dn",o))
            QN.append(dn*dn)
        p["partition_method"] = "quantum_numbers"
        p["quantum_numbers"] = QN

    mpi.report("Constructing the solver...")

    # Construct the solver
    S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

    mpi.report("Preparing the hybridization function...")

    # Set hybridization function
    delta_w = GfImFreq(indices = orb_names, beta=beta, n_points=n_iw)
    delta_w_part = delta_w.copy()
    for e, v in zip(epsilon,V):
        delta_w_part << inverse(iOmega_n - e)
        delta_w_part.from_L_G_R(np.transpose(v),delta_w_part,v)
        delta_w += delta_w_part

    S.G0_iw << inverse(iOmega_n + mu - delta_w)

    mpi.report("Running the simulation...")

    # Solve the problem
    S.solve(h_int=H, **p)

    # Save the results
    if mpi.is_master_node():
        with HDFArchive(results_file_name,'w') as Results:
            Results['G_tau'] = S.G_tau
