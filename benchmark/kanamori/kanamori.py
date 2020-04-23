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
V = 1.0
epsilon = 2.3

spin_names = ("up","dn")
orb_names = range(num_orbitals)

for use_qn in (True, False):

    n_iw = 1025
    n_tau = 10001

    p = {}
    p["max_time"] = -1
    p["random_name"] = ""
    p["random_seed"] = 123 * mpi.rank + 567
    p["length_cycle"] = 50
    p["n_warmup_cycles"] = 50000
    p["n_cycles"] = 3000000

    results_file_name = "kanamori" + (".qn" if use_qn else "") + ".h5"

    mpi.report("Welcome to Kanamori benchmark.")

    gf_struct = set_operator_structure(spin_names,orb_names,False)
    mkind = get_mkind(False,None)

    ## Hamiltonian
    H = h_int_kanamori(spin_names,orb_names,
                       np.array([[0,U-3*J],[U-3*J,0]]),
                       np.array([[U,U-2*J],[U-2*J,U]]),
                       J,False)

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
    delta_w = GfImFreq(indices = [0], beta=beta, n_points=n_iw)
    delta_w << (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
    S.G0_iw << inverse(iOmega_n + mu - delta_w)

    mpi.report("Running the simulation...")

    # Solve the problem
    S.solve(h_int=H, **p)

    # Save the results  
    if mpi.is_master_node():
        with HDFArchive(results_file_name,'w') as Results:
            Results['G_tau'] = S.G_tau
