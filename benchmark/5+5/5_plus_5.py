#!/bin/env pytriqs

import inspect
import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.U_matrix import cubic_names, U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from triqs_cthyb import SolverCore
import triqs_cthyb.version as version
from pytriqs.gf import GfImFreq, iOmega_n, inverse
from itertools import product

def five_plus_five(use_interaction=True):

    results_file_name = "5_plus_5.h5"

    # Block structure of GF
    L = 2 # d-orbital
    spin_names = ("up","dn")
    orb_names = cubic_names(L)

    # Input parameters
    beta = 40.
    mu = 26

    U = 4.0
    J = 0.7
    F0 = U
    F2 = J*(14.0/(1.0 + 0.63))
    F4 = F2*0.63

    # Dump the local Hamiltonian to a text file (set to None to disable dumping)
    H_dump = "H.txt"
    # Dump Delta parameters to a text file (set to None to disable dumping)
    Delta_dump = "Delta_params.txt"

    # Hybridization function parameters
    # Delta(\tau) is diagonal in the basis of cubic harmonics
    # Each component of Delta(\tau) is represented as a list of single-particle
    # terms parametrized by pairs (V_k,\epsilon_k).
    delta_params={"xy" : {'V':0.2,'e':-0.2},
                  "yz" : {'V':0.2,'e':-0.15},
                  "z^2" : {'V':0.2,'e':-0.1},
                  "xz" : {'V':0.2,'e':0.05},
                  "x^2-y^2" : {'V':0.2,'e':0.4}}

    atomic_levels = {('up_xy',0) :        -0.2,
                     ('dn_xy',0) :        -0.2,
                     ('up_yz',0) :        -0.15,
                     ('dn_yz',0) :        -0.15,
                     ('up_z^2',0) :       -0.1,
                     ('dn_z^2',0) :       -0.1,
                     ('up_xz',0) :        0.05,
                     ('dn_xz',0) :        0.05,
                     ('up_x^2-y^2',0) :   0.4,
                     ('dn_x^2-y^2',0) :   0.4}

    n_iw = 1025
    n_tau = 10001

    p = {}
    p["max_time"] = -1
    p["random_name"] = ""
    p["random_seed"] = 123 * mpi.rank + 567
    p["length_cycle"] = 50
    p["n_warmup_cycles"] = 1000
    p["n_cycles"] = 30000
    p["partition_method"] = "autopartition"
    p["measure_G_tau"] = True
    p["move_shift"] = True
    p["measure_pert_order"] = False
    p["performance_analysis"] = False
    p["use_trace_estimator"] = False


    mpi.report("Welcome to 5+5 (5 orbitals + 5 bath sites) test.")

    gf_struct = set_operator_structure(spin_names,orb_names,False)
    mkind = get_mkind(False,None)

    # Local Hamiltonian
    U_mat = U_matrix(L,[F0,F2,F4],basis='cubic')
    H = h_int_slater(spin_names,orb_names,U_mat,False,H_dump=H_dump)

    p["h_int"] = H

    # Quantum numbers (N_up and N_down)
    QN=[Operator(),Operator()]
    for cn in orb_names:
        for i, sn in enumerate(spin_names):
            QN[i] += n(*mkind(sn,cn))
    if p["partition_method"] == "quantum_numbers": p["quantum_numbers"] = QN

    mpi.report("Constructing the solver...")

    # Construct the solver
    S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

    mpi.report("Preparing the hybridization function...")

    # Set hybridization function
    if Delta_dump: Delta_dump_file = open(Delta_dump,'w')
    for sn, cn in product(spin_names,orb_names):
        bn, i = mkind(sn,cn)
        V = delta_params[cn]['V']
        e = delta_params[cn]['e']

        delta_w = GfImFreq(indices = [i], beta=beta)
        delta_w << (V**2) * inverse(iOmega_n - e)

        S.G0_iw[bn][i,i] << inverse(iOmega_n +mu - atomic_levels[(bn,i)] - delta_w)

        # Dump Delta parameters
        if Delta_dump:
            Delta_dump_file.write(bn + '\t')
            Delta_dump_file.write(str(V) + '\t')
            Delta_dump_file.write(str(e) + '\n')

    mpi.report("Running the simulation...")

    # Solve the problem
    S.solve(**p)

    # Save the results
    if mpi.is_master_node():
        Results = HDFArchive(results_file_name,'w')
        Results['G_tau'] = S.G_tau
        Results['use_interaction'] = use_interaction
        Results['delta_params'] = delta_params
        Results['spin_names'] = spin_names
        Results['orb_names'] = orb_names

        import __main__
        Results.create_group("log")
        log = Results["log"]
        log["version"] = version.version
        log["triqs_hash"] = version.triqs_hash
        log["cthyb_hash"] = version.cthyb_hash
        log["script"] = inspect.getsource(__main__)

if __name__ == '__main__':

    for use_interaction in [True, False]:
        five_plus_five(use_interaction=use_interaction)
