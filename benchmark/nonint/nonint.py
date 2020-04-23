#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs_cthyb import SolverCore
from pytriqs.operators import Operator, n
from pytriqs.gf import Gf, inverse, iOmega_n

mpi.report("Welcome to nonint (non-interacting many-band systems) test.")
mpi.report("This test is aimed to reveal excessive state truncation issues.")

beta = 40.0
N_max = 10

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 1200000

for modes in range(1,N_max+1):
    V = [0.2]*modes
    e = [-0.2]*modes

    #gf_struct = {str(n):[0] for n in range(0,len(V))}
    gf_struct = [ [str(bidx), [0]] for bidx in range(0,len(V)) ]

    # Local Hamiltonian
    H = Operator()

    # Quantum numbers (N_up and N_down)
    QN = []
    for b, idxs in gf_struct: QN.append(n(b,0))
    p["partition_method"] = "quantum_numbers"
    p["quantum_numbers"] = QN
    
    mpi.report("Constructing the solver...")

    # Construct the solver
    S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

    mpi.report("Preparing the hybridization function...")

    # Set hybridization function
    #for m, b in enumerate(sorted(gf_struct.keys())):
    for m, (b, idxs) in enumerate(gf_struct):
        delta_w = Gf(mesh=S.G0_iw.mesh, target_shape=[])
        delta_w << (V[m]**2) * inverse(iOmega_n - e[m])
        S.G0_iw[b][0,0] << inverse(iOmega_n - e[m] - delta_w)

    mpi.report("Running the simulation...")

    # Solve the problem
    S.solve(h_int=H, **p)

    # Save results
    if mpi.is_master_node():
        with HDFArchive('nonint.h5','a') as Results:
            Results[str(modes)] = {'G_tau':S.G_tau,'V':V,'e':e}
