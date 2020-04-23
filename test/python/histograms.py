#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from h5 import HDFArchive
from pytriqs.operators import *
from triqs_cthyb import *
from pytriqs.statistics.histograms import Histogram
from pytriqs.gf import *
import numpy as np

spin_names = ("up","dn")
mkind = lambda sn: (sn,0)
gf_struct = [["dn",[0]], ["up",[0]]]

# Input parameters
beta = 10.0
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
p["n_warmup_cycles"] = 5000
p["n_cycles"] = 50000
p["move_shift"] = True
p["move_double"] = True
p["measure_pert_order"] = True
p["performance_analysis"] = True

H = U*n("up",0)*n("dn",0) -mu*(n("up",0) + n("dn",0))

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w << (V**2)*(inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon))
for sn in spin_names: S.G0_iw[sn] << inverse(iOmega_n - delta_w)

# Solve the problem
S.solve(h_int=H, **p)

from pytriqs.utility.comparison_tests import *

def assert_histograms_are_close(hi1, hi2):
    assert hi1.n_data_pts == hi2.n_data_pts
    assert hi1.n_lost_pts == hi2.n_lost_pts
    assert hi1.limits == hi2.limits
    assert_arrays_are_close(hi1.data, hi2.data)

if mpi.is_master_node():
    with HDFArchive('histograms.out.h5','w') as ar:
        ar['perturbation_order'] = S.perturbation_order
        ar['perturbation_order_total'] = S.perturbation_order_total
        ar['performance_analysis'] = S.performance_analysis

    with HDFArchive('histograms.ref.h5','r') as ar:
        assert_histograms_are_close(ar['perturbation_order_total'], S.perturbation_order_total)
        for block, h in ar['perturbation_order'].items():
            assert_histograms_are_close(h, S.perturbation_order[block])
        for name, h in ar['performance_analysis'].items():
            assert_histograms_are_close(h, S.performance_analysis[name])
