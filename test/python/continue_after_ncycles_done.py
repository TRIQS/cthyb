# Checks continue_after_ncycles_done under an MPI launcher, where mc_generic installs the cycle monitor.
# The ranks get different n_cycles, so with the default (true) the faster ranks keep measuring until
# the slowest one is done. With false every rank stops after exactly its own n_cycles, and a repeated
# solve with the same seeds is bit-identical.
import numpy as np
import triqs.utility.mpi as mpi
from triqs.gfs import GfImFreq, SemiCircular, inverse, iOmega_n
from triqs.operators import n
from triqs_cthyb import SolverCore, ConstrParametersT, SolveParametersT

cp = dict(beta=10.0, gf_struct=[['up', 1], ['do', 1]], n_iw=200, n_tau=401)
n_cycles = 500 * (mpi.rank + 1)
n_cycles_tot = mpi.all_reduce(n_cycles)

def solve(continue_after_ncycles_done):
    S = SolverCore(ConstrParametersT(**cp))
    delta_w = GfImFreq(beta=cp['beta'], target_shape=[1, 1], n_points=cp['n_iw'])
    delta_w << 0.25 * SemiCircular(1.0)
    for name, g0 in S.G0_iw:
        g0 << inverse(iOmega_n + 0.5 - delta_w)
    S.solve(SolveParametersT(h_int=2.0 * n('up', 0) * n('do', 0),
                             n_cycles=n_cycles,
                             n_warmup_cycles=100,
                             length_cycle=50,
                             random_seed=123 * mpi.rank + 567,
                             measure_pert_order=True,
                             continue_after_ncycles_done=continue_after_ncycles_done))
    h = S.perturbation_order_total
    return S, h.n_data_pts + h.n_lost_pts

S1, n_meas1 = solve(False)
S2, n_meas2 = solve(False)
_, n_meas_cont = solve(True)

mpi.report(f"measures: n_cycles_tot = {n_cycles_tot}, stop = {n_meas1}, {n_meas2}, continue = {n_meas_cont}")

assert S1.last_solve_parameters.continue_after_ncycles_done is False
assert n_meas1 == n_cycles_tot and n_meas2 == n_cycles_tot
assert n_meas_cont >= n_cycles_tot
for (name, g1), (_, g2) in zip(S1.G_tau, S2.G_tau):
    assert np.array_equal(g1.data, g2.data), f"G_tau[{name}] differs between identical runs"
