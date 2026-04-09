"""
Test measure_densities against the density_matrix + trace_rho_op approach.

Runs the solver with both measure_density_matrix=True and measure_densities=True,
then checks that the per-orbital densities from the two methods agree within
the statistical error bars of the trace-rho-op measurement.
"""

import triqs.utility.mpi as mpi
from triqs.operators import *
from triqs.atom_diag import trace_rho_op
from triqs_cthyb import Solver
from triqs.gf import *
import numpy as np

# --- Parameters ---
beta = 10.0
U = 2.0
mu = 1.0
h = 0.1       # small magnetic field to break spin symmetry
V = 1.0
epsilon = 2.3

gf_struct = [["dn", 2], ["up", 2]]

p = {}
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 5000
p["n_cycles"] = 50000
p["move_double"] = False
p["measure_G_tau"] = False
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = True
p["measure_densities"] = True

H = U * n("up", 0) * n("dn", 0) + U * n("up", 1) * n("dn", 1)
H += 0.5 * h * (n("up", 0) - n("dn", 0)) + 0.5 * h * (n("up", 1) - n("dn", 1))

# --- Solve ---
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=10001, n_iw=1025)

delta_w = GfImFreq(beta=beta, target_shape=(2, 2))
delta_w << (V**2) * (inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon))
for bn, g in S.G0_iw:
    g << inverse(iOmega_n - np.array([[-mu, 0.1], [0.1, -mu]]) - delta_w)

S.solve(h_int=H, **p)

if mpi.is_master_node():
    # --- density_matrix + trace_rho_op reference ---
    dm = S.density_matrix
    ref = {}
    for bl, bl_size in gf_struct:
        ref[bl] = np.array([trace_rho_op(dm, n(bl, a), S.h_loc_diagonalization).real
                            for a in range(bl_size)])

    # --- measure_densities result ---
    meas = S.densities
    meas_err = S.densities_errors
    assert meas is not None, "S.densities is None — measurement did not run"
    assert meas_err is not None, "S.densities_errors is None"

    print("Comparison: measure_densities vs density_matrix + trace_rho_op")
    all_ok = True
    for bl, bl_size in gf_struct:
        for a in range(bl_size):
            d = meas[bl][a]
            r = ref[bl][a]
            e = meas_err[bl][a]
            diff = abs(d - r)
            n_sigma = diff / e if e > 0 else 0.0
            status = "OK" if n_sigma < 4.0 else "FAIL"
            if status == "FAIL":
                all_ok = False
            print(f"  {bl}[{a}]:  densities={d:.6f}  ref={r:.6f}  err={e:.6f}  diff={diff:.6f}  ({n_sigma:.1f} sigma)  {status}")

    assert all_ok, "measure_densities and density_matrix + trace_rho_op disagree by more than 4 sigma"
    print("\nAll densities agree within 4 sigma.")
