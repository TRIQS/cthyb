# %load run_single_band.py
from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import *
from triqs_cthyb import Solver
import numpy as np

# Parameters of the model
t = 1.0
U = 5.0
beta = 10.0
n_loops = 10

# Construct the impurity solver
S = Solver(beta = beta, gf_struct = [('up',[0]), ('down',[0])] )

# This is a first guess for G
S.G_iw << SemiCircular(2*t)

tail_fit = dict(
    perform_tail_fit = True,
    fit_max_moment = 3,
    fit_min_w = 1.2,
    fit_max_w = 3.0,
    )

# DMFT loop with self-consistency
for do_tail in [False, True]:
    for i in range(n_loops):

        print "\n\nIteration = %i / %i" % (i+1, n_loops)

        # Symmetrize the Green's function and use self-consistency
        g = 0.5 * ( S.G_iw['up'] + S.G_iw['down'] )
        for name, g0 in S.G0_iw:
            g0 << inverse( iOmega_n + U/2.0 - t**2 * g )

        # Solve the impurity problem
        opts = tail_fit if do_tail else {}
        S.solve(h_int = U * n('up',0) * n('down',0),   # Local Hamiltonian 
            n_cycles  = 1000000,                       # Number of QMC cycles
            n_warmup_cycles = 5000,                    # Warmup cycles
            **opts
            )

        # Save iteration in archive
        with HDFArchive("data-U%.2f_tail_fit_%s.h5" % (U, do_tail)) as A:

            A['niter'] = i

            A['U'] = U
            A['beta'] = beta

            A['G0_iw-%i'%i] = S.G0_iw
            A['Delta_tau-%i'%i] = S.Delta_tau

            A['G_iw-%i'%i] = S.G_iw
            A['G_iw_raw-%i'%i] = S.G_iw_raw

            A['Sigma_iw-%i'%i] = S.Sigma_iw
            A['Sigma_iw_raw-%i'%i] = S.Sigma_iw_raw
