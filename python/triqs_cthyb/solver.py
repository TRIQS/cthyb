################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2017 by H. UR Strand, P. Seth, I. Krivenko,
#                       M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from solver_core  import SolverCore
from pytriqs.gf import *
import pytriqs.utility.mpi as mpi
import numpy as np

from tail_fit import tail_fit as cthyb_tail_fit

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001, n_l=30):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ (str,[int,...]), ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and a list of integer
                    indices.
                    For example: ``[ ('up', [0, 1, 2]), ('down', [0, 1, 2]) ]``.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_tau : integer, optional
               Number of imaginary time points used for the Green's functions.
        n_l : integer, optional
             Number of legendre polynomials to use in accumulations of the Green's functions.
        """
        if isinstance(gf_struct,dict):
            if mpi.is_master_node(): print "WARNING: gf_struct should be a list of pairs [ (str,[int,...]), ...], not a dict"
            gf_struct = [ [k, v] for k, v in gf_struct.iteritems() ]

        # Initialise the core solver
        SolverCore.__init__(self, beta=beta, gf_struct=gf_struct, 
                            n_iw=n_iw, n_tau=n_tau, n_l=n_l)

        self.Sigma_iw = self.G0_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw = self.G0_iw.copy()
        self.G_iw.zero()
        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_tau = n_tau

    def solve(self, **params_kw):
        """
        Solve the impurity problem.
        If ``measure_G_tau`` (default = ``True``), ``G_iw`` and ``Sigma_iw`` will be calculated and their tails fitted.
        In addition to the solver parameters, parameters to control the tail fitting can be provided.

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core solver.
                     Two required :ref:`parameters <solve_parameters>` are
                        * `h_int` (:ref:`Operator object <triqslibs:operators>`): the local Hamiltonian of the impurity problem to be solved,
                        * `n_cycles` (int): number of measurements to be made.
        perform_post_proc : boolean, optional, default = ``True``
                            Should ``G_iw`` and ``Sigma_iw`` be calculated?
        perform_tail_fit : boolean, optional, default = ``False``
                           Should the tails of ``Sigma_iw`` and ``G_iw`` be fitted?
        fit_max_moment : integer, optional, default = 3
                         Highest moment to fit in the tail of ``Sigma_iw``.
        fit_known_moments : ``ndarray.shape[order, Sigma_iw[0].target_shape]``, optional, default = None
                            Known moments of Sigma_iw, given as an numpy ndarray
        fit_min_n : integer, optional, default = ``int(0.8 * self.n_iw)``
                    Index of ``iw`` from which to start fitting.
        fit_max_n : integer, optional, default = ``n_iw``
                    Index of ``iw`` to fit until.
        """

        # -- Deprecation checks for measure parameters
        
        depr_params = dict(
            measure_g_tau='measure_G_tau',
            measure_g_l='measure_G_l',
            )
        
        for key in depr_params.keys():
            if key in params_kw.keys():
                print 'WARNING: cthyb.solve parameter %s is deprecated use %s.' % \
                    (key, depr_params[key])
                val = params_kw.pop(key)
                params_kw[depr_params[key]] = val

        # -- Tail post proc flags
                
        perform_post_proc = params_kw.pop("perform_post_proc", True)
        perform_tail_fit = params_kw.pop("perform_tail_fit", False)
        if perform_post_proc and perform_tail_fit:
            # If tail parameters provided for Sigma_iw fitting, use them, otherwise use defaults
            if not (("fit_min_n" in params_kw) or ("fit_max_n" in params_kw) or ("fit_max_w" in params_kw) or ("fit_min_w" in params_kw)):
	        if mpi.is_master_node():
                    warning = ("!------------------------------------------------------------------------------------!\n"
                               "! WARNING: Using default high-frequency tail fitting parameters in the CTHYB solver. !\n"
                               "! You should check that the fitting range is suitable for your calculation!          !\n"
                               "!------------------------------------------------------------------------------------!")
                    print warning
            fit_min_n = params_kw.pop("fit_min_n", None)
            fit_max_n = params_kw.pop("fit_max_n", None)
            fit_min_w = params_kw.pop("fit_min_w", None)
            fit_max_w = params_kw.pop("fit_max_w", None)
            fit_max_moment = params_kw.pop("fit_max_moment", None)
            fit_known_moments = params_kw.pop("fit_known_moments", None)

        # Check fundamental Green function property G(iw)[i,j] = G(-iw)*[j,i]
        if not is_gf_hermitian(self.G0_iw):
            if mpi.is_master_node():
                warning = ( "!-------------------------------------------------------------------------------------------!\n"
                            "! WARNING: S.G0_iw violates fundamental Green Function property G0(iw)[i,j] = G0(-iw)*[j,i] !\n"
                            "! Symmetrizing S.G0_iw ...                                                                  !\n"
                            "!-------------------------------------------------------------------------------------------!")
                print warning
            self.G0_iw = make_hermitian(self.G0_iw)

        # Call the core solver's solve routine
        solve_status = SolverCore.solve(self, **params_kw)

        # Post-processing:
        # (only supported for G_tau, to permit compatibility with dft_tools)
        if perform_post_proc and (self.last_solve_parameters["measure_G_tau"] == True):
            # Fourier transform G_tau to obtain G_iw
            for name, g in self.G_tau:
                bl_size = g.target_shape[0]
                known_moments = make_zero_tail(g, 4)
                known_moments[1,...] = np.eye(bl_size)
                self.G_iw[name].set_from_fourier(g, known_moments)

            # We enforce the fundamental Green function property G(iw)[i,j] = G(-iw)*[j,i]
            # for the output Green function and store the symmetry violation to self.hermiticity_violation_G_iw
            G_iw_hermitized = make_hermitian(self.G_iw)
            self.hermiticity_violation_G_iw = self.G_iw - G_iw_hermitized
            self.G_iw = G_iw_hermitized
            self.G_iw_raw = self.G_iw.copy()

            # Solve Dyson's eq to obtain Sigma_iw and G_iw and fit the tail
            self.Sigma_iw = dyson(G0_iw=self.G0_iw, G_iw=self.G_iw)
            self.Sigma_iw_raw = self.Sigma_iw.copy()

            if perform_tail_fit:
                
                cthyb_tail_fit(
                    Sigma_iw=self.Sigma_iw,
                    fit_min_n = fit_min_n, fit_max_n = fit_max_n,
                    fit_min_w = fit_min_w, fit_max_w = fit_max_w,
                    fit_max_moment = fit_max_moment,
                    fit_known_moments = fit_known_moments,
                    )

                # Recompute G_iw with the fitted Sigma_iw
                self.G_iw = dyson(G0_iw=self.G0_iw, Sigma_iw=self.Sigma_iw)
            else:

                # Enforce 1/w behavior of G_iw in the tail fit window
                # and recompute Sigma_iw
                for name, g in self.G_iw:
                    tail = np.zeros([2] + list(g.target_shape), dtype=np.complex)
                    tail[1] = np.eye(g.target_shape[0])
                    g.replace_by_tail_in_fit_window(tail)

                self.Sigma_iw = dyson(G0_iw=self.G0_iw, G_iw=self.G_iw)

        return solve_status
