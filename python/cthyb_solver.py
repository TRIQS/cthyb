from cthyb import SolverCore
from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi
import numpy as np

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001, n_l=30):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : dict{str:list}
                    Structure of the Green's functions. It must be a
                    dictionary which maps the name of each block of the
                    Green's function as a string to a list of integer
                    indices.
                    For example: { 'up': [1,2,3], 'down', [1,2,3] }.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_tau: integer, optional
               Number of imaginary time points used for the Green's functions.
        n_l: integer, optional
             Number of legendre polynomials to use in accumulations of the Green's functions.
        """
        # Initialise the core solver
        SolverCore.__init__(self, beta, gf_struct, n_iw=n_iw, n_tau=n_tau, n_l=n_l)

        self.Sigma_iw = self.G0_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw = self.G0_iw.copy()
        self.G_iw.zero()
        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_tau = n_tau

    def solve(self, h_loc, **params_kw):
        """
        Solve the impurity problem.
        If measure_g_tau (default = True), G_iw and Sigma_iw will be calculated and their tails fitted.
        In addition to the solver parameters, parameters to control the tail fitting can be provided.

        Parameters
        ----------
        h_loc : Operator object
                The local Hamiltonian of the impurity problem to be solved.
        perform_tail_fit : boolean, optional, default = False
                           Should the tails of Sigma_iw and G_iw be fitted?
        fit_n_moments : integer, optional, default = 3
                        Number of moments to fit in the tail of Sigma_iw.
        fit_known_moments : dict{str:TailGf object}, optional, default = {block_name: TailGf(dim1, dim2, n_moments, order_min)}
                            Known moments of Sigma_iw, given as a TailGf object.
        fit_min_n : integer, optional, default = int(0.8 * self.n_iw)
                    Index of iw from which to start fitting.
        fit_max_n : integer, optional, default = n_iw
                    Index of iw to fit until.
        """

        perform_tail_fit = params_kw.pop("perform_tail_fit", False)
        if perform_tail_fit:
            # If tail parameters provided for Sigma_iw fitting, use them, otherwise use defaults
            if not (("fit_min_n" in params_kw) or ("fit_max_n" in params_kw)):
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
            fit_n_moments = params_kw.pop("fit_n_moments", None)
            fit_known_moments = params_kw.pop("fit_known_moments", None)

        print_warning = False
        for name, indices in self.gf_struct.items():
            dim = len(indices)
            if (self.G0_iw[name].tail[1] != np.eye(dim)).any(): print_warning = True
	if print_warning and mpi.is_master_node():
            warning = ("!--------------------------------------------------------------------------------------!\n"
                       "! WARNING: Some components of your G0_iw do not decay as 1/iw. Continuing nonetheless. !\n"
                       "!--------------------------------------------------------------------------------------!")
            print warning

        # Call the core solver's solve routine
        SolverCore.solve(self, h_loc = h_loc, **params_kw)

        # Post-processing:
        # (only supported for G_tau, to permit compatibility with dft_tools)
        if self.last_solve_parameters["measure_g_tau"] == True:
            # Fourier transform G_tau to obtain G_iw
            for name, g in self.G_tau: self.G_iw[name] << Fourier(g)
            # Solve Dyson's eq to obtain Sigma_iw and G_iw and fit the tail
            self.Sigma_iw = self.dyson(G0_iw=self.G0_iw,G_iw=self.G_iw)
            if perform_tail_fit: self.tail_fit(Sigma_iw=self.Sigma_iw,G0_iw=self.G0_iw,G_iw=self.G_iw,\
                                               fit_min_n=fit_min_n,fit_max_n=fit_max_n,fit_min_w=fit_min_w,fit_max_w=fit_max_w,\
                                               fit_n_moments=fit_n_moments,fit_known_moments=fit_known_moments)

    # Determine one of G0_iw, G_iw and Sigma_iw from other two using Dyson's equation
    def dyson(self,**kwargs):
        """
        Solve Dyson's equation for given two of G0_iw, G_iw and Sigma_iw to yield the third.

        Parameters
        ----------
        G0_iw : Gf, optional
                Non-interacting Green's function.
        G_iw : Gf, optional
               Interacting Green's function.
        Sigma_iw : Gf, optional
                   Self-energy.
        """
        if not (len(kwargs)==2 and set(kwargs.keys())<set(['G0_iw','G_iw', 'Sigma_iw'])):
            raise ValueError, 'dyson: Two (and only two) of G0_iw, G_iw and Sigma_iw must be provided to determine the third.'
        if 'G0_iw' not in kwargs:
            G0_iw = inverse(kwargs['Sigma_iw'] - inverse(kwargs['G_iw']))
            return G0_iw
        elif 'G_iw' not in kwargs:
            G_iw = inverse(inverse(kwargs['G0_iw']) - kwargs['Sigma_iw'])
            return G_iw
        elif 'Sigma_iw' not in kwargs:
            Sigma_iw = inverse(kwargs['G0_iw']) - inverse(kwargs['G_iw'])
            return Sigma_iw

    # Fit tails for Sigma_iw and G_iw.
    # Either give window to fix in terms of matsubara frequencies index (fit_min_n/fit_max_n) or value (fit_min_w/fit_max_w).
    def tail_fit(self,Sigma_iw,G0_iw=None,G_iw=None,fit_min_n=None,fit_max_n=None,fit_min_w=None,fit_max_w=None,fit_n_moments=None,fit_known_moments=None):
        """
        Fit the tails of Sigma_iw and optionally G_iw.

        Parameters
        ----------
        Sigma_iw : Gf
                   Self-energy.
        G0_iw : Gf, optional
                Non-interacting Green's function.
        G_iw : Gf, optional
               Interacting Green's function.
               If G0_iw and G_iw are provided, the tails of G_iw will also be fitted.
        fit_min_n : int, optional, default=int(0.8*len(Sigma_iw.mesh))
                    Matsubara frequency index from which tail fitting should start.
        fit_max_n : int, optional, default=int(len(Sigma_iw.mesh))
                    Matsubara frequency index at which tail fitting should end.
        fit_min_w : float, optional
                    Matsubara frequency from which tail fitting should start.
        fit_max_w : float, optional
                    Matsubara frequency at which tail fitting should end.
        fit_n_moments : int, optional
                        Number of moments to fit in the tail of Sigma_iw.
        fit_known_moments : dict{str:TailGf object}, optional, default = {block_name: TailGf(dim1, dim2, n_moments, order_min)}
                            Known moments of Sigma_iw, given as a TailGf object.
        """

        # Define default tail quantities
        if fit_known_moments is None:
            fit_known_moments = {name:TailGf(sig.N1,sig.N2,0,0) for name, sig in Sigma_iw} # TailGf(dim1, dim2, n_moments, order_min)
        if fit_min_w is not None: fit_min_n = int(0.5*(fit_min_w*Sigma_iw.mesh.beta/np.pi - 1.0))
        if fit_max_w is not None: fit_max_n = int(0.5*(fit_max_w*Sigma_iw.mesh.beta/np.pi - 1.0))
        if fit_min_n is None: fit_min_n = int(0.8*len(Sigma_iw.mesh))
        if fit_max_n is None: fit_max_n = int(len(Sigma_iw.mesh))
        if fit_n_moments is None: fit_n_moments = 3

        # Now fit the tails of Sigma_iw and G_iw
        for name, sig in Sigma_iw: sig.fit_tail(fit_known_moments[name], fit_n_moments, fit_min_n, fit_max_n)
        if (G_iw is not None) and (G0_iw is not None):
            for name, g in G_iw: g.tail = inverse( inverse(G0_iw[name].tail) - Sigma_iw[name].tail )

        return Sigma_iw, G_iw
