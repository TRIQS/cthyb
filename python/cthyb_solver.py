from cthyb import SolverCore
from pytriqs.gf.local import *
from numpy import identity
import pytriqs.utility.mpi as mpi

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001, n_l=30):
        """
        :param beta: Inverse temperature.
        :param gf_struct: Structure of the Green's functions. It must be a
                          dictionary which maps the name of each block of the
                          Green's function as a string to a list of integer
                          indices.
                          For example: { 'up': [1,2,3], 'down', [1,2,3] }.
        :param n_iw: (optional, default = 1025) Number of Matsubara frequencies
                     used for the Green's functions.
        :param n_tau: (optional, default = 10001) Number of imaginary time points
                     used for the Green's functions.
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
        """ Solve the impurity problem """

        # Default tail fitting parameters
        fit_known_moments = {}
        for name, indices in self.gf_struct.items():
            dim = len(indices)
            fit_known_moments[name] = TailGf(dim,dim,1,1) # TailGf(dim1, dim2, n_moments, order_min)
            fit_known_moments[name][1] = identity(dim) # 1/w behaviour
        fit_n_moments = 3
        fit_min_n = int(0.8 * self.n_iw) # Fit last 80% of frequencies
        fit_max_n = self.n_iw

        # Call the core solver's core routine
        SolverCore.solve(self, h_loc = h_loc, **params_kw)

        # Post-processing: 
        # (only supported for G_tau, to permit compatibility with dft_tools)
        if params_kw.get("measure_g_tau",True):
         # Fourier transform G_tau to obtain G_iw and fit the tail
         for name, g in self.G_tau:
             self.G_iw[name] <<= Fourier(g)
             self.G_iw[name].fit_tail(fit_known_moments[name], fit_n_moments, fit_min_n, fit_max_n)
 
         # Solve Dyson's eq to obtain Sigma_iw
         self.Sigma_iw = inverse(self.G0_iw) - inverse(self.G_iw)
