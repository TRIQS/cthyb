from cthyb import SolverCore
import pytriqs.utility.mpi as mpi

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001):
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
        SolverCore.__init__(self, beta, gf_struct, n_iw=n_iw, n_tau=n_tau)

        self.Sigma_iw = self.G0_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw = self.G0_iw.copy()
        self.G_iw.zero()

    def solve(self, h_loc, params=None, **params_kw):
        """ Solve the impurity problem """

        # FIXME Default tail fitting parameters
        #known_moments =
        #n_moments =
        #min_n_iw = int(0.8 * n_iw)
        #max_n_iw = n_iw
            
        if params==None:
            if mpi.rank == 0: print "Using keyword arguments provided as parameters in the solver."
            params = SolverCore.solve_parameters()
            for i in params_kw:
                params[i] = params_kw[i]
        else:
            if mpi.rank == 0: print "Using parameter list in the solver."

        # Call the core solver's core routine
        SolverCore.solve(self, h_loc, params)

        # Post-processing:
        # Fourier transform G_tau to obtain G_iw and fit the tail
        for name, g in self.G_tau:
	    self.G_iw[name] <<= Fourier(g)

        # FIXME fit tails bit

        # Solve Dyson's eq to obtain Sigma_iw
        self.Sigma_iw = inverse(self.G0_iw) - inverse(self.G_iw)
