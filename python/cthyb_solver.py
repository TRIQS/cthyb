from cthyb import SolverCore
import pytriqs.utility.mpi as mpi

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau_g0=1000, n_tau_g=1000):
        """
        :param beta: Inverse temperature.
        :param gf_struct: Structure of the Green's functions. It must be a
                          dictionary which maps the name of each block of the
                          Green's function as a string to a list of integer
                          indices.
                          For example: { 'up': [1,2,3], 'down', [1,2,3] }.
        :param n_iw: (optional, default = 1025) Number of Matsubara frequencies
                     used for the Green's functions.
        """

        # Initialise the core solver
        self.S = SolverCore(beta, gf_struct, n_tau_g0=n_tau_g0, n_tau_g=n_tau_g)

    def solve(self, h_loc, params=None, **params_kw):
        """ Solve the impurity problem """
            

        if params==None:
            if mpi.rank == 0: print "Using keyword arguments provided as parameters."
            p = self.S.solve_parameters()
            for i in params_kw:
                params[i] = params_kw[i]
        else:
            if mpi.rank == 0: print "Using parameter list."

        # Call the core solver's core routine
        self.S.solve(h_loc=H, params=p)
