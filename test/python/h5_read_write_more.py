
""" CTHYB test for storing h_loc

Author: Hugo U.R. Strand (2019) """

# ----------------------------------------------------------------------    

import pytriqs.utility.mpi as mpi
from pytriqs.gf import inverse, iOmega_n, SemiCircular, BlockGf
from pytriqs.operators import n
from h5 import HDFArchive

from triqs_cthyb import Solver
 
# ----------------------------------------------------------------------
if __name__ == '__main__':

    solv = Solver(beta = 10., gf_struct = [['0',[0]]])
    solv.G0_iw['0'][0,0] << inverse( iOmega_n - 0.5 - SemiCircular(2.) )
    
    solv.solve(
        h_int = 1.0 * n('0', 0),
        n_warmup_cycles = int(1e0),
        n_cycles = int(1e1 / mpi.size),        
        )

    if mpi.is_master_node():

        filename = 'data.h5'

        with HDFArchive(filename, 'w') as res:
            res['solv'] = solv

        with HDFArchive(filename, 'r') as res:
            ref_solv = res['solv']


        cf_attr = [
            'Delta_infty', 'Delta_tau', 'G0_iw', 'G2_iw', 'G2_iw_nfft', 'G2_iw_ph', 'G2_iw_ph_nfft', 'G2_iw_pp', 'G2_iw_pp_nfft', 'G2_iwll_ph', 'G2_iwll_pp', 'G2_tau', 'G_l', 'G_tau', 'G_tau_accum', 'O_tau', 'average_sign', 'constr_parameters', 'density_matrix', 'h_loc', 'last_constr_parameters', 'last_solve_parameters', 'performance_analysis', 'solve_parameters'
            ]

        success = True
        
        for attr in cf_attr:

            val = getattr(solv, attr)
            ref_val = getattr(ref_solv, attr)

            if not (val == ref_val):

                if type(val) == BlockGf:
                    for b, g in val:
                        assert (g.data == ref_val[b].data).all()

                else:
                    print('-'*72)
                    print('Attrib: ', attr)
                    print('Type:   ', type(val))
                
                    print('Value:          ', val)
                    print('Value from hdf5:', ref_val)
                    print('Error: Values differ!')
                    print('-'*72)

                    success = False

        assert( success )
