# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.operators import Operator
from pytriqs.archive import HDFArchive
from pytriqs.gf import Gf, inverse, iOmega_n, InverseFourier

from pyed.ParameterCollection import ParameterCollection
from pyed.GfUtils import g4_single_particle_transform

# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_cthyb_g2_Ht.h5', 'r') as A: 
        p = A['p']

    p.g4_ph = g4_single_particle_transform(p.g2_iw_ph[('0', '0')], p.T)
    
    chi = np.sum(p.g4_ph.data, axis=(0, 1, 2)) / p.beta**2

    p.chi_re_abs_max = np.max(np.abs(p.g4_ph.data.real), axis=(0, 1, 2))
    p.chi_im_abs_max = np.max(np.abs(p.g4_ph.data.imag), axis=(0, 1, 2))

    print p.chi_re_abs_max.shape
    print p.chi_im_abs_max.shape

    print np.max(p.chi_re_abs_max)
    print np.max(p.chi_im_abs_max)

    #print chi.imag
    #print chi.real
    
    np.testing.assert_array_almost_equal(chi.imag, np.zeros_like(chi.imag))
    p.chi = chi.real
    
    with HDFArchive('data_static_cthyb_t.h5', 'w') as A: 
        A['p'] = p

