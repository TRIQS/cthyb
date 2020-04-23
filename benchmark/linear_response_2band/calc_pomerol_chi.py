
""" Test calculation for two-band Hubbard atom with two bath sites.

Construct the generalized susceptibility from the Green's functions
and compute the static generalized susceptibility by summing over
fermionic Matsubara frequencies. 

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.operators import Operator
from h5 import HDFArchive
from pyed.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, inverse, iOmega_n, Fourier

# ----------------------------------------------------------------------

from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.chi_from_gg2 import chi_from_gg2_PH

# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_pomerol.h5', 'r') as A: 
        p = A['p']

    p.chi0_nn = chi0_from_gg2_PH(p.g_iw['0'], p.g4_ph)
    p.chi_nn = chi_from_gg2_PH(p.g_iw['0'], p.g4_ph)

    chi = np.sum(p.chi_nn.data, axis=(0, 1, 2)) / p.beta**2
    
    np.testing.assert_array_almost_equal(chi.imag, np.zeros_like(chi.imag))
    p.chi = chi.real
    
    with HDFArchive('data_pomerol_chi.h5', 'w') as A: 
        A['p'] = p

