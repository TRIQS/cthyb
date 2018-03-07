
  """ Test calculation for two-band Hubbard atom with two bath sites.

Use cthyb one- and two-particle Green's function to costruct the 
generalized susceptibility trasformed back to the original basis.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.operators import Operator
from pytriqs.archive import HDFArchive
from pytriqs.gf import Gf, inverse, iOmega_n, InverseFourier, Fourier, GfImFreq

from pyed.ParameterCollection import ParameterCollection
from pyed.GfUtils import g2_single_particle_transform
from pyed.GfUtils import g4_single_particle_transform

# ----------------------------------------------------------------------

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt     

from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.chi_from_gg2 import chi_from_gg2_PH

# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_cthyb.h5', 'r') as A: 
        p = A['p']

    p.g_tau = g2_single_particle_transform(p.g_tau['0'], p.T)

    p.g_iw = GfImFreq(
        name=r'$g$', beta=p.beta,
        statistic='Fermion', n_points=500,
        target_shape=(4, 4))
    
    p.g_iw << Fourier(p.g_tau)

    p.g4_ph = g4_single_particle_transform(p.g4_ph, p.T)

    p.chi0_nn = chi0_from_gg2_PH(p.g_iw, p.g4_ph)
    p.chi_nn = chi_from_gg2_PH(p.g_iw, p.g4_ph)
    
    p.chi = np.sum(p.chi_nn.data, axis=(0, 1, 2)) / p.beta**2

    with HDFArchive('data_cthyb_chi.h5', 'w') as A: 
        A['p'] = p
