
""" Test calculation for two-band Hubbard atom with two bath sites.

Plot the static suceptibility calculated in three different ways from
exact diagonalization (pyed and pomerol) and from cthyb.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """
  
# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from h5 import HDFArchive
from pytriqs.gf import GfImTime, GfImFreq
from pytriqs.gf import MeshImTime, MeshProduct, Gf
from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization
from pyed.GfUtils import g2_single_particle_transform
from pyed.GfUtils import g4_single_particle_transform

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    with HDFArchive('data_model.h5','r') as A: model = A["p"]
    with HDFArchive('data_pomerol_chi.h5','r') as A: dynamic = A["p"]
    with HDFArchive('data_cthyb_chi.h5','r') as A: cthyb = A["p"]
    with HDFArchive('data_pyed_field.h5','r') as A: pyed = A["p"]

    shape = (16, 16)
    
    import matplotlib.pyplot as plt

    print('min/max chi =', np.min(pyed.chi.real), np.max(pyed.chi.real))
    print('min/max chi_field =', np.min(pyed.chi_field.real), np.max(pyed.chi_field.real))
    print('min/max chi_static =', np.min(pyed.chi_static.real), np.max(pyed.chi_static.real))

    #exit()
        
    opt = dict(
        origin='upper',
        cmap=plt.get_cmap('RdBu'),
        vmax = +0.15,
        vmin = -0.15,
        )

    def plot_chi(data, title, transpose=False):

        plt_data = data.reshape(shape)
        if transpose:
            plt_data = plt_data.T
            
        plt.title(title, fontsize=8)
        plt.imshow(plt_data, **opt)
        plt.axis('equal')
        plt.xticks(range(0, 16, 2))
        plt.yticks(range(0, 16, 2))

    data = [

        (pyed.chi_static.real,
         r'PYED $\langle \bar{a} b \bar{c} d \rangle$'),

        (pyed.chi.real,
         r'PYED $\chi^{(PH)}_{\bar{a}b\bar{c}d} = \frac{1}{\beta}\int_0^\beta \tau \langle [\bar{a} b](\tau) [\bar{c} d] \rangle - \langle \bar{a} b \rangle \langle \bar{c} d \rangle$'),

        (dynamic.chi,
         r'Pomerol $\chi^{(PH)}_{abcd} = \frac{1}{\beta^2} \sum_{nm}\chi^{(PH)}_{\bar{a}b\bar{c}d}(0, \nu_n, \nu_m)$'),

        (cthyb.chi.real,
         r'CTHYB $\chi^{(PH)}_{abcd} = \frac{1}{\beta^2} \sum_{nm}\chi^{(PH)}_{\bar{a}b\bar{c}d}(0, \nu_n, \nu_m)$'),

        (pyed.chi_field.real,
         r'PYED $\partial_{F_{\bar{a}b}} \langle \bar{c} d \rangle$'),

        (pyed.chi.real + pyed.chi.real.swapaxes(0, 1),
         r'PYED $R_{\bar{a}b\bar{c}d} = \chi^{(PH)}_{\bar{a}b\bar{c}d} + \chi^{(PH)}_{\bar{b}a\bar{c}d}$'),
        
        (dynamic.chi + dynamic.chi.swapaxes(0, 1),
         r'Pomerol $R_{\bar{a}b\bar{c}d} = \chi^{(PH)}_{\bar{a}b\bar{c}d} + \chi^{(PH)}_{\bar{b}a\bar{c}d} $'),

        (cthyb.chi.real  + cthyb.chi.real.swapaxes(0, 1),
         r'CTHYB $R_{\bar{a}b\bar{c}d} = \chi^{(PH)}_{\bar{a}b\bar{c}d} + \chi^{(PH)}_{\bar{b}a\bar{c}d} $'),
        ]

    plt.figure(figsize=(6, 10))
    subp = [4, 2, 1]
    for value, title in data:
        plt.subplot(*subp); subp[-1] += 1
        plot_chi(value, title)

        if subp[-1] % 2 == 1:
            plt.colorbar()

    plt.tight_layout()
    plt.savefig('figure_w0_response_2band_kanamori.pdf')
    plt.show()
        
# ----------------------------------------------------------------------
