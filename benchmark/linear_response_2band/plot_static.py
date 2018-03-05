  
""" Test calculation for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

 """ 

# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from pytriqs.archive import HDFArchive
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
    with HDFArchive('data_static_pyed.h5','r') as A: static = A["p"]
    with HDFArchive('data_static_pomerol.h5','r') as A: dynamic = A["p"]
    with HDFArchive('data_static_cthyb_t.h5','r') as A: cthyb = A["p"]

    tau = np.array([float(tau) for tau in static.chi_tau.mesh])
    static.chi_tau_int = np.zeros_like(static.chi)
    
    for i1, i2, i3, i4, in itertools.product(range(4), repeat=4):
        print i1, i2, i3, i4
        c_tau = static.chi_tau[i1, i2, i3, i4].data
        static.chi_tau_int[i1, i2, i3, i4] = -np.trapz(c_tau, x=tau)/model.beta

    shape = (16, 16)
    
    import matplotlib.pyplot as plt

    if False:
        plt.figure(figsize=(10, 10))

        cmin = np.min(static.chi_tau.data)
        cmax = np.max(static.chi_tau.data)

        subp = [16, 16, 1]
        for i1, i2, i3, i4, in itertools.product(range(4), repeat=4):
            print i1, i2, i3, i4
            plt.subplot(*subp); subp[-1] += 1
            plt.plot(tau, static.chi_tau[i1, i2, i3, i4].data)
            plt.ylim([cmin, cmax])
            plt.xticks([])
            plt.yticks([])

        #plt.show();exit()

    opt = dict(
        origin='upper',
        cmap=plt.get_cmap('gist_earth_r'),
        vmax = 0.5,
        vmin = 0.0,
        )

    def plot_chi(data, title):
        plt.figure(figsize=(6, 5))
        plt.title(title)
        plt.imshow(data.reshape(shape), **opt)
        plt.colorbar()
        plt.axis('equal')

    plot_chi(
        static.chi,
        r'Static $\chi^{(s)}_{abcd} = \langle c^\dagger_a c_b c^\dagger_c c_d \rangle$'
        )

    plot_chi(
        static.chi_tau_int,
        r'Static integr $\chi^{(s)}_{abcd} = \int_0^\beta \tau \langle (c^\dagger_a c_b)(\tau) c^\dagger_c c_d \rangle$'
        )

    plot_chi(
        dynamic.chi,
        r'From dynamic $\chi^{(d)}_{abcd} = \sum_{nm}\chi_{abcd}(\Omega=0, \nu_n, \nu_m)$'
        )

    plot_chi(
        cthyb.chi,
        r'From dynamic cthyb $\chi^{(d)}_{abcd} = \sum_{nm}\chi_{abcd}(\Omega=0, \nu_n, \nu_m)$'
        )
    
    plt.show()    
        
# ----------------------------------------------------------------------
