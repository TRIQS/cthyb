# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

import matplotlib.pyplot as plt

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_field_pyed.h5','r') as A: p = A["p"]

    plt.plot(p.H_vec, p.Sz_vec, '.-')
    plt.plot(p.H_vec, -p.chi * p.H_vec, '-r', lw=0.5)
    plt.xlabel(r'Magnetic field H')
    plt.ylabel(r'$\langle S_z \rangle$')
    plt.show()

