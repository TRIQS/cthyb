# ----------------------------------------------------------------------

import os
import glob
import numpy as np

# ----------------------------------------------------------------------
from pytriqs.gf import *
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------    

from pytriqs.operators import n
from pytriqs.operators.util.op_struct import set_operator_structure

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def calc_dynamic(plot=True):

    filenames = glob.glob('data_pomerol*.h5')
    assert( len(filenames) == 1 )
    filename = filenames[0]

    print '--> Loading:', filename
    with HDFArchive(filename, 'r') as s:
        p = s['p']

    p.chi_m = p.G2_iw_ph[('up', 'up')] - p.G2_iw_ph[('up', 'dn')]
    p.chi = np.sum(p.chi_m.data) / p.beta**2

    with HDFArchive(filename, 'w') as s:
        s['p'] = p
    
    print 'beta, chi =', p.beta, p.chi

    if plot: plot_dynamic(p)

# ----------------------------------------------------------------------
def plot_dynamic(p):

    plt.figure(figsize=(3.25*2, 8))

    subp = [3, 2, 1]

    G2_iw_ph = p.G2_iw_ph

    d = np.squeeze(p.G2_iw_ph[('up', 'up')].data)

    lim = np.max([np.abs(d.real), np.abs(d.imag)])
    opt = dict(vmin=-lim, vmax=lim, cmap='PuOr')

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.real, **opt)

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.imag, **opt)

    d = np.squeeze(p.G2_iw_ph[('up', 'dn')].data)

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.real, **opt)

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.imag, **opt)

    d = np.squeeze(p.chi_m.data)
    lim = np.max([np.abs(d.real), np.abs(d.imag)])

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.real, **opt)

    ax = plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(d.imag, **opt)
        
    plt.tight_layout()
    plt.savefig('figure_dynamic.pdf')
        
# ----------------------------------------------------------------------
if __name__ == '__main__':

    paths = glob.glob('pomerol_*')
    #paths = glob.glob('pomerol_nw8_beta*')
    #paths = glob.glob('pomerol_nw8_beta0.1000*')

    for path in paths:
        print '--> path:', path
        os.chdir(path)
        calc_dynamic(plot=False)
        os.chdir('../')

    plt.show()
    
