
  """ Test calculation for two-band Hubbard atom with two bath sites.

Plotter for the different types of two particle propagators in the 
calculation. Nb, slow!, plots 16x16 subplots per propagator. Used for 
debugging purposes...

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import itertools
import numpy as np

import matplotlib.pyplot as plt

# ----------------------------------------------------------------------

from pytriqs.operators import Operator
from pytriqs.archive import HDFArchive
from pyed.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, inverse, iOmega_n, InverseFourier

# ----------------------------------------------------------------------
def plot_g2(g, **kwargs):

    plt.figure(figsize=(10, 10))
    shape = g.target_shape
    N = np.product(shape)
    n = np.ceil(np.sqrt(N))
    subp = [n, n, 1]

    title = kwargs.pop('title')
    component = kwargs.pop('component')

    lab = [r'0\uparrow', r'0\downarrow', r'1\uparrow', r'1\downarrow']
    ca_op = [ '$c_{%s}$' % lab[idx] for idx in range(4) ]
    cc_op = [ '$c^\dagger_{%s}$' % lab[idx] for idx in range(4) ]
    
    for i1, i2, i3, i4 in itertools.product(*[range(s) for s in shape]):         

        print i1, i2, i3, i4
        
        ax = plt.subplot(*subp); subp[-1] += 1
        if (i1, i2, i3, i4) == (0, 0, 0, 0) and title is not None:
            plt.title(title, fontsize=4)
        else:
            op_title = '%s%s%s%s' % (cc_op[i1], ca_op[i2], cc_op[i3], ca_op[i4])
            #plt.title('%i,%i,%i,%i' % (i1, i2, i3, i4), fontsize=5)
            plt.title(op_title, fontsize=5)
            
        ax.set_xticks([])
        ax.set_yticks([])
        
        if component == 'Re':
            data = g.data.real
        elif component == 'Im':
            data = g.data.imag
        else:
            raise NotImplementedError
        d = data[0, :, :, i1, i2, i3, i4]
        #if np.max(np.abs(d)) > 1e-8:
        if True:
            #plt.pcolormesh(d, vmin=vmin, vmax=vmax)
            plt.pcolormesh(d, **kwargs)
                
    plt.tight_layout()
    #plt.savefig('figure_%s.pdf' % title.replace(' ', '_'))
    
# ----------------------------------------------------------------------
def plot_g2_r(g, **kwargs):    
    plot_g2(g, component='Re', **kwargs)

# ----------------------------------------------------------------------
def plot_g2_i(g, **kwargs):
    plot_g2(g, component='Im', **kwargs)
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_pomerol.h5', 'r') as A: 
        pom = A['p']

    with HDFArchive('data_pomerol_chi.h5', 'r') as A: 
        pomc = A['p']
        
    with HDFArchive('data_cthyb.h5', 'r') as A: 
        cthyb = A['p']

    with HDFArchive('data_cthyb_chi.h5', 'r') as A: 
        cthybc = A['p']

    cut = 0.01
    data = pom.g4_ph.data
    vmin_r = np.min(data.real) * cut
    vmax_r = np.max(data.real) * cut
    vmin_i = np.min(data.imag) * cut
    vmax_i = np.max(data.imag) * cut

    vm_r = np.max([np.abs(vmin_r), np.abs(vmax_r)])
    vm_i = np.max([np.abs(vmin_i), np.abs(vmax_i)])

    #opt_r = dict(vmin=vmin_r, vmax=vmax_r)
    #opt_i = dict(vmin=vmin_i, vmax=vmax_i)

    opt_r = dict(vmin=-vm_r, vmax=vm_r, cmap=plt.get_cmap('RdBu'))
    opt_i = dict(vmin=-vm_i, vmax=vm_i, cmap=plt.get_cmap('RdBu'))
    
    print vmin_r, vmax_r
    print vmin_i, vmax_i

    #exit()

    #plot_g2_r(cthyb.g2_iw_ph[('0', '0')], title='g2 cthyb Re', **opt_r)
    #plot_g2_i(cthyb.g2_iw_ph[('0', '0')], title='g2 cthyb Im', **opt_i)

    plot_g2_r(cthybc.g4_ph, title='g2 cthyb Re', **opt_r)
    plot_g2_i(cthybc.g4_ph, title='g2 cthyb Im', **opt_i)
    
    plt.show(); exit()

    plot_g2_r(p.g4_ph, title='g2 pom Re', **opt_r)
    plot_g2_i(p.g4_ph, title='g2 pom Im', **opt_i)

    plt.show(); exit()
    
    plot_g2_r(s.chi0_nn, title='chi0 pom Re', **opt_r)
    plot_g2_i(s.chi0_nn, title='chi0 pom Im', **opt_i)

    plot_g2_r(s.chi_nn, title='chi pom Re', **opt_r)
    plot_g2_i(s.chi_nn, title='chi pom Im', **opt_i)
    
    plt.show()
