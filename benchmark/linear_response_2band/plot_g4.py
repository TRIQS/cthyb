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
def plot_g2(g, title=None, component='Re', vmin=None, vmax=None):

    plt.figure(figsize=(10, 10))
    shape = g.target_shape
    N = np.product(shape)
    n = np.ceil(np.sqrt(N))
    subp = [n, n, 1]

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
            plt.pcolormesh(d, vmin=vmin, vmax=vmax)
                
    plt.tight_layout()
    #plt.savefig('figure_%s.pdf' % title.replace(' ', '_'))
    
# ----------------------------------------------------------------------
def plot_g2_r(g, title=None, vmin=None, vmax=None):
    plot_g2(g, title=title, component='Re', vmin=vmin, vmax=vmax)

# ----------------------------------------------------------------------
def plot_g2_i(g, title=None, vmin=None, vmax=None):
    plot_g2(g, title=title, component='Im', vmin=vmin, vmax=vmax)
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive('data_g4_pomerol.h5', 'r') as A: 
        p = A['p']

    cut = 1.0
    data = p.g4_ph.data
    vmin_r = np.min(data.real) * cut
    vmax_r = np.max(data.real) * cut
    vmin_i = np.min(data.imag) * cut
    vmax_i = np.max(data.imag) * cut
    opt_r = dict(vmin=vmin_r, vmax=vmax_r)
    opt_i = dict(vmin=vmin_i, vmax=vmax_i)

    print vmin_r, vmax_r
    print vmin_i, vmax_i

    #exit()

    plot_g2_r(p.g4_ph, title='g2 pom Re', **opt_r)
    plot_g2_i(p.g4_ph, title='g2 pom Im', **opt_i)
    plt.show()
