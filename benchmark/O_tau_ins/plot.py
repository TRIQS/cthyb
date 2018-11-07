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
if __name__ == '__main__':

    filename = 'data_cthyb.h5'
    print '--> Loading:', filename
    with HDFArchive(filename, 'r') as s:
        cthyb = s['p']

    if False:
        filename = '../../test/python/O_tau_ins.ref.h5'
        print '--> Loading:', filename
        with HDFArchive(filename, 'r') as s:
            O_tau_regr = s['O_tau']
        
    filename = 'data_pyed_h_field_0.0000.h5'
    print '--> Loading:', filename
    with HDFArchive(filename, 'r') as s:
        pyed = s['p']
        
    plt.figure(figsize=(3.25*2, 8))
    
    subp = [3, 1, 1]
    plt.subplot(*subp); subp[-1] += 1
    
    #oplotr(O_tau_regr, '-', label=r'cthyb regr $-\langle n_\uparrow(\tau) n_\downarrow \rangle$ (should be bad)', alpha=1.0, lw=1.0, zorder=100)
    oplotr(cthyb.O_tau, '-', label=r'cthyb $-\langle n_\uparrow(\tau) n_\downarrow \rangle$', alpha=1.0, lw=1.0, zorder=100)
    oplotr(pyed.O_tau, label=r'pyed $-\langle n_\uparrow(\tau) n_\downarrow \rangle$')

    plt.subplot(*subp); subp[-1] += 1
    tau = np.array([ float(t) for t in cthyb.O_tau.mesh])

    tau_ref = np.array([ float(t) for t in pyed.O_tau.mesh])
    O_ref = pyed.O_tau.data.copy()
    O_interp = np.interp(tau, tau_ref, O_ref)

    O_tau_diff = cthyb.O_tau.copy()
    O_tau_diff.data[:] -= O_interp
    O_tau_diff.name = r'$\Delta$' + 'O_tau'

    oplotr(O_tau_diff, '-', label=r'cthyb - pyed', alpha=0.75)
    
    plt.subplot(*subp); subp[-1] += 1

    oplotr(cthyb.G_tau, '.', label='cthyb G_tau', alpha=0.2)
    oplotr(pyed.G_tau, label='pyed G_tau')
    #oplotr(cthyb.Gl_tau, '-', label='cthyb Gl_tau', alpha=0.75)

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('figure_n_n_by_insertion_cthyb.pdf')
    plt.show()
    
