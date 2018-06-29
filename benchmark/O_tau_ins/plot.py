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

    filename = 'data_pyed_h_field_0.0000.h5'
    print '--> Loading:', filename
    with HDFArchive(filename, 'r') as s:
        pyed = s['p']
        
    #oplotr(cthyb.O_tau, label='cthyb shifted O_tau')
    plt.figure(figsize=(3.25, 5))
    
    subp = [2, 1, 1]
    plt.subplot(*subp); subp[-1] += 1
    
    oplotr(pyed.O_tau, label=r'pyed $-\langle n_\uparrow(\tau) n_\downarrow \rangle$')
    oplotr(cthyb.O_tau, '.', label=r'cthyb $-\langle n_\uparrow(\tau) n_\downarrow \rangle$', alpha=0.25)

    plt.subplot(*subp); subp[-1] += 1

    oplotr(pyed.G_tau, label='pyed G_tau')
    oplotr(cthyb.G_tau, '.', label='cthyb G_tau', alpha=0.25)

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('figure_n_n_by_insertion_cthyb.pdf')
    plt.show()
    
