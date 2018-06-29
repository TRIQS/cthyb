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
        
    oplotr((cthyb.O_tau + 1)/2., label='cthyb O_tau')
    oplotr(pyed.O_tau, label='pyed O_tau')

    oplotr(cthyb.G_tau, label='cthyb G_tau')
    oplotr(pyed.G_tau, label='pyed G_tau')
    plt.legend(loc='best')

    plt.show()
    
