# ----------------------------------------------------------------------

import os
import glob
import numpy as np

# ----------------------------------------------------------------------
from pytriqs.gf import *
from h5 import HDFArchive

# ----------------------------------------------------------------------    

from pytriqs.operators import n
from pytriqs.operators.util.op_struct import set_operator_structure

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt     

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
if __name__ == '__main__':

    plt.figure(figsize=(3.25, 2*3))
    plt.title('Magnetic Susceptibility\n Hubbard atom + 2 bath sites.')
    
    # ------------------------------------------------------------------
    # -- external field pyed
    
    filenames = glob.glob('pyed_beta*/data_pyed_extrap*.h5')

    style = 'sk'
    for filename in filenames:
        print '--> Loading:', filename

        with HDFArchive(filename, 'r') as s:
            field = s['field']

        plt.plot(1./field.beta, field.chi, style, alpha=0.25)

    plt.plot([], [], style, alpha=0.25, label='ed field')

    # ------------------------------------------------------------------
    # -- dynamic pomerol

    filenames = glob.glob('pomerol_*/*.h5')

    styles = { 2:'.m', 4:'.c', 8:'.y' }

    for filename in filenames:
        print '--> Loading:', filename
        with HDFArchive(filename, 'r') as s:
            p = s['p']
            
        if hasattr(p, 'chi'):
            style = styles[p.n_inu]
            plt.plot(1./p.beta, p.chi, style, alpha=0.5)
        
    for nw, style in styles.items():
        plt.plot([], [], style, alpha=0.5, label='ed dynamic $n_w=%i$' % nw)

    # ------------------------------------------------------------------
    # -- dynamic cthyb

    filenames = glob.glob('cthyb_*/*.h5')

    styles = {
        #(8, 6):'+m',
        (8, 7):'_b',
        (8, 8):'_g',
        (8, 9):'_r',
        }

    for filename in filenames:
        print '--> Loading:', filename
        with HDFArchive(filename, 'r') as s:
            p = s['p']
            
        if hasattr(p, 'chi'):

            nw = p.solve.measure_G2_n_fermionic
            nc = int(np.log10(p.n_cycles))
            
            style = styles[(nw, nc)]

            plt.plot(1./p.beta, p.chi, style, alpha=0.5)

    for (nw, nc), style in styles.items():
        plt.plot([], [], style, alpha=0.5, label='cthyb dynamic $n_w=%i$, $\log n_c=%i$' % (nw, nc))

    # ------------------------------------------------------------------

    plt.legend(loc='best', fontsize=6, frameon=False)
    plt.xlabel(r'$T$')
    plt.ylabel(r'$\chi_m$')

    plt.tight_layout()
    plt.savefig('figure_summary.pdf')
    plt.show()
    
