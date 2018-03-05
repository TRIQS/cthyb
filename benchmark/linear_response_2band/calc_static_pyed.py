  
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
def calc_field(H, field_op, exp_op, p):

    H += field_op
    
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: m = A["p"]
    else: m = None
    m = mpi.bcast(m)

    p = ParameterCollection()
    ed = TriqsExactDiagonalization(m.H, m.op_full, m.beta)

    p.chi = np.zeros((4, 4, 4, 4))
    p.chi_tau = GfImTime(name=r'$g$', beta=m.beta,
                         statistic='Boson', n_points=50,
                         target_shape=(4, 4, 4, 4))

    print m.op_imp.reverse()

    for i1, i2, i3, i4, in itertools.product(range(4), repeat=4):

        print i1, i2, i3, i4
        
        o1, o2, o3, o4 = m.op_imp[i1], m.op_imp[i2], m.op_imp[i3], m.op_imp[i4]

        O1 = dagger(o1) * o2
        O2 = dagger(o3) * o4

        p.chi[i1, i2, i3, i4] = ed.get_expectation_value(O1 * O2)
        ed.set_g2_tau(p.chi_tau[i1, i2, i3, i4], O1, O2)
        
    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():
        with HDFArchive('data_static_pyed.h5','w') as res:
            res['p'] = p
        
# ----------------------------------------------------------------------
