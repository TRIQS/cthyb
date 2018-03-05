  
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

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization
from pyed.GfUtils import g2_single_particle_transform
from pyed.GfUtils import g4_single_particle_transform

# ----------------------------------------------------------------------
def calc_field(H, field_op, exp_op, p):

    H += field_op
    ed = TriqsExactDiagonalization(H, p.op_full, p.beta)
    return ed.get_expectation_value(exp_op)
    
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: m = A["p"]
    else: m = None
    m = mpi.bcast(m)

    p = ParameterCollection()
    
    H_max = 10. / m.beta
    p.H_vec = np.linspace(-H_max, H_max, num=50)
    p.Sz_vec = np.zeros_like(p.H_vec)
    
    for idx, H in enumerate(p.H_vec):
        p.Sz_vec[idx] = calc_field(m.H, H * m.Sz, m.Sz, m)
        print H, p.Sz_vec[idx]

    # ------------------------------------------------------------------
    from scipy.interpolate import InterpolatedUnivariateSpline as IUS

    spl = IUS(p.H_vec, p.Sz_vec)
    p.chi = -spl(0, nu=1) # Linear response
    print 'chi =', p.chi
    
    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():
        with HDFArchive('data_field_pyed.h5','w') as res:
            res['p'] = p
        
# ----------------------------------------------------------------------
