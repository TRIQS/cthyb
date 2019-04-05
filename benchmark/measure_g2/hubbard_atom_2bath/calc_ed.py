""" Exact diagonalization test calculation 
    for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

 """ 

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf, BlockGf, Block2Gf
from pytriqs.gf import MeshImTime, MeshImFreq, MeshProduct

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

#from pytriqs.applications.susceptibility.Dummy import Dummy
#from pytriqs.applications.susceptibility.fourier import chi4_iw_from_tau
#from pytriqs.applications.susceptibility.fourier import chi3_iw_from_tau
#from pytriqs.applications.susceptibility.fourier import chi2_iw_from_tau
#from pytriqs.applications.susceptibility.fourier import g_iw_from_tau

# ----------------------------------------------------------------------

from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# ----------------------------------------------------------------------
def make_calc(U=10):
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    params = dict(
        beta = 2.0,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.00,
        epsilon2 = 4.00,
        mu = 2.0,
        U = U,
        ntau = 40,
        niw = 15,
        )

    # ------------------------------------------------------------------

    class Dummy():
        def __init__(self):
            pass

    d = Dummy() # storage space
    d.params = params

    print '--> Solving SIAM with parameters'
    for key, value in params.items():
        print '%10s = %-10s' % (key, str(value))
        globals()[key] = value # populate global namespace
    
    # ------------------------------------------------------------------

    up, do = 0, 1
    docc = c_dag(up,0) * c(up,0) * c_dag(do,0) * c(do,0)
    nA = c_dag(up,0) * c(up,0) + c_dag(do,0) * c(do,0)
    nB = c_dag(up,1) * c(up,1) + c_dag(do,1) * c(do,1)
    nC = c_dag(up,2) * c(up,2) + c_dag(do,2) * c(do,2)

    d.H = -mu * nA + epsilon1 * nB + epsilon2 * nC + U * docc + \
        V1 * (c_dag(up,0)*c(up,1) + c_dag(up,1)*c(up,0) + \
              c_dag(do,0)*c(do,1) + c_dag(do,1)*c(do,0) ) + \
        V2 * (c_dag(up,0)*c(up,2) + c_dag(up,2)*c(up,0) + \
              c_dag(do,0)*c(do,2) + c_dag(do,2)*c(do,0) )
    
    # ------------------------------------------------------------------
    # -- Exact diagonalization

    fundamental_operators = [
        c(up,0), c(do,0), c(up,1), c(do,1), c(up,2), c(do,2)]
    
    ed = TriqsExactDiagonalization(d.H, fundamental_operators, beta)
    
    # ------------------------------------------------------------------
    # -- Single-particle Green's functions

    d.G_tau = Gf(mesh=MeshImTime(beta, 'Fermion', ntau), target_shape=[])
    d.G_iw = Gf(mesh=MeshImFreq(beta, 'Fermion', niw), target_shape=[])
    
    ed.set_g2_tau(d.G_tau, c(up,0), c_dag(up,0))
    ed.set_g2_iwn(d.G_iw, c(up,0), c_dag(up,0))
    
    # ------------------------------------------------------------------
    # -- Two particle Green's functions
    
    imtime = MeshImTime(beta, 'Fermion', ntau)
    prodmesh = MeshProduct(imtime, imtime, imtime)
    G2opt = dict(mesh=prodmesh, target_shape=[])

    G02_tau = Gf(name='$G^{(2)}_0(\tau_1, \tau_2, \tau_3)$', **G2opt)
    ed.set_g40_tau(G02_tau, d.G_tau)

    G2_tau = Gf(name='$G^{(2)}(\tau_1, \tau_2, \tau_3)$', **G2opt)
    ed.set_g4_tau(G2_tau, c_dag(up,0), c(up,0), c_dag(do,0), c(do,0))

    G2opt_1111 = dict(mesh=prodmesh, target_shape=[1, 1, 1, 1])
    
    d.G02_tau = Gf(**G2opt_1111)
    d.G2_tau = Gf(**G2opt_1111)

    d.G02_tau.data[:, :, :, 0, 0, 0, 0] = G02_tau.data
    d.G2_tau.data[:, :, :, 0, 0, 0, 0] = G2_tau.data
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_ed.h5'
    with HDFArchive(filename,'w') as res:
        for key, value in d.__dict__.items():
            res[key] = value
        
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc(U=5)
