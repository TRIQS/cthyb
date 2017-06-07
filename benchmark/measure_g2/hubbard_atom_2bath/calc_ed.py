  
""" Exact diagonalization test calculation 
    for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

 """ 

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf, BlockGf, Block2Gf
from pytriqs.gf import MeshImTime, MeshProduct
from pytriqs.gf import GfImTime, GfImFreq

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from pytriqs.applications.susceptibility.Dummy import Dummy
from pytriqs.applications.susceptibility.fourier import chi4_iw_from_tau
from pytriqs.applications.susceptibility.fourier import chi3_iw_from_tau
from pytriqs.applications.susceptibility.fourier import chi2_iw_from_tau
from pytriqs.applications.susceptibility.fourier import g_iw_from_tau

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

    Gopt = dict(beta=beta, statistic='Fermion', indices=[1])    
    d.G_tau = GfImTime(name=r'$G(\tau)$', n_points=ntau, **Gopt)
    d.G_iw = GfImFreq(name='$G(i\omega_n)$', n_points=niw, **Gopt)
    
    ed.set_g2_tau(d.G_tau, c(up,0), c_dag(up,0))
    ed.set_g2_iwn(d.G_iw, c(up,0), c_dag(up,0))

    # chi2pp = + < c^+_u(\tau^+) c_u(0^+) c^+_d(\tau) c_d(0) >
    #        = - < c^+_u(\tau^+) c^+_d(\tau) c_u(0^+) c_d(0) >

    chi2opt = dict(beta=beta, statistic='Fermion', indices=[1], n_points=ntau)    
    d.chi2pp_tau = GfImTime(name=r'$\chi^{(2)}_{PP}(\tau)$', **chi2opt)
    ed.set_g2_tau(d.chi2pp_tau, c_dag(up,0)*c_dag(do,0), c(up,0)*c(do,0))
    d.chi2pp_tau *= -1.0 * -1.0 # commutation sign and gf sign
    d.chi2pp_iw = g_iw_from_tau(d.chi2pp_tau, niw)

    # chi2ph = < c^+_u(\tau^+) c_u(\tau) c^+_d(0^+) c_d(0) >
    
    d.chi2ph_tau = GfImTime(name=r'$\chi^{(2)}_{PH}(\tau)$', **chi2opt)
    #d.chi2ph_tau = Gf(name=r'$\chi^{(2)}_{PH}(\tau)$', **chi2opt)    
    ed.set_g2_tau(d.chi2ph_tau, c_dag(up,0)*c(up,0), c_dag(do,0)*c(do,0))
    d.chi2ph_tau *= -1.0 # gf sign
    d.chi2ph_iw = g_iw_from_tau(d.chi2ph_tau, niw)
    
    # ------------------------------------------------------------------
    # -- Two particle Green's functions
    
    imtime = MeshImTime(beta, 'Fermion', ntau)
    prodmesh = MeshProduct(imtime, imtime, imtime)
    G2opt = dict(mesh=prodmesh, target_shape=[1, 1, 1, 1])

    d.G02_tau = Gf(name='$G^{(2)}_0(\tau_1, \tau_2, \tau_3)$', **G2opt)
    ed.set_g40_tau(d.G02_tau, d.G_tau)
    d.G02_iw = chi4_iw_from_tau(d.G02_tau, niw)

    d.G2_tau = Gf(name='$G^{(2)}(\tau_1, \tau_2, \tau_3)$', **G2opt)
    ed.set_g4_tau(d.G2_tau, c_dag(up,0), c(up,0), c_dag(do,0), c(do,0))
    #ed.set_g4_tau(d.G2_tau, c(up,0), c_dag(up,0), c(do,0), c_dag(do,0)) # <cc^+cc^+>
    d.G2_iw = chi4_iw_from_tau(d.G2_tau, niw)

    # -- trying to fix the bug in the fft for w2

    d.G02_iw.data[:] = d.G02_iw.data[:, ::-1, ...].conj()
    d.G2_iw.data[:] = d.G2_iw.data[:, ::-1, ...].conj()
    
    # ------------------------------------------------------------------
    # -- 3/2-particle Green's functions (equal times)

    prodmesh = MeshProduct(imtime, imtime)
    chi3opt = dict(mesh=prodmesh, target_shape=[1, 1, 1, 1])

    # chi3pp = <c^+_u(\tau) c_u(0^+) c^+_d(\tau') c_d(0) >
    #        = - <c^+_u(\tau) c^+_d(\tau') c_u(0^+) c_d(0) >
    
    d.chi3pp_tau = Gf(name='$\Chi^{(3)}_{PP}(\tau_1, \tau_2, \tau_3)$', **chi3opt)
    ed.set_g3_tau(d.chi3pp_tau, c_dag(up,0), c_dag(do,0), c(up,0)*c(do,0))
    d.chi3pp_tau *= -1.0 # from commutation
    d.chi3pp_iw = chi3_iw_from_tau(d.chi3pp_tau, niw)

    # chi3ph = <c^+_u(\tau) c_u(\tau') c^+_d(0^+) c_d(0) >
    
    d.chi3ph_tau = Gf(name='$\Chi^{(3)}_{PH}(\tau_1, \tau_2, \tau_3)$', **chi3opt)
    ed.set_g3_tau(d.chi3ph_tau, c_dag(up,0), c(up,0), c_dag(do,0)*c(do,0))
    d.chi3ph_iw = chi3_iw_from_tau(d.chi3ph_tau, niw)
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_ed.h5'
    with HDFArchive(filename,'w') as res:
        for key, value in d.__dict__.items():
            res[key] = value
        
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc(U=5)
