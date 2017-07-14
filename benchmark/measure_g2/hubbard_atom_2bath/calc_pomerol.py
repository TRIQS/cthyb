# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation 
    for a Hubbard atom with two bath sites. 

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED

# ----------------------------------------------------------------------
def make_calc():
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    params = dict(
        beta = 2.0,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.00,
        epsilon2 = 4.00,
        mu = 2.0,
        U = 5.0,
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

    up, do = 'up', 'dn'
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

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        (up, 0) : ('loc', 0, 'up'),
        (do, 0) : ('loc', 0, 'down'),
        (up, 1) : ('loc', 1, 'up'),
        (do, 1) : ('loc', 1, 'down'),
        (up, 2) : ('loc', 2, 'up'),
        (do, 2) : ('loc', 2, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(d.H) # -- Diagonalize H

    gf_struct = {up : [0], do : [0]}

    # -- Single-particle Green's functions
    G_iw = ed.G_iw(gf_struct, beta, n_iw=niw)
    G_tau = ed.G_tau(gf_struct, beta, n_tau=ntau)
    G_w = ed.G_w(gf_struct, beta, energy_window=(-2.5, 2.5), n_w=100, im_shift=0.01)

    d.G_iw = G_iw['up']
    d.G_tau = G_tau['up']
    d.G_w = G_w['up']

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        block_order='AABB',
        beta=beta, gf_struct=gf_struct,
        blocks=set([("up", "dn")]),
        n_iw=niw, n_inu=niw)
    
    G2_iw = ed.G2_iw_inu_inup(channel='AllFermionic', **opt)
    d.G2_iw = G2_iw['up', 'dn']

    G2_iw_pp = ed.G2_iw_inu_inup(channel='PP', **opt)
    d.G2_iw_pp = G2_iw_pp['up', 'dn']

    G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)
    d.G2_iw_ph = G2_iw_ph['up', 'dn']
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pomerol.h5'
    with HDFArchive(filename,'w') as res:
        for key, value in d.__dict__.items():
            res[key] = value
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
