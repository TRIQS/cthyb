""" 
Test the tail fit on Matsubara frequency Green's fucntions with 
noise growing quadratically with the frequency.

Author: Hugo U.R. Strand """

import itertools
import numpy as np

from pytriqs.archive import HDFArchive
from pytriqs.gf import Gf, MeshImFreq, MeshImTime, iOmega_n, inverse, Fourier

nw = 512
beta = 50.0
target_shape = [2, 2]

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
Delta_iw = Gf(mesh=wmesh, target_shape=target_shape)

Ek = np.array([
    [1.00,  0.5],
    [0.5, -1.20],
    ])

E_loc = np.array([
    [0.33,  0.5],
    [0.5, -0.1337],
    ])

V = np.array([
    [1.0, 0.25],
    [0.25, -1.0],
    ])

Delta_iw << inverse( iOmega_n - Ek ) + inverse( iOmega_n + Ek ) - E_loc
Delta_iw.from_L_G_R(V, Delta_iw, V)

Delta_iw_ref = Delta_iw.copy()

noise_ampl = 1e-3

def get_noise():
    return 2.*np.random.random() - 1.

for w in wmesh:
    re = Delta_iw[w].real + noise_ampl * get_noise() * np.abs(w.value)**2
    im = Delta_iw[w].imag + noise_ampl * get_noise() * np.abs(w.value)**2
    
    Delta_iw[w] = re + 1.j * im

order_max = 3
n_min = 20
n_max = 100

w_vec = np.array([ w for w in wmesh])[nw:]
print len(w_vec)

w_min = w_vec[n_min].value
w_max = w_vec[n_max].value

print w_min, w_max
#exit()

# -- triqs/unstable implementation

if False:
    from pytriqs.gf import BlockGf
    from pytriqs.gf.tools import tail_fit
    Delta_iw_fit = Delta_iw.copy()
    Delta_iw_fit_bgf = BlockGf(name_list=[0], block_list=[Delta_iw_fit])
    tail_fit(Delta_iw_fit_bgf, fit_min_n=n_min, fit_max_n=n_max, fit_max_moment=order_max)
    filename = 'data_tail_fit_example_unstable.h5'
    figure_filename = 'figure_tail_fit_example_unstable.pdf'

# -- triqs/new_tail implementation
elif False:
    from pytriqs.gf.gf_fnt import fit_tail_on_window, replace_by_tail

    known_moments = np.zeros((0, 2, 2), dtype=np.complex) # no known moments
    print 'known_moments.shape =', known_moments.shape
    
    tail, err = fit_tail_on_window(
        Delta_iw,
        n_min = n_min,
        n_max = n_max,
        known_moments = known_moments,
        n_tail_max = 10 * len(Delta_iw.mesh),
        expansion_order=order_max
        )

    print 'tail =', tail

    Delta_iw_fit = Delta_iw.copy()
    replace_by_tail(Delta_iw_fit, tail, n_min=n_min)
    filename = 'data_tail_fit_example_new_tail.h5'
    figure_filename = 'figure_tail_fit_example_new_tail.pdf'
    print tail.shape
    print err

else:
    from pytriqs.gf import BlockGf
    Delta_iw_fit = Delta_iw.copy()
    Delta_iw_fit_bgf = BlockGf(name_list=['foo'], block_list=[Delta_iw_fit])

    #from triqs_cthyb.tail_fit import tail_fit as cthyb_tail_fit
    from cthyb.tail_fit import tail_fit as cthyb_tail_fit
    cthyb_tail_fit(Delta_iw_fit_bgf, fit_min_n=n_min, fit_max_n=n_max, fit_max_moment=order_max)
    
    filename = 'data_tail_fit_example_cthyb.h5'
    figure_filename = 'figure_tail_fit_example_cthyb.pdf'

# -- Store to hdf5

with HDFArchive(filename, 'w') as a:
    a['order_max'] = order_max
    a['n_min'], a['n_max'] = n_min, n_max
    a['Delta_iw'] = Delta_iw
    a['Delta_iw_ref'] = Delta_iw_ref
    a['Delta_iw_fit'] = Delta_iw_fit

# -- Plot 

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt
plt.figure(figsize=(10, 10))

ylim = [-4, 4]

plt.plot([w_min.imag]*2, ylim, 'sr-', lw=0.5)
plt.plot([w_max.imag]*2, ylim, 'sr-', lw=0.5)
         
for i1, i2 in itertools.product(range(2), repeat=2):
    oplotr(Delta_iw[i1, i2], alpha=0.1)
    oploti(Delta_iw[i1, i2], alpha=0.1)

    oplotr(Delta_iw_ref[i1, i2], lw=4, alpha=0.25)
    oploti(Delta_iw_ref[i1, i2], lw=4, alpha=0.25)

    oplotr(Delta_iw_fit[i1, i2])
    oploti(Delta_iw_fit[i1, i2])

ax = plt.gca()
ax.legend_ = None

plt.ylim([-1, 1])
plt.tight_layout()
plt.savefig(figure_filename)
plt.show(); exit()

