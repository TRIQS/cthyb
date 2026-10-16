# Check sigma_high_frequency_moments against the exact self-energy moments of
# an impurity with bath, solved by exact diagonalization. The orbital basis is
# rotated so that Sigma_0 has sizable off-diagonal elements.
#
# Exact moments, G(iw) = sum_n G_n / (iw)^n, with H the full Hamiltonian:
#   G_2 = <{[c,H],c+}>,  G_3 = <{[c,H],[H,c+]}>
#   non-interacting: G0_2 = E,  G0_3 = E E + V V+
#   Sigma_0 = G_2 - G0_2,  Sigma_1 = G_3 - G0_3 + G0_2 G0_2 - G_2 G_2

import numpy as np
from itertools import product

from triqs.operators import c, c_dag, Operator
from triqs.atom_diag import AtomDiagComplex, atomic_density_matrix, trace_rho_op
from triqs_cthyb.tail_fit import sigma_high_frequency_moments

rng = np.random.default_rng(7)
beta = 3.0
norb, nbath = 2, 2
U, J = 3.0, 0.6
spins = ('up', 'dn')

# random unitary: c_i (Kanamori basis) = sum_a rot[a, i] d_a
q, r = np.linalg.qr(rng.normal(size=(norb, norb)) + 1j * rng.normal(size=(norb, norb)))
rot = q @ np.diag(np.diag(r) / abs(np.diag(r)))

def ck(s, i):  return sum(rot[a, i] * c(s, a) for a in range(norb))
def ckd(s, i): return sum(np.conj(rot[a, i]) * c_dag(s, a) for a in range(norb))
def nk(s, i):  return ckd(s, i) * ck(s, i)

h_int = Operator()
for i in range(norb):
    h_int += U * nk('up', i) * nk('dn', i)
for i, j in product(range(norb), repeat=2):
    if i == j: continue
    h_int += (U - 2 * J) * nk('up', i) * nk('dn', j)
    if i < j:
        for s in spins: h_int += (U - 3 * J) * nk(s, i) * nk(s, j)
    h_int += -J * ckd('up', i) * ck('dn', i) * ckd('dn', j) * ck('up', j)
    h_int += J * ckd('up', i) * ckd('dn', i) * ck('up', j) * ck('dn', j)

E = rng.normal(size=(norb, norb)) + 1j * rng.normal(size=(norb, norb))
E = 0.5 * (E + E.conj().T)
V = rng.normal(size=(norb, nbath)) + 1j * rng.normal(size=(norb, nbath))
eps = rng.normal(size=nbath)

h_0 = Operator()
for s in spins:
    for a, b in product(range(norb), repeat=2):
        h_0 += E[a, b] * c_dag(s, a) * c(s, b)
    for k in range(nbath):
        h_0 += eps[k] * c_dag('b' + s, k) * c('b' + s, k)
        for a in range(norb):
            h_0 += V[a, k] * c_dag(s, a) * c('b' + s, k) + np.conj(V[a, k]) * c_dag('b' + s, k) * c(s, a)

h = h_0 + h_int
fops = [(s, a) for s in spins for a in range(norb)] + [('b' + s, k) for s in spins for k in range(nbath)]
ad = AtomDiagComplex(h, fops)
rho = atomic_density_matrix(ad, beta)

def comm(A, B): return A * B - B * A
def moment(op_fn, bl):
    return np.array([[trace_rho_op(rho, op_fn(c(bl, a), c_dag(bl, b)), ad)
                      for b in range(norb)] for a in range(norb)])

gf_struct = [(s, norb) for s in spins]
sigma_moments = sigma_high_frequency_moments(rho, ad, gf_struct, h_int)

for bl in spins:
    G_2 = moment(lambda x, y: comm(x, h) * y + y * comm(x, h), bl)
    G_3 = moment(lambda x, y: comm(x, h) * comm(h, y) + comm(h, y) * comm(x, h), bl)
    G0_2, G0_3 = E, E @ E + V @ V.conj().T
    sigma_0 = G_2 - G0_2
    sigma_1 = G_3 - G0_3 + G0_2 @ G0_2 - G_2 @ G_2

    # the test is only meaningful with off-diagonal Sigma_0
    assert np.abs(sigma_0 - np.diag(np.diag(sigma_0))).max() > 1e-2

    np.testing.assert_allclose(sigma_moments[bl][0], sigma_0, atol=1e-10)
    np.testing.assert_allclose(sigma_moments[bl][1], sigma_1, atol=1e-10)
