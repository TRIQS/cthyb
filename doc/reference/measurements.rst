.. _measurements:

Measurements: definitions
=========================

Here we list all the observables that can be measured by the solver along with their definitions.
Each measurement can be turned on or off via the corresponding :doc:`solve() parameters <parameters>`.

Green's function
----------------

The single-particle Green's function can be measured in imaginary time or the Legendre basis
(see `L. Boehnke et al., Phys. Rev. B 84, 075145 (2011) <http://link.aps.org/doi/10.1103/PhysRevB.84.075145>`_).

Imaginary time
**************

The imaginary-time Green's function is defined as

.. math::

    G^A_{ij}(\tau) = -\langle \mathcal{T}_\tau c_{Ai}(\tau)c_{Aj}^\dagger(0) \rangle,

where :math:`A` is a block index and :math:`i,j` are inner indices within the block.
The block structure of :math:`G^A_{ij}(\tau)`, i.e. valid values of the indices, is fixed by
constructor's parameter ``gf_struct``.

This measurement is turned on by setting ``measure_g_tau`` to ``True``.
The result of accumulation is accessible as ``G_tau`` attribute of the solver object.
The number of time points on the grid is specified through constructor's parameter ``n_tau`` .

.. note::

    The imaginary-time measurement is most efficient. The performance of the algorithm does not scale
    with the number of points in the grid on which it is measured, so this number can and should be
    chosen large. By Nyquist's theorem, the Fourier transform will correctly reproduce the function
    in the frequency domain on the first :math:`N_\omega\approx N_\tau/4\pi` frequencies.

Legendre polynomial basis
*************************

The Legendre coefficients of the Green's function are defined as

.. math::

    G^A_{ij}(l) = \sqrt{2l+1}\int_0^\beta d\tau\, P_l[x(\tau)] G^A_{ij}(\tau),

where :math:`x(\tau)=2\tau/\beta-1`, and :math:`P_l(x)` are Legendre polynomials of order :math:`l`.

This measurement is controlled through the switch ``measure_g_l``.
The result of accumulation is accessible as ``G_l`` attribute of the solver object.
The number of Legendre coefficients to be measured is specified through constructor's parameter ``n_l``.

Two-particle Green's functions
------------------------------

CTHYB implements several ways to measure the two-particle Green's functions.
We define :math:`G^{(2)}` in the imaginary time as

.. math::

    G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4) =
    \langle\mathcal{T}_\tau c^\dagger_\alpha(\tau_1) c_\beta(\tau_2) c^\dagger_\gamma(\tau_3) c_\delta(\tau_4)\rangle.

.. NOTE::

    Are we actually going to measure this one?

Depending on the value of the ``measure_g2_block_order`` parameter the combined block/inner indices
:math:`\alpha,\beta,\gamma,\delta` are interpreted differently.

* ``measure_g2_block_order = 'AABB'``: :math:`\{\alpha,\beta,\gamma,\delta\} = \{(A,a),(A,b),(B,c),(B,d)\}`;
* ``measure_g2_block_order = 'ABBA'``: :math:`\{\alpha,\beta,\gamma,\delta\} = \{(A,a),(B,b),(B,c),(A,d)\}`.

These two combinations exhaust all the possibilities for non-vanishing contributions to :math:`G^{(2)}`.

Measuring :math:`G^{(2)}` is in general expensive, and in some cases time can be saved by skipping
particular block pairs :math:`(A,B)`. Block pairs to be measured are selected via the ``measure_g2_blocks``
parameter (all possible pairs are selected by default):

``measure_g2_blocks = set([("A1","B1"), ("A2","B2"), ...])``

One bosonic and two fermionic Matsubara frequencies
***************************************************

This measurement is activated by setting ``measure_g2_inu = True``.
It includes two sub-measurements, one in the particle-hole channel, and one in the particle-particle channel.

* Particle-hole channel, switched on by ``measure_g2_ph = True``.

    .. math::

        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega;\nu,\nu') =
        \frac{1}{\beta}\int_0^\beta d\tau_1d\tau_2d\tau_3d\tau_4\
        e^{-i\nu\tau_1} e^{i(\nu+\omega)\tau_2} e^{-i(\nu'+\omega)\tau_3} e^{i\nu'\tau_4}
        G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4).

    Accumulation result is available via ``G2_iw_inu_inup_ph`` solver attribute.

* Particle-particle channel, switched on by ``measure_g2_pp = True``.

    .. math::

        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega;\nu,\nu') =
        \frac{1}{\beta}\int_0^\beta d\tau_1d\tau_2d\tau_3d\tau_4\
        e^{-i\nu\tau_1} e^{i(\omega-\nu')\tau_2} e^{-i(\omega-\nu)\tau_3} e^{i\nu'\tau_4}
        G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4).

    Accumulation result is available via ``G2_iw_inu_inup_pp`` solver attribute.

Both sub-measurements can be enabled at the same time.
Numbers of bosonic and fermionic Matsubara frequencies are set by ``measure_g2_n_iw``
and ``measure_g2_n_inu`` parameters respectively.

One bosonic Matsubara frequency and two Legendre coefficients
*************************************************************

This measurement is activated by setting ``measure_g2_legendre = True``.
As in the previous paragraph, it includes two sub-measurements.

* Particle-hole channel, switched on by ``measure_g2_ph = True``.

    .. math::

        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega_m;\ell,\ell') \equiv \sum_{n,n'\in\mathbb{Z}}
        \bar T_{2n+m+1,\ell}
        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega_m;\nu_n,\nu_{n'})
        \bar T^*_{2n'+m+1,\ell'}.

    Accumulation result is available via ``G2_iw_l_lp_ph`` solver attribute.


* Particle-particle channel, switched on by ``measure_g2_pp = True``.

    .. math::

        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega_m;\ell,\ell') \equiv \sum_{n,n'\in\mathbb{Z}}
        \bar T_{2n+m+1,\ell}
        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega_m;\nu_n,\nu_{n'})
        \bar T^*_{2n'+m+1,\ell'}.

    Accumulation result is available via ``G2_iw_l_lp_pp`` solver attribute.

Here we have introduced transformation matrices from the Matsubara frequency domain to the
Legendre polynomial basis:

.. math::

    \bar T_{o,\ell} \equiv \frac{\sqrt{2\ell+1}}{\beta}
    \int_0^\beta d\tau e^{io\pi\frac{\tau}{\beta}} P_\ell[x(\tau)] =
    \sqrt{2\ell+1}i^o i^\ell j_\ell\left(\frac{o\pi}{2}\right).

Both sub-measurements can be enabled at the same time.
Numbers of bosonic Matsubara frequencies and Legendre coefficients are set by ``measure_g2_n_iw``
and ``measure_g2_n_l`` parameters respectively.


Impurity density matrix
-----------------------

The impurity density matrix (a.k.a. reduced density matrix) is the full density matrix of the system
with the bath degrees of freedom traced out.

.. math::

    \hat\rho_\mathrm{imp} = \mathrm{Tr}_\mathrm{bath}[e^{-\beta\hat H}/Z].

One can use this object to :ref:`estimate average values <static>`
of the static (:math:`\tau`-independent) impurity observables,

.. math::

    \langle\hat O\rangle = \mathrm{Tr}_\mathrm{at}[\hat O\hat\rho_\mathrm{imp}].

This measurement is activated by setting ``measure_density_matrix`` to ``True``. It also requires
enabling ``use_norm_as_weight`` parameter.

The impurity density matrix is accessible as ``density_matrix`` attribute of the solver object.

.. warning::
    Presently the density matrix is treated as block-diagonal with the same block structure as
    :math:`\hat H_\mathrm{loc}`. The block-offdiagonal matrix elements are not accumulated,
    so results can only be reliably used with static observables of the same block structure.

    The ``density_matrix`` attribute returns a list of matrices, one matrix per diagonal block.

Perturbation order histograms
-----------------------------

The perturbation order within a block :math:`A` is defined as a half of the number of
operators with the block index :math:`A` in the dynamical trace.
The total perturbation order is similarly related to the total number of operators in the dynamical trace.

Statistical histograms of the block-wise, as well as total perturbation orders will be measured if
``measure_pert_order`` is set to ``True``.

.. note::

    These two kinds of histograms are independent measurements. The total perturbation order histogram
    is expressed as a convolution of the block-wise histograms solely for the non-interacting systems.

For each block, the corresponding partial histogram is accessible as ``perturbation_order[block_name]``.
The ``perturbation_order_total`` attribute holds the total perturbation order histogram.

Average sign
------------

The average sign is defined as a ratio of two Monte-Carlo averages

.. math::

    \langle\mathrm{sign}\rangle = \frac
    {\langle\mathrm{sign}(W)(|\mathrm{Tr}_{at}[\ldots]|/W_{at})\rangle_{MC}}
    {\langle|\mathrm{Tr}_{at}[\ldots]|/W_{at}\rangle_{MC}},

where :math:`\mathrm{sign}(W)` is the sign of the total weight of a configuration,
and :math:`|\mathrm{Tr}_{at}[\ldots]|/W_{at}` is the atomic reweighting factor.

If ``use_norm_as_weight = False`` (no reweighting of the atomic problem), the reweighting
factor equals 1, and our definition of the average sign coincides with the usual one for
fermionic QMC algorithms. Otherwise, the denominator ensures the correct normalization
of the observable.

Result of this measurement is always available as ``average_sign`` attribute of the solver.
