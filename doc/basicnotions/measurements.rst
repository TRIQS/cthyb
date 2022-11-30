.. _measurements:

Measurements: definitions
=========================

Here we list all the observables that can be measured by the solver along with their definitions.
Each measurement can be turned on or off via the corresponding :ref:`solve() parameters <solve_parameters>`.

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

This measurement is turned on by setting ``measure_G_tau`` to ``True``.
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

This measurement is controlled through the switch ``measure_G_l``.
The result of accumulation is accessible as ``G_l`` attribute of the solver object.
The number of Legendre coefficients to be measured is specified through constructor's parameter ``n_l``.

Local susceptibility and other operator pairs
---------------------------------------------

Higher order response functions like local susceptibilities

.. math::
   \chi_{\hat{O}_1, \hat{O}_2}(\tau) \equiv \langle \hat{O}_1(\tau) \hat{O}_2 \rangle

can be sampled (by insertion) in the trace. This approach is limited to operators :math:`\hat{O}_1` and :math:`\hat{O}_2` that commutes with each other and the local Hamiltonian, :math:`[\hat{O}_1, \hat{O}_2] = 0`, :math:`[\hat{O}_i, H_{loc}]=0`.

This measurement is controlled by the argument ``measure_O_tau``. To enable the measurement pass a tuple of the operators to be sampled, e.g. ``measure_O_tau = (n('up',0), n('do',0))`` will measure the response function :math:`\langle \hat{n}_\uparrow(\tau) \hat{n}_\downarrow \rangle`. The resulting response function is accessible as the ``O_tau`` attribute of the solver object. The number of operator insertions is by default taken to be the square of the perturbation order of each configuration, however, for cases with perturbation order lower than the parameter ``measure_O_tau_min_ins`` the number of insertions is kept fixed to this minimum value, the default value is a minimum of 10 insertions.

Two-particle Green's functions
------------------------------

CTHYB implements several ways to measure the two-particle Green's functions.
We define :math:`G^{(2)}` in the imaginary time as

.. math::

    G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4) =
    \langle\mathcal{T}_\tau c^\dagger_\alpha(\tau_1) c_\beta(\tau_2) c^\dagger_\gamma(\tau_3) c_\delta(\tau_4)\rangle.
    
Depending on the value of the ``measure_G2_block_order`` parameter the combined block/inner indices
:math:`\alpha,\beta,\gamma,\delta` are interpreted differently.

* ``measure_g2_block_order = 'AABB'``: :math:`\{\alpha,\beta,\gamma,\delta\} = \{(A,a),(A,b),(B,c),(B,d)\}`;
* ``measure_g2_block_order = 'ABBA'``: :math:`\{\alpha,\beta,\gamma,\delta\} = \{(A,a),(B,b),(B,c),(A,d)\}`.

These two combinations exhaust all the possibilities for non-vanishing contributions to :math:`G^{(2)}`.

Measuring :math:`G^{(2)}` is in general expensive, and in some cases time can be saved by skipping
particular block pairs :math:`(A,B)`. Block pairs to be measured are selected via the ``measure_G2_blocks``
parameter (all possible pairs are selected by default):

``measure_G2_blocks = set([("A1","B1"), ("A2","B2"), ...])``

See the `PhD thesis of L. Boehnke <http://ediss.sub.uni-hamburg.de/volltexte/2015/7325/pdf/Dissertation.pdf>`_
for an in-depth discussion of these measurements.

Imaginary time binning
**********************

The simplest measure is to bin :math:`G^{(2)}` directly in three imaginary times, after using time-tranlsational invariance to fix :math:`\tau_4 = 0`.

* The imaginary time binning is switched on by ``measure_G2_tau = True``.

  .. math::

     G^{(2)}_{\alpha\beta\gamma\delta}( \tau_1, \tau_2, \tau_3) \equiv
     G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3, 0)

  Accumulation result is available via the ``G2_tau`` solver attribute.

The number of bins in each imaginary time dimension is set by ``measure_G2_n_tau``. 

Matsubara frequency measurements
********************************

In Matsubara frequency space there are three different measurements:

* Three Fermionic Matsubara frequency measurement, switched on by ``measure_G2_iw = True``.

    .. math::

        G^{(2)}_{\alpha\beta\gamma\delta}(\nu_1, \nu_2, \nu_3) =
        \int_0^\beta d\tau_1 d\tau_2 d\tau_3 \,
        e^{-i\nu_1 \tau_1 + i\nu_2 \tau_2 - i\nu_3 \tau_3}
        G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,0) \, .

  Accumulation result is available via the ``G2_inu`` solver attribute.
  
* Particle-hole channel measurement,
  with one Bosonic (:math:`\omega`) and two Fermionic (:math:`\nu, \nu'`) frequencies,
  switched on by ``measure_G2_iw_ph = True``.

    .. math::

        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega;\nu,\nu') =
        \frac{1}{\beta}\int_0^\beta d\tau_1d\tau_2d\tau_3d\tau_4 \,
        e^{-i\nu\tau_1 + i(\nu+\omega)\tau_2 - i(\nu'+\omega)\tau_3 + i\nu'\tau_4}
        G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4) \, .

  Accumulation result is available via ``G2_iw_ph`` solver attribute.

* Particle-particle channel measurement,
  with one Bosonic (:math:`\Omega`) and two Fermionic (:math:`\nu, \nu'`) frequencies,
  switched on by ``measure_G2_iw_pp = True``.

    .. math::

        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega;\nu,\nu') =
        \frac{1}{\beta}\int_0^\beta d\tau_1d\tau_2d\tau_3d\tau_4\
        e^{-i\nu\tau_1 + i(\omega-\nu')\tau_2 - i(\omega-\nu)\tau_3 + i\nu'\tau_4}
        G^{(2)}_{\alpha\beta\gamma\delta}(\tau_1,\tau_2,\tau_3,\tau_4).

  Accumulation result is available via ``G2_iw_pp`` solver attribute.

The number of Bosonic and Fermionic Matsubara frequencies are set by the ``measure_G2_n_bosonic``
and ``measure_G2_n_fermionic`` parameters respectively.

All these frequency measurements use direct evaluation in frequency space, using an optimized frequency spreader that avoids the re-evaluation of exponential functions. This implementation performs well for low perturbation order and small number of sampled Matsubara frequencies.

As an alternative all frequency measurements are also implemented using non-equidistant fast fourier transform (NFFT) to speed up the sampling procedure. The corresponding flags and attributes are:

* (``measure_G2_iw_nfft``, ``G2_iw_nfft``),

* (``measure_G2_iw_ph_nfft``, ``G2_iw_ph_nfft``), and

* (``measure_G2_pp_nfft``, ``G2_iw_ph_nfft``).

Depending on the impurity model the NFFT buffer can be adjusted for maximum performance by setting ``nfft_buf_sizes``.

Whether the direct frequency evaluation or NFFT performs better is problem dependent and has to be tested case by case.

Mixed Matsubara Frequency and Legendre measurements
***************************************************

* Particle-hole channel, switched on by ``measure_G2_iwll_ph = True``.

    .. math::

        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega_m;\ell,\ell') \equiv \sum_{n,n'\in\mathbb{Z}}
        \bar T_{2n+m+1,\ell}
        G^{(2)ph}_{\alpha\beta\gamma\delta}(\omega_m;\nu_n,\nu_{n'})
        \bar T^*_{2n'+m+1,\ell'}.

    Accumulation result is available via ``G2_iwll_ph`` solver attribute.


* Particle-particle channel, switched on by ``measure_G2_iwll_pp = True``.

    .. math::

        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega_m;\ell,\ell') \equiv \sum_{n,n'\in\mathbb{Z}}
        \bar T_{2n+m+1,\ell}
        G^{(2)pp}_{\alpha\beta\gamma\delta}(\omega_m;\nu_n,\nu_{n'})
        \bar T^*_{2n'+m+1,\ell'}.

    Accumulation result is available via ``G2_iwll_pp`` solver attribute.

Numbers of bosonic Matsubara frequencies and Legendre coefficients are set by the ``measure_G2_n_iw``
and ``measure_G2_n_l`` parameters respectively.

The one Bosonic Matsubara frequency is treated using non-equidistant fast fourier tranform (NFFT) and the NFFT buffer size is set by ``measure_G2_iwll_nfft_buf_size``.
    
The transformation matrices :math:`\bar{T}_{o, \ell}` introduced above transforms from the Matsubara frequency domain to the Legendre polynomial basis:

.. math::

    \bar T_{o,\ell} \equiv \frac{\sqrt{2\ell+1}}{\beta}
    \int_0^\beta d\tau e^{io\pi\frac{\tau}{\beta}} P_\ell[x(\tau)] =
    \sqrt{2\ell+1}i^o i^\ell j_\ell\left(\frac{o\pi}{2}\right).


Impurity density matrix
-----------------------

The impurity density matrix (a.k.a. reduced density matrix) is the full density matrix of the system
with the bath degrees of freedom traced out.

.. math::

    \hat\rho_\mathrm{imp} = \mathrm{Tr}_\mathrm{bath}[e^{-\beta\hat H}/Z].

One can use this object to `estimate average values <../guide/static_observables_notebook.ipynb>`_
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
