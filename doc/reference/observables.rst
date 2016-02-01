.. module:: pytriqs.applications.impurity_solvers.cthyb

.. _options:

Measurements: definitions
=========================

Here we list all the observables that can be measured by the solver along with their definition.
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

Impurity density matrix
-----------------------

The impurity density matrix (a.k.a. reduced density matrix) is the full density matrix of the system 
with the bath degrees of freedom traced out.

.. math::

    \hat\rho_\mathrm{imp} = \mathrm{Tr}_\mathrm{bath}[e^{-\beta\hat H}/Z].

One can use this object to estimate average values of the static (:math:`\tau`-independent)
impurity observables,

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

Average perturbation order
--------------------------

The perturbation order within a block :math:`A` is defined as a half of the number of
operators with the block index :math:`A` in the dynamical trace.
The total perturbation order is similarly related to the total number of operators in the dynamical trace.

Statistical histograms of the block-wise, as well as total perturbation orders will be measured if
``measure_pert_order`` is set to ``True``.

.. note::

    These two kinds of histograms are independent measurements. The total perturbation order histogram
    is expressed as a convolution of the block-wise histograms solely for the non-interacting systems.

.. warning::
    Update this section when the histograms are accessible through solver's attributes.

For each block, the corresponding partial histogram is dumped to a text file
`histo_pert_order_<block_name>.dat` upon completion of the simulation. The total histogram is
dumped to `histo_pert_order.dat`.

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