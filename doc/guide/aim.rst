.. _aim:

An example: the Anderson impurity model
=======================================

To illustrate how the CTQMC solver works in practice, we show the example of a
one-orbital Anderson impurity embedded in a flat (Wilson) conduction bath. The
interacting part of the local Hamiltonian of this problem is simply:

.. math::

  \mathcal{H}_\mathrm{int} = U n_\uparrow n_\downarrow,

and the non-interacting Green's function is:

.. math::

  G^{-1}_{0,\sigma} (i \omega_n) = i \omega_n - \epsilon_f - V^2 \Gamma_\sigma(i \omega_n).

In this example, there is a Coulomb interaction :math:`U` on the impurity level,
which is at an energy :math:`\epsilon_f`. The bath Green's function is :math:`\Gamma(i
\omega_n)`, and it has a flat density of states over the interval
:math:`[-1,1]`.  Finally, :math:`V` is the hybridization amplitude between the
impurity and the bath. Let us solve this problem with the CTQMC solver. Here is
the python :download:`script <aim.py>`:

.. literalinclude:: aim.py

Running this script on a single processor takes about 5 minutes and generates
an HDF5 archive file called :file:`aim_solution.h5`. This file contains the Green's
function in imaginary time, in imaginary frequencies, and in Legendre polynomial basis
found by the solver.
Let us plot the Green's function:

.. plot:: guide/aim_plot.py
   :include-source:
   :scale: 70

As expected the result shows a particle-hole symmetric impurity Green's
function (the real part vanishes up to the statistical noise).

Let us now go through the script in some more details.

.. literalinclude:: aim.py
  :lines: 1-5

These lines import the classes to manipulate Green's functions, fermionic
operators, and the CTQMC solver.

.. literalinclude:: aim.py
  :lines: 7-9

This just sets the parameters of the problem.

.. literalinclude:: aim.py
  :lines: 11-13

This is the construction of the Solver object. The class is described
in more detail in the section :ref:`ctqmc_ref`. Basically, the constructor
of the Solver needs two keywords:

- ``beta``: the inverse temperature,
- ``gf_struct``: a list of pairs [ (str, int), ...] describing the block structure of the Green's function.

This ensures that all quantities within Solver are correctly initialised,
and in particular that the block structure of all Green's function objects is consistent.

After the Solver is constructed it needs to know what the non-interacting Green's function
of the impurity is. From this information, the Solver will deduce the hybridization function
which is used in the algorithm. The non-interacting Green's function must be put in the
class member ``S.G0_iw``:

.. literalinclude:: aim.py
  :lines: 15-16

At this stage, everything is ready for the Solver and we just run it calling its member
function ``solve``:

.. literalinclude:: aim.py
  :lines: 18-23

The run is controlled by the parameters of ``solve()``:

- ``h_int``: The interacting part of the local Hamiltonian written with :ref:`TRIQS operators <triqslibs:operators>`.
- ``n_cycles``: The number of Monte Carlo cycles.
- ``length_cycle``: The number of Monte Carlo moves in a cycle.
- ``random_name``: The name of the random number generator.
- ``measure_g_l``: We want to accumulate the Green's function in a basis of Legendre polynomials.

When ``solve()`` has finished, it puts the result for the interacting Green's
function in imaginary time in its member ``S.G_tau`` and the Green's function
in imaginary frequencies in ``S.G_iw``. The last lines of the script save the
Green's function in the HDF archive.
