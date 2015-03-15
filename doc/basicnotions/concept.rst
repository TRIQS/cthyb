.. _solver_concept:

The solver concept
===================

A solver class must obey a specific concept so that inputs don't depend on
which specific solver you are using. The solver class must contain the members
`G_tau` and `G0_iw` corresponding to the interacting and non-interacting Green's
function, and a member function `solve()` that uses `G0_iw` to compute `G_iw`.
Here's the concept explicitly:

.. class:: GenericSolver

  The generic concept of a solver class

  .. attribute:: G_tau

    The interacting Green's function. This is a :ref:`full Green's function
    <fullgreen>` with each block being a :ref:`GfImTime`. G is set after
    the solve() method has been called.

  .. attribute:: G0_iw

    The non-interacting Green's function. This is a :ref:`full Green's function
    <fullgreen>` with each block being a :ref:`GfImFreq`. G0_iw is needed
    by the solver and should be set before solve() is called.

  .. function:: solve()

    This function computes the interacting Green's function of the problem
    using G0_iw as an input. The solution is stored in G_iw.
