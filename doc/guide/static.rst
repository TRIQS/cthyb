.. _static:

Measuring static observables
============================

In addition to the interacting Green's functions one is often
interested in such observables, as orbital occupation numbers and
double occupancy on the impurity. Expectation values of these
time-independent operators are readily expressed in terms of
the reduced density matrix of the system,

.. math::

    \hat\rho_\mathrm{imp} = \mathrm{Tr}_\mathrm{bath}[e^{-\beta\hat H}/Z].

Here, :math:`e^{-\beta\hat H}/Z` is the density matrix of the full system
(impurity + bath) described by the Hamiltonian :math:`\hat H`.

Let us consider a single-orbital Anderson impurity in a weak
external magnetic field [:download:`script <static.py>`].

At first, we import all necessary modules, define some input parameters
and construct a ``Solver`` object in the usual way.

.. literalinclude:: static.py
    :lines: 1-21

We instruct the ``solve`` function to accumulate the density matrix by passing
``measure_density_matrix = True`` and ``use_norm_as_weight = True``. The
latter parameter tells the solver to employ a reweighting scheme -- use a
spherical norm instead of the trace to calculate the atomic weight.
The reweighting allows to take into account important contributions to
:math:`\hat\rho_\mathrm{imp}`, which otherwise would be missed.

.. literalinclude:: static.py
    :lines: 23-29

Results of density matrix accumulation are accessible via `density_matrix` attribute.
We also need information about the structure of the local Hilbert space to compute
expectation values of static observables. This information is stored as a special
object in `h_loc_diagonalization`.

.. literalinclude:: static.py
    :lines: 31-35

According to a well known formula, an expectation value of an operator :math:`\hat O`
acting on the impurity degrees of freedom is given by

.. math::

    \langle\hat O\rangle = \mathrm{Tr}_\mathrm{at}[\hat\rho_\mathrm{imp} \hat O].

There is a function called `trace_rho_op()`, which does this trace and returns the
expectation value of :math:`\hat O`.

.. literalinclude:: static.py
    :lines: 37-44

Typical output of the script may look like::

    <N_up> = 0.620679103675
    <N_down> =  0.379442861526
    <N_up*N_down> = 0.00369097179335

