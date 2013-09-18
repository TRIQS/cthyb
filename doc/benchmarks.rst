
Benchmarks
==========

To benchmark the solver we run a multi-orbital test calculation called "5+5".
The Green's function is measured for a correlated impurity atom with five orbitals.
Each orbital (cubic harmonic) is hybridized with one bath site.

Parameters of the bath and local atomic levels are shown in the table below:

+--------------------+----------------+----------------+-----------------+----------------+---------------------+
|                    | :math:`d_{xy}` | :math:`d_{yz}` | :math:`d_{z^2}` | :math:`d_{xz}` | :math:`d_{x^2-y^2}` |
+====================+================+================+=================+================+=====================+
| :math:`V`          | 0.2            | 0.2            | 0.2             | 0.2            | 0.2                 |
+--------------------+----------------+----------------+-----------------+----------------+---------------------+
| :math:`\epsilon_k` | -0.2           | -0.15          | -0.1            | 0.05           | 0.4                 |
+--------------------+----------------+----------------+-----------------+----------------+---------------------+
| :math:`\epsilon_d` | -0.2           | -0.15          | -0.1            | 0.05           | 0.4                 |
+--------------------+----------------+----------------+-----------------+----------------+---------------------+

The interaction term of the atomic Hamiltonian is given by a full rotationally-invariant :math:`U`-matrix.
It is parametrized by three radial integrals :math:`F_0, F_2, F_4` or alternatively by constants :math:`U,J,\alpha` with relations 
:math:`U=F_0`, :math:`J=F_2(1+\alpha)/14`, :math:`F_4=\alpha F_2`.
Those parameters are chosen as :math:`U=4.0, J=0.7, \alpha=0.63`.

The chemical potential is set to :math:`\mu=26` (7 electrons in the ground state) and the inverse temperature is :math:`\beta=40`.

Only those boundary states are included into the trace, for which :math:`\exp(-\beta E_{at})/Z_{at} > 10^{-10}` (16 states).

5-band results
--------------

Here must be a plot of a :math:`G(\tau)`. But it is still too noisy. Such a shame...


Histogram of subspace dimensions
--------------------------------

A typical distribution of subspace dimensions in one "5+5" run.

.. image:: dims.stats.5+5.png

This histogram is produced for

- ``length_cycle = 50``
- ``n_warmup_cycles = 300``
- ``n_cycles = 10000``

The hybrid mode
---------------

The solver can be run in a hybrid mode, when two different algorithms are used
to calculate :math:`\exp(-\tau \hat H)` for smaller and larger dimensions of considered subspaces.
This mode is controlled by input parameter ``krylov_small_matrix_size``. For matrices larger
than ``krylov_small_matrix_size`` the Krylov subspace reduction algorithm will be used to
construct the evolution operator. Otherwise the exponential will be computed as a product
of 3 matrices (one diagonal and two precomputed unitary matrices).

Here are some performance benchmark data obtained for different values of ``krylov_small_matrix_size``
(``n_warmup_cycles = 100``, ``n_cycles = 3000``).


- "5+5" test, quantum numbers :math:`N_{\uparrow}` and :math:`N_{\downarrow}`

+------------------------------+--------------------------+
| ``krylov_small_matrix_size`` | duration, seconds        |
+==============================+==========================+
| 100                          | 591                      |
+------------------------------+--------------------------+
| 50                           | 999                      |
+------------------------------+--------------------------+
| 25                           | 2450                     |
+------------------------------+--------------------------+
| 10                           | 2559                     |
+------------------------------+--------------------------+
| 5                            | 2666                     |
+------------------------------+--------------------------+
| 1                            | 2660                     |
+------------------------------+--------------------------+
| 0                            | 2668                     |
+------------------------------+--------------------------+
      
- "5+5" test, no quantum numbers

+------------------------------+--------------------------+
| ``krylov_small_matrix_size`` | duration, seconds        |
+==============================+==========================+
| 1024                         | at least 4 days          |
+------------------------------+--------------------------+
| 0                            | 66984                    |
+------------------------------+--------------------------+

- "5+5" test, cthyb_matrix, quantum numbers :math:`N_{\uparrow}` and :math:`N_{\downarrow}`

One run of the old solver for the same number of steps lasted in average 2871 seconds.