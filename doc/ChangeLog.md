Version 3.0.0
=============

Fundamental Green Function Symmetries
-------------------------------------

We now enforce the fundamental Green function properties `G[iw](i,j) = G[-iw](j,i)*`
and `G[tau](i,j) = G[tau](j,i)*` for both the solver-input S.G0_iw as well as the
solver outputs S.G_tau and S.G_iw.
If the input violates the relation a warning will be issued and the input will
be automatically symmetrized. The output will always be symmetrized
in the collect_results section of the G_tau measurement.


Version 2.2.1
=============

CTHYB Version 2.2.1 makes the application available
through the Anaconda package manager. We adjust
the install pages of the documentation accordingly.
It further introduces two minor fixes.
We provide a more detailed description of the changes below.

doc
---
* Add a section on the Anaconda package to the install page

General
-------
* Correct error in solve params deprecation warning Fix #124
* Remove redundant inclusion of triqs/utility/serialization.hpp


Version 2.2.0
=============

CTHYB version 2.2.0 is a compatibility release
for TRIQS version 2.2.0. It provides improvements to
the documentation and fixes various compiler warnings.

We provide a more detailed description of the changes below.

cmake
-----
* Install additional triqs namespace headers provided by cthyb
* Synchronize project structure with app4triqs
* Bump cthyb version to 2.2.0 and adjust triqs version requirement
* Minor adjustments in cmake to restore cmake version 3.0.2 compatibility

doc
---
* Update install instructions for nfft
* Install svg files used by tutorials
* Add tutorial for pert order histogram
* Add tutoria for dynamic suscept
* Fix import statements -> static observables (import from atom_diag)

General
-------
* Use std::vector over triqs::arrays::vector for density containers
* Fix various compiler warnings and deprecation messages
* Consistently use std:: namespace instead of std::c14::
* Remove remainders of separate cpp2py install from modulefile
* Bump cmake version requirement to 3.0.2
* Port cthyb after latest changes to triqs mpi/itertools functionality
* Add measurement O_tau_min_ins

h5
--
* Store additional member variables when storing solver object, FIX #115

packaging
---------
* update package name to triqs_cthyb and conflict with cthyb


Version 2.1.0
-------------

Contributors
~~~~~~~~~~~~

* **Hugo U. R. Strand**, new functionality, two-particle Green's function sampling in imaginary time and frequency, and sampling by operator insertion.
* **Igor Krivenko**, occupation number constraints in the trace, two-particle Green's function sampling with NFFT and mixed basis frequency*(legendre)^2.
* **Lewin Boehnke**, (indirect contribution) through the first mixed basis sampling implementation in an early version of ``cthyb``.
* **Nils Wentzell**, nfft_buf_t, cmake, and install.

Testers
~~~~~~~

Manuel Zingl, Daniel Mantadakis, and Jonathan Karp.

Changes in installation and cmake files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Headers are now installed to ``${CMAKE_INSTALL_PREFIX}/include/triqs_cthyb``.
* The c++ namespace is now a single namespace called ``triqs_cthyb``.
* The python module is named ``triqs_cthyb``.

Solver Interface
~~~~~~~~~~~~~~~~

* The type of gf_struct has changed from a dict to a list of pairs

  ``{ 'up' : [0, 1], 'dn' : [0, 1] }    -->   [ ('up', [0,1]), ('dn', [0,1]) ]``

* measure_g_l has been renamed to measure_G_l
* measure_g_tau has been renamed to measure_G_tau

* Truncation of the local Hilbert space in total density using ``loc_n_min`` and ``loc_n_max``
* Two particle Greens function sampling with a whole set of measurements see ``measure_G2_*``
* Two particle susceptibilities sampling by insertion using `measure_O_tau`
* Tunable parameters for the determinant regularization see `det_*` parameters
* Solver object can be stored and loaded directly to/from h5 archive
* Tail fit functionality is rewritten for the Triqs 2.0 change in how tail coefficients are handled

Changes in behavior
~~~~~~~~~~~~~~~~~~~

* ``move_double`` is now enabled by default.
* The sampled single particle Greens function is modified when tail fitting is enabled.
* The ``gf_struct`` should be a list of list of block index and a list of target space indices, e.g.

  ``gf_struct = [['up', [0, 1]], ['down', [0, 1]]]``.

  When passing a dictionary a warning is printed.

Dependencies
~~~~~~~~~~~~

* The NFFT3 library by default, to compile without the library use the cmake flag ``-DMeasureG2=OFF``.
