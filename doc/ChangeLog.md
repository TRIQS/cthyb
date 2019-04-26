Version 2.2.0
-------------

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
