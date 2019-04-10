Version 2.1
-----------

Contributors
~~~~~~~~~~~~

Hugo U. R. Strand, new functionality, G2 sampl, sampl by ins.
Nils Wentzell, cmake and install

Changes in installation and cmake files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Headers are now installed to ${CMAKE_INSTALL_PREFIX}/include/triqs_cthyb
* The c++ namespace is now a single namespace called triqs_cthyb

Solver Interface
~~~~~~~~~~~~~~~~

* The type of gf_struct has changed from a dict to a list of pairs
  { 'up' : [0, 1], 'dn' : [0, 1] }    -->   [ ('up', [0,1]), ('dn', [0,1]) ]
* measure_g_l has been renamed to measure_G_l
* measure_g_tau has been renamed to measure_G_tau

* Truncation of the local Hilbert space in total density using `loc_n_min` and `loc_n_max`
* Two particle Greens function sampling with a whole set of measurements see `measure_G2_*`
* Two particle susceptibilities sampling by insertion using `measure_O_tau`
* Tunable parameters for the determinant regularization see `det_*` parameters
* Solver object can be stored and loaded directly to/from h5 archive
* Tail fit functionality is rewritten for the Triqs 2.0 change in how tail coefficients are handled

Changes in behavior
~~~~~~~~~~~~~~~~~~~

* `move_double` is now enabled by default.
* The sampled single particle Greens function is modified when tail fitting is enabled.
* The `gf_struct` should be a list of list of block index and a list of target space indices, e.g. `gf_struct = [['up', [0, 1]], ['down', [0, 1]]]`. When passing a dictionary a warning is printed.

Quality Assurance
~~~~~~~~~~~~~~~~~

* Checking the benchmarks TODO (HS)
* Checking interoperability with triqs/DFTTools (MZ, AG) TODO
* Testing hermicity of the tail fit, Jonathan Karp, Columbia University
