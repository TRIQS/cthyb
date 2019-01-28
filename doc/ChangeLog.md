Version 2.1
-----------



Version 1.5
-----------

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
