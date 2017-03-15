.. highlight:: bash

.. _install:

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources of the solver from github::

     $ git clone https://github.com/TRIQS/cthyb.git cthyb.src

#. Create an empty build directory where you will compile the code::

     $ mkdir cthyb.build && cd cthyb.build

#. In the build directory call cmake specifying where the TRIQS library is installed::

     $ cmake -DTRIQS_PATH=path_to_triqs ../cthyb.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

.. note:: Be careful with the cmake command above: set TRIQS_PATH, not CMAKE_INSTALL_PREFIX (this variable is only for the TRIQS library)!

Version compatibility
---------------------

Be careful that the version of the TRIQS library and of the solver must be
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
If you want to use a version of
the solver that is not the latest one, go into the directory with the sources
and look at all available versions::

     $ cd cthyb.src && git tag

Checkout the version of the code that you want::

     $ git checkout 1.0.0

Then follow the steps 2 to 4 described above to compile the code.

Custom CMake options
--------------------

Functionality of ``cthyb`` can be tweaked using extra compile-time options passed to CMake::

    cmake -DTRIQS_PATH=path_to_triqs -DOPTION1=value1 -DOPTION2=value2 ... ../cthyb.src

+---------------------------------------------------------------+-----------------------------------+
| Options                                                       | Syntax                            |
+===============================================================+===================================+
| Disable testing (not recommended)                             | -DTests=OFF                       |
+---------------------------------------------------------------+-----------------------------------+
| Build the documentation locally                               | -DBUILD_DOC=ON                    |
+---------------------------------------------------------------+-----------------------------------+
| Allow the hybridization \Delta(tau) to be complex             | -DHYBRIDISATION_IS_COMPLEX=ON     |
+---------------------------------------------------------------+-----------------------------------+
| Allow the local Hamiltonian H_loc to be complex               | -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON |
+---------------------------------------------------------------+-----------------------------------+
| Enable extended debugging output (*developers only*)          | -DEXT_DEBUG=ON                    |
+---------------------------------------------------------------+-----------------------------------+
| Save visited configurations to configs.h5 (*developers only*) | -DSAVE_CONFIGS=ON                 |
+---------------------------------------------------------------+-----------------------------------+

.. note::

    * Combination of options ``HYBRIDISATION_IS_COMPLEX=ON`` and ``LOCAL_HAMILTONIAN_IS_COMPLEX=OFF``
      is not supported.

    * The two-particle Green's function measurement requires the TRIQS library to be built with NFFT support.
