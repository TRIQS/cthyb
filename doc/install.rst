.. highlight:: bash

.. _install:


Packaged Versions of CTHyb
==========================

A Docker image including the latest version of CTHyb is available `here <https://hub.docker.com/r/flatironinstitute/triqs>`_. For more information, please see the page on :ref:`TRIQS Docker <triqslibs:triqs_docker>`.

We further provide a Debian package for the Ubuntu LTS Versions 16.04 (xenial) and 18.04 (bionic), which can be installed by following the steps outlined :ref:`here <triqslibs:triqs_debian>`, and the subsequent command::

        sudo apt-get install -y cthyb


Compiling CTHyb from source
===========================

Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox and Cpp2Py (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that both are installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources of the solver from github::

     $ git clone https://github.com/TRIQS/cthyb.git cthyb.src

#. Create an empty build directory where you will compile the code::

     $ mkdir cthyb.build && cd cthyb.build

#. Make sure that you have added the TRIQS and Cpp2Py installation to your environment variables::

     $ source path_to_triqs/share/cpp2pyvarsh.sh
     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake::

     $ cmake ../cthyb.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

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

    cmake -DOPTION1=value1 -DOPTION2=value2 ... ../cthyb.src

+---------------------------------------------------------------+-----------------------------------------------+
| Options                                                       | Syntax                                        |
+===============================================================+===============================================+
| Disable testing (not recommended)                             | -DBuild_Tests=OFF                             |
+---------------------------------------------------------------+-----------------------------------------------+
| Specify an installation path other than path_to_triqs         | -DCMAKE_INSTALL_PREFIX=path_to_cthyb          |
+---------------------------------------------------------------+-----------------------------------------------+
| Build the documentation locally                               | -DBuild_Documentation=ON                      |
+---------------------------------------------------------------+-----------------------------------------------+
| Allow the hybridization \Delta(tau) to be complex             | -DHybridisation_is_complex=ON                 |
+---------------------------------------------------------------+-----------------------------------------------+
| Allow the local Hamiltonian H_loc to be complex               | -DLocal_hamiltonian_is_complex=ON             |
+---------------------------------------------------------------+-----------------------------------------------+
| Measure the two particle object (requires the NFFT library)   | -DMeasureG2=ON                                |
+---------------------------------------------------------------+-----------------------------------------------+
| Enable extended debugging output (*developers only*)          | -DEXT_DEBUG=ON                                |
+---------------------------------------------------------------+-----------------------------------------------+
| Save visited configurations to configs.h5 (*developers only*) | -DSAVE_CONFIGS=ON                             |
+---------------------------------------------------------------+-----------------------------------------------+

.. note::

    * Combination of options ``HYBRIDISATION_IS_COMPLEX=ON`` and ``LOCAL_HAMILTONIAN_IS_COMPLEX=OFF``
      is not supported.

    * The two-particle Green's function measurement requires the TRIQS library to be built with NFFT support.
