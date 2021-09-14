.. highlight:: bash

.. _install:

Install CTHYB
*************

Packaged Versions of CTHyb
==========================

.. _ubuntu_debian:
Ubuntu Debian packages
----------------------

We provide a Debian package for the Ubuntu LTS Versions 16.04 (xenial) and 18.04 (bionic), which can be installed by following the steps outlined :ref:`here <triqslibs:triqs_debian>`, and the subsequent command::

        sudo apt-get install -y triqs_cthyb

.. _anaconda:
Anaconda (experimental)
-----------------------

We provide Linux and OSX packages for the `Anaconda <https://www.anaconda.com/>`_ distribution. The packages are provided through the `conda-forge <https://conda-forge.org/>`_ repositories. After `installing conda <https://docs.conda.io/en/latest/miniconda.html>`_ you can install CTHYB with::

        conda install -c conda-forge triqs_cthyb

See also `github.com/conda-forge/triqs_cthyb-feedstock <https://github.com/conda-forge/triqs_cthyb-feedstock/>`_.

.. _docker:
Docker
------

A Docker image including the latest version of CTHyb is available `here <https://hub.docker.com/r/flatironinstitute/triqs>`_. For more information, please see the page on :ref:`TRIQS Docker <triqslibs:triqs_docker>`.


Compiling CTHYB from source
===========================

.. note:: To guarantee reproducibility in scientific calculations we strongly recommend the use of a stable `release <https://github.com/TRIQS/triqs/releases>`_ of both TRIQS and its applications.

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

#. The NFFT3 library for non-uniform Fourier transformations https://www-user.tu-chemnitz.de/~potts/nfft/.
   To compile without NFFT3 the two-particle measurements has to be disabled by passing the cmake flag: ``-DMeasureG2=OFF``.
   
Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/cthyb`` repository from GitHub::

     $ git clone https://github.com/TRIQS/cthyb cthyb.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir cthyb.build && cd cthyb.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../cthyb.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``cthyb`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd cthyb.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.2.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``cthyb`` can be configured using CMake-options::

    cmake ../cthyb.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_cthyb          |
+-----------------------------------------------------------------+-----------------------------------------------+
| Allow the hybridization \Delta(tau) to be complex               | -DHybridisation_is_complex=ON                 |
+-----------------------------------------------------------------+-----------------------------------------------+
| Allow the local Hamiltonian H_loc to be complex                 | -DLocal_hamiltonian_is_complex=ON             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Measure the two particle object (requires the NFFT library)     | -DMeasureG2=ON                                |
+-----------------------------------------------------------------+-----------------------------------------------+
| Save visited configurations to configs.h5 (*developers only*)   | -DSAVE_CONFIGS=ON                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Enable extended debugging output (*developers only*)            | -DEXT_DEBUG=ON                                |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+

.. note::

    * Combination of options ``HYBRIDISATION_IS_COMPLEX=ON`` and ``LOCAL_HAMILTONIAN_IS_COMPLEX=OFF``
      is not supported.

    * The two-particle Green's function measurement requires the NFFT library. To build ``cthyb`` without NFFT pass ``-DMeasureG2=OFF`` to cmake.
