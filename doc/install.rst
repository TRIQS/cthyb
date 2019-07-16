.. highlight:: bash

.. _install:

Compiling triqs_cthyb from source
===============================


Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/triqs_cthyb`` repository from GitHub::

     $ git clone https://github.com/TRIQS/triqs_cthyb triqs_cthyb.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir triqs_cthyb.build && cd triqs_cthyb.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../triqs_cthyb.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``triqs_cthyb`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd triqs_cthyb.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``triqs_cthyb`` can be configured using CMake-options::

    cmake ../triqs_cthyb.src -DOPTION1=value1 -DOPTION2=value2 ... ../triqs_cthyb.src

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_triqs_cthyb      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
