.. _bootstrap:

=============================================
Quickstart Using Bootstrap
=============================================

Introduction
------------

The Amanzi code is built using several required and optional Third Party Libraries (TPLs).  The script ``bootstrap.sh`` located in the ``config`` directory is designed to build most of the TPLs and the Amanzi source code on most UNIX-like OS (including Mac OSX).  It will also execute the included test suite.  A full set of options and defaults using ``<amanzi source root>/config/bootstrap.sh --help``.

The following pacakges are required before executing ``bootstrap.sh``

    * A recent and stable GCC, Clang, or Intel compiler. We have successfully built with
      GNU versions >= 7.5, Clang >= 5, and Intel >= 19.
    * Git version 1.8 or higher.
    * OpenSSL (required to build CURL) installation. This is available on
      most UNIX flavored systems and Macs.
    * BLAS/LAPACK built and tuned to the target architecture. See the 
      Trilinos Configuration section for more information.
    * A location to install all the TPLs. The entire software stack
      will be installed in directories rooted to this directory. Since the
      download--patch-build-install cycle for each TPL is tied to a single  
      target, the user must have read and write permission for this location.
    * The build directory will need approximately 2.1 Gb of space (for staticlly linked build)
      and the install directory tree requires 300 Mb.

The Amanzi build system utilizes CMake (version 3.17.0 or higher). 
If ``boostrap.sh`` does not find the cmake executable, it will download and install CMake 3.17.0.

Amanzi is a parallel project; therefore, a working MPI installation is required.  MPI installation directory and compiler wrapper locations can be specified through ``bootstrap.sh`` options.  If the MPI location is not specified the script will search common locations.  If working mpi wrappers are not found, the script will download and install OpenMPI 3.1.4.  See the :ref:`MPI` section on how to obtain and install other mpi implementations.

At this time, the script does not build compilers, ``LAPACK``, or ``BLAS``.  See the :ref:`compilers` section on how to obtaion a supported compiler and the :ref:`LAPACK` section for ``LAPACK`` and ``BLAS`` suggestions.


Bootstrap.sh Execution
--------------------------

After obtaining the Amanzi source code, the TPLs and source can be built utilizing the script ``bootstrap.sh`` located in ``<amanzi source root>/config/``.  The script will create two directories in the current working directory.  The first is ``TPL_BUILD`` for downloading and building the TPLs.  The second is ``amanzi-build``, where the Amanzi source code will be built and the test suite executed, if enabled.  The TPLs and Amanzi will be installed in the location specified by the option ``--prefix``, $HOME/amanzi by default.  All options and defaults can be printed using the ``--help`` option.


Building TPLs and Amanzi
++++++++++++++++++++++++

Choose a location to build the TPLs and Amanzi.  From that location execute ``<amanzi source root>/config/bootstrap.sh``.  The full list of options and defaults is printed using the ``--help`` option.

Below are a selection of common and useful options.

Specify TPL and Amanzi installation directory::

    --prefix=<install dir>

Specify MPI install location::

    --with-mpi=<mpi install dir>

Specify MPI compilers::

    --with-c-compiler=<mpi install dir>/bin/mpicc
    --with-cxx-compiler=<mpi install dir>/bin/mpicxx
    --with-fort-compiler=<mpi install dir>/bin/mpif90

Specify processor count for parallel building::

    --parallel=n

Enable/Disable capabilities, mesh toolkits, and the test suite::

    --enable-structured or disable-structured
    --enable-unstructured or disable-unstructured

Disabling the test suite is recommended on HPC clusters.  The test suite can be manually run on commute nodes using ``make test`` in the ``amanzi-build`` directory.::

    --disable-test_suite

Building Amanzi Only
++++++++++++++++++++

TPLs are updated less frequently than the Amanzi source code.  Therefore, users may wish to build Amanzi without rebuilding the TPLs.  The ``bootstrap.sh`` script installs the CMake configuration file in the install directory under tpls/share/cmake.  To skip rebuilding the TPLs use the following option.::

  --tpl-config-file=<install dir>/tpls/share/cmake/amanzi-tpl-config.cmake

Example ``bootstrap.sh`` Command Line
+++++++++++++++++++++++++++++++++++++

The following example command line is from a Mac OS X system using openmpi installed through MacPorts::

  sh <amanzi_source_path>/config/bootstrap.sh --parallel=4 \
                                              --with-c-compiler=/opt/local/bin/openmpicc \
                                              --with-cxx-compiler=/opt/local/bin/openmpicxx \
                                              --with-fort-compiler=/opt/local/bin/openmpif90 \
                                              --prefix=<ascem_path>/install \
                                              --with-cmake=/usr/bin/cmake \
                                              --with-mpi=/opt/local \
                                              --enable-test_suite \
                                              --enable-geochemistry

