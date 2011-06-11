=====================================
 Amanzi Extern Software Requirements
=====================================

System Software Environment
=====================================

The libraries and tools listed in this section are assumed to be part
of the target platform. Typically, a user would not build these
libraries from source.  The Amanzi team recommends that users use
pre-built binaries or tools such as MacPorts, yum, rpm to install
these libraries.

At a minimum the build environment must have the following system
software libraries and tools present.

* A UNIX-like OS. Amanzi has been built on several flavors of Linux
  (Red Hat/CentOS, Ubuntu, Fedora) and Mac OS 10.5 and 10.6. Native
  Windows builds are not supported at this time.
* C and C++ compiler.
* Python 2.6 version or higher.
* MPI. Amanzi has been built against OpenMPI and MPICH. Amanzi does
  not support non-MPI builds at this time.
* BLAS/LAPACK built and tuned to the target architecture.

Required External Software
=====================================

Amanzi leverages existing software libraries whenever possible. Amanzi
can not built until all the required libraries are built
successfully. In some cases (netCDF, Exodus, HDF5 and Trilinos), the
libraries must be configured in a particular way to work correctly
with Amanzi and other TPLs. The Amanzi source directory contains a
script that automatically builds the TPLs correctly. Users are
strongly encouraged to use this tool to build the TPLs. The current
documentation for a correct TPL configuration is tracked in this
tool. Advanced users creating their own build process should refer to
this script for supported configurations and TPL build order.

Amanzi builds the source located in
https://software.lanl.gov/ascem/trac/browser every hour on several OS
types.  The build and test suite results for each OS is recorded and
tracked on the ASCEM Trac site at
https://software.lanl.gov/ascem/trac/buildbot.

