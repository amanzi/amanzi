## Getting Started:  Building / Installing / Running

The User Guide is available at

  https://amanzi.github.io/amanzi/UserGuide

and provides the easiest way to click through what you would like to
read about. Nevertheless, all of the building/installing and running
content is available in this source tree.

### Building and Installing

Amanzi can be built with standard C++/C/Fortran compiler families
(GCC, Intel, Clang + gfortran) and a standard MPI library (Mpich,
OpenMPI).  It uses CMake for its build system.

In addition, Amanzi, depends on several Third Party Libraries (TPLs),
including Trilinos and PETSc, and so special instructions for our 
TPL metabuilder "bootstrap.sh" can be found in

[doc/build_instructions/building_bootstrap.rst](doc/building_instructions/building_bootstrap.rst)

Also, you can get help by typing

./bootstrap.sh --help

Note the bootstrap metabuilder will download the currently supported
version of all the TPLs for you, and hence requires about 3.5GB of 
available diskspae (350MB for the downloads, then upacking and building),
and the build of Amanzi requires about another 2.2G of diskspace.

### Running

The instructions for getting started and running Amanzi are included
in the User Guide (link included above), and can be found in

[doc/user_guide/background/getting_started.rst](/doc/user_guide/background/getting_started.rst)

Amanzi comes with a number of tests that are used to generate the User
Guide.  These are in testing/verification and testing/benchmarking
users are encouraged to play around with them.  These tests include
python scripts that run Amanzi, and plot the results.
