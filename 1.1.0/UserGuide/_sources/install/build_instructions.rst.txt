==============================================
Building and Installing Amanzi Using CMake
==============================================

Building Amanzi
---------------

Amanzi depends on several external software packages or third party libraries
(TPL) to build. Before building Amanzi a user must either build each of these 
TPLs, provide an installation path for each TPL or define the include directories
and libraries needed for each TPL. We recommend that users build Amanzi and all the
TPLs using the ``boostrap.sh`` shell script found in the ``config`` directory. 
The script can be run from the source directory. On most
UNIX flavored systems, this script will build and install Amanzi and all the required
TPLs with little information from the user. The script has prerequisites for
compilers, ``LAPACK``/``BLAS`` and ``MPI``. Please read 
:ref:`TPLs` section for more information and :ref:`bootstrap` for installing using the ``bootstrap.sh``.

We recommended that users install the TPLs in separate directory, since a full TPL build
should be an infrequent task. This is the default behavior of the ``bootstrap.sh`` 
script. See ``bootstrap.sh --help`` for more information on how to
control the installation paths. The TPL build will generate a CMake configuration file
that defines all the required TPL directories and CMake variables to build Amanzi.
Once the user has installed all the required TPLs, the process to build
Amanzi from the command line is:

  1. Generate build files (Makefile, XCode project files, etc.) using CMake ::

       cmake -C <TPL install prefix>/share/cmake/amanzi-tpl-config.cmake \
             -D CMAKE_C_COMPILER=/full/path/mpicc \
	     -D CMAKE_CXX_COMPILER=/full/path/mpicxx \
	     -D CMAKE_Fortran_COMPILER=/full/path/mpif90 \
             [Addtional CMake variable definitions]
             <Amanzi source directory>

  2. Use the build files to build Amanzi, i.e. ``make``, ::
   
       make [-j n]

or the user can run ``bootstrap.sh`` defining a previously built TPL configuration with ::

  bootstrap.sh --tpl-config-file=<TPL install prefix>/share/cmake/amanzi-tpl-config.cmake


Amanzi TPL Configuration Settings
+++++++++++++++++++++++++++++++++

We recommend that users build all the TPLs using the Amanzi SuperBuild project.
This will create a CMake file that will initialize the cache settings for all 
the TPLs. Instructions in this section are designed for advanced users. 

The Amanzi CMake files search for each TPL through the CMake find_package function.
Each TPL has a search path defined by the variable ``<TPL_Name>_DIR`` where
TPL_Name is the name of the package. The Amanzi CMake files search for the
appropriate include directories and libraries for this package. The search
also includes searching for dependent packages. For example, searching for
ASCEM-IO also triggers a search for HDF5 since ASCEM-IO depends on HDF5. 
The exceptions to this variable naming convention are Boost and HDF5. For these
packages, Amanzi uses the CMake FindBoost module found in the CMake installation
and the FindHDF5 module found in the ``tools/cmake/Modules`` directory.

Users can bypass the package search by defining the ``<TPL_Name>_INCLUDE_DIRS``
and ``<TPL_Name>_LIBRARIES`` variables. These variables should contain the full
path name of the include directories and libraries for that TPL and ANY
dependent package the TPL requires to build and link against. For example,
if you want to bypass the ExodusII search logic, then define
``ExodusII_INCLUDE_DIRS`` and ``ExodusII_LIBRARIES`` and these variables must
also contain the NetCDF include directories and libraries. These variables
are CMake list types. Each directory or filename should be separated with 
a semicolon.

Amanzi Build Configuration Settings
+++++++++++++++++++++++++++++++++++

The naming convention for Amanzi build options is ENABLE_<Feature_Name>
and is a boolean type (ON,OFF). The current build options with the default values
are documented here.


ENABLE_Structured:

      - Default: OFF
      - Description: Build structured mesh capability. 
      - Dependencies: CCSE


ENABLE_Unstructured:

      - Default: OFF
      - Description: Build unstructured mesh capability.
      - Dependencies: At leat one of the mesh frame works, STK, MSTK or MOAB.


ENABLE_DBC:

      - Default: ON
      - Description: Enable design by contract build.
      - Dependencies:


ENABLE_Config_Report:

      - Default: ON
      - Description: Print out configuration report to the terminal.
      - Dependencies:

ENABLE_MSTK_Mesh:

      - Default: OFF
      - Description: Build the MSTK mesh frame work.
      - Dependencies: MSTK


ENABLE_MOAB_Mesh:

       - Default: OFF
       - Description: Build the MOAB mesh frame work.
       - Dependencies: MOAB, requires a specific version. See Software Requires for more information.

ENABLE_CCSE_TOOLS:

       - Default: OFF
       - Description: Build structured AMR post processing tools 
       - Dependencies: python, f2py, gfortran compatible with mpif90

ENABLE_UnitTest:

       - Default: ON
       - Description: Build the unit test test suite.
       - Dependencies: UnitTest++


ENABLE_OpenMP:

       - Default: OFF
       - Description: Build Amanzi executables with OpenMP support.
       - Dependencies: OpenMP


Installing Amanzi
-----------------

CMake will generate an ``install`` target in build files. For Makefiles,
``make install`` will install Amanzi under the directory defined by 
``CMAKE_INSTALL_PREFIX``. The default install location is ``/usr/local``.

Once installed, other CMake software projects can build and link against Amanzi
as a library. 
See https://software.lanl.gov/ascem/trac/wiki/Amanzi/BuildSystemIntegration for a simple example.

