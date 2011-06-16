=====================================
 Amanzi Extern Software Requirements
=====================================

System Software Environment
===========================

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

Required Third-Party Software
=============================

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

Required Libraries for Building Amanzi
--------------------------------------

CMake:
        :Versions: 2.8.3 required.
        :Description: Cross-platform software build system
        :Role: Forms the basis of the Amanzi build and testing tools
        :Dependencies: Stand-alone
        :Information: http://www.cmake.org/


Boost:
        :Versions: 1.46.0
        :Description: C++ Library which extends the capabilibies of the standard library.
        :Role: Specific Boost libraries are used in various parts of the Amanzi code base.

        - filesystem. Used to create directories, test directory existence and other
          file system activities.
        - graph. Used in the mesh audit tool.
        - mpl. Meta-programming library used in the mesh factory.  

        :Dependencies: A good C++ compiler and reasonably modern platform.
        :Information: http://www.boost.org/


zlib:
        :Versions: 1.2.5
        :Description: Compression library
        :Role: Used by HDF5 for input and output.
        :Dependencies: 
        :Information: http://zlib.net/


cURL:
        :Versions: 7.21.2
        :Description: Download tool.
        :Role: Used by netCDF  
        :Dependencies:
        :Information: http://curl.haxx.se/


HDF5:
        :Versions: 1.8.7
        :Description: File format libary
        :Role: Used by Amanzi for input and output of problem data.
        :Dependencies: zlib
        :Information: http://www.hdfgroup.org/HDF5/


netCDF:
        :Versions: 4.1.1
        :Description: File format libary
        :Role: Used by Amanzi for input and output of mesh data.
        :Dependencies:  cURL, HDF5
        :Information: http://www.unidata.ucar.edu/software/netcdf/


ExodusII:
        :Versions: 4.98
        :Description: Mesh data base libary
        :Role: Used by Amanzi to describe mesh geometry for import.
        :Dependencies: netCDF
        :Information: http://sourceforge.net/projects/exodusii/


Trilinos:
        :Versions: 10.6.4
        :Description: Library collection of tools for numberic computing in C++
        :Role: Used throughout Amanzi for data structures and algotithms
        
        - Eperta, a parallel-aware array libarary
        - STKmesh, a mesh database libary (optional)
        - NOX. Nonlinear Object-Oriented Solutions package  

        :Dependencies: ExodusII, (if STKmesh used) CMake, MPI
        :Information: http://trilinos.sandia.gov/


Optional Third-Party Libraries used in Amanzi
---------------------------------------------

These libaries are not required to build Amanzi, but will provide it
with additional capabilities.

Note that, while each one of the mesh database libraries is listed as
optional: STKMesh (above, in Trilinos) MOAB and MSTK, at least one of
these is required for Amanzi to function.

UnitTest++:
        :Versions: 1.4
        :Description: C++ Unit test creation framework
        :Role: Used to build Amanzi unit tests
        :Dependencies: 
        :Information: http://sourceforge.net/projects/unittest-cpp/


MOAB:
        :Versions: Revision 4225 from the SVN repository
        :Description: A Meshing database library
        :Role: An optional backend for Amanzi meshing
        :Dependencies: ExodusII
        :Information: 

MSTK:
        :Versions: 1.80
        :Description: A meshing database library
        :Role: An optional backend for Amanzi meshing
        :Dependencies: ExodusII
        :Information:  https://software.lanl.gov/MeshTools/trac

ASCEM-IO:
        :Versions: 1.0
        :Description: A parallel I/O library
        :Role: Basis library used by the output class to create problem output.
        :Dependencies: HDF5, netCDF
        :Information: http://ascem-io.secure-water.org 



Optional Third-Party Software Tools
-----------------------------------

SWIG:
        :Versions:
        :Description:
        :Role: 
        :Dependencies:
        :Information: 

XDMF:
        :Versions:
        :Description:
        :Role: 
        :Dependencies:
        :Information: 

VisIt:
        :Versions: 2.2.2
        :Description: Visualization program
        :Role: Tool to view output data. 
        :Dependencies: Pre-built binaries available (VERY difficult to build)
        :Information: https://wci.llnl.gov/codes/visit/home.html

Doxygen:
        :Versions:
        :Description: A source-code to documentation tool.
        :Role: Used to create the Amanzi code documentation and test descriptions.
        :Dependencies:
        :Information: 


