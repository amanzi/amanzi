==============================================
 Required Environment and Software for Amanzi
==============================================




System Software Environment
===========================

The libraries and tools listed in this section are assumed to be part
of the target platform. Typically, a user would not build these
libraries from source.  The Amanzi team recommends that users use
pre-built binaries or platform-specific pacakge management tools to
install these libraries.

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

Amanzi builds the source located in
https://software.lanl.gov/ascem/trac/browser every hour on several OS
types.  The build and test suite results for each OS is recorded and
tracked on the ASCEM Trac site at
https://software.lanl.gov/ascem/trac/buildbot.


Required Third-Party Software
=============================

Amanzi uses existing software libraries whenever possible and can not
built unless all the required libraries are present. In some cases
(netCDF, Exodus, HDF5 and Trilinos), the libraries must be configured
in a particular way to work correctly with Amanzi and other TPLs. 

The Amanzi source directory (TODO: link) contains scripts that build
these librarues correctly. Users are strongly encouraged to use the
appropiate tool for their platform. Advanced users creating their own
build process should refer to this script for supported configurations
and TPL build order.


Required Libraries
------------------

Boost:
        :Versions: 1.46.0
        :Description: C++ Library which extends the capabilibies of the standard library.
        :Role: Specific Boost libraries are used in various parts of the Amanzi code base.

        - ``filesystem``: Used to create directories, test directory existence and other
          file system activities.
        - ``graph``: Used in the mesh audit tool.
        - ``mpl``: Meta-programming library used in the mesh factory.  

        :Dependencies: A modern C++ compiler. See the Boost website for specifics.
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
        
        - ``Epetra``: a parallel-aware array libarary
        - ``STKmesh``: a mesh database libary (optional)
        - ``NOX``: Nonlinear Object-Oriented Solutions package  
        - TODO: Complete?

        :Dependencies: ExodusII, (if STKmesh used) CMake, MPI
        :Information: http://trilinos.sandia.gov/

CCSE:
        :Version: 0.1.2
        :Description: Base library for structured-mesh objects
        :Role: Used by Amanzi to implement structured-grid adaptive integrator
        :Dependencies: MPI, OpenMPI (if enabled)
        :Information: https://ccse.lbl.gov/Software/ccse_core.html


Required Software Tools
-----------------------

CMake:
        :Versions: 2.8.3 required.
        :Description: Cross-platform software build system
        :Role: Forms the basis of the Amanzi build and testing tools
        :Dependencies: A suitable build backand. GNU Make is standard.
        :Information: http://www.cmake.org/


Optional Third-Party Software
=============================

These tools and libraries are not essential to create a working Amanzi
installation, but will enable additional Amanzi features, or provide
useful when using Amanzi.


Optional Libraries
------------------

Note that, while each of the mesh database libraries is optional:
STKMesh (above, in Trilinos) MOAB and MSTK, *at least one* of these is
*required* for Amanzi to function.

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
        :Versions: 1.83rc3
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



Optional Software Tools
-----------------------

Mercurial:
        :Versions: TODO: Versions
        :Description: A dirtributed version control system
        :Role: Used by Amanzi to record and track changes to the software, and coordinate developer contributions. Required in order to obtain development versions of the Amanzi source.
        :Dependencies: Python 2.6
        :Information: http://mercurial.selenic.com/

SWIG:  Wait, is this a tool or a library?
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
        :Dependencies: Stand-alone binaries available.
        :Information: 


