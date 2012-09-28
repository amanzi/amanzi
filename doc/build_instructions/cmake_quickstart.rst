============================================================
Amanzi CMake Quick Start Guide
============================================================

Introduction
++++++++++++

This document is an initial introduction to the CMake software package for users and
developers who are not familiar with CMake. It does not document all of the available
features of CMake, but provides an overview of some of the common build settings and the
general CMake approach to create a build system for Amanzi. Links to other web resources 
are listed at the end for tutorials, FAQs and other resources related to CMake.

What Is CMake?
++++++++++++++
CMake, is a cross-platform (Windows, Unix, MacOSX,developing Android and iOS) build system
generator and from this one sentence description of CMake two axioms that any developer more
familiar with UNIX Makefiles should always remember,

**CMake generates build systems for many different OS types and development environments. It does not
exclusively generate UNIX Makefiles.**

and 

**CMake only controls the generation of the build system, and once that is complete, the development
environment controls the software build.**.

To see how these constraints could impact build development consider the ``LD_LIBRARY_PATH`` environment
variable, a list of paths that any compiler on UNIX platform searches for libraries. This
environment variable is not defined on MacOSX or a Windows, and hence CMake does not have a method
to alter or set this environment variable in a build system file. However if this variable is set
when the developer executes ``make`` on a CMake generated Makefile the compiler will still use
``LD_LIBRARY_PATH`` to search for libraries. CMake requires that the
developer think of the build as process creating objects (library, script and stand alone binary)
and defining rules to create these objects. These rules must transcend UNIX platform conventions.  

Why Use CMake?
++++++++++++++

Some of the advantages of CMake are:

#. The ability to create build systems for different development environments such XCode, Eclipse,
   Visual Studio and UNIX Makefiles.

#. Intuitive scripting language to generate libraries, executables and integrated test code.

#. Extensive modules/functions to find and configure a project in the CMake distribution and many more
   available on-line.

#. Easy to create CMake configuration files for a project that can then be imported into other CMake
   files for other projects. 


How To Generate A Build System
++++++++++++++++++++++++++++++

To generate a build system, the source directory must have at least one ``CMakeLists.txt`` which
contains the CMake commands. From the command line, generation begins with ``cmake [OPTIONS] <source>``
where ``<source>`` is the root source directory of the project containing a ``CMakeLists.txt`` file. CMake
will deleve into subdirectories of this source directory if it finds a ``ADD_SUBDIRECTORY(dir1 dir2 dir3 ..)`` 
command. CMake will continue down through each directory until it no longer detects a subdirectory to add,
and then return to the root directory and continue processing the remaining directories. For example, if a
``CMakeLists.txt`` file had ``ADD_SUBDIRECTORY(dir1 dir2 dir3)`` and the ``CMakeLists.txt`` file in ``dir1``
had ``ADD_SUBDIRECTORY(dirA dirB)``, the processing order would be: ``dir1,dir1/dirA,dir1/dirB,dir2,dir3``. CMake 
provides documentation of all the available commands from the command using the ``--help-html``, ``--help-full``
and ``--help-man`` options. The table below list a few of the most common commands:

+-----------------------------------------------+------------------------------------------------------------------------+
| CMake Command                                 | Description                                                            |
+===============================================+========================================================================+
| PROJECT(project_name [lang1 lang2 ...])       |   Set the project name to **project_name**, will also define           |
|                                               |   ``<project_name>_SOURCE_DIR`` and ``<project_name>_BINARY_DIR``      |
|                                               |   which are the full path names of ther source and build directories.  |
+-----------------------------------------------+------------------------------------------------------------------------+
| SET(<var_name> <var_value> ... )              |   Set CMake variable **var_name** to **var_value** additional options  |
|                                               |   control the scope of this variable and it's definition as a          |
|                                               |   cache variable. See the CMake documentation for details.             |
+-----------------------------------------------+------------------------------------------------------------------------+
| ADD_LIBRARY(<name> [OPTIONS] file1  ... )     |  Add a library to the project. Files ``file1 file2 ...`` are           |
|                                               |  project source files or headers required by the library. Creates a    | 
|                                               |  target **name** that builds the library. The name of the              |
|                                               |  library follows the platform conventions. This function is also       |
|                                               |  use to create imported library targets.                               |
+-----------------------------------------------+------------------------------------------------------------------------+
| ADD_EXECUTABLE(<name> [OPTIONS] source1 ... ) |  Add executable to the project. Creates target **name** that           |
|                                               |  builds a stand-alone binary. Final executable name follows platform   |
|                                               |  conventions.                                                          |
+-----------------------------------------------+------------------------------------------------------------------------+
| INCLUDE_DIRECTORIES(dir1 dir2 dir ... )       |  Add directories ``dir1``, ``dir2``, etc. as include directories       |
|                                               |  to the build, i.e. ``-Idir1``, ``-Idir2``, etc.                       |
+-----------------------------------------------+------------------------------------------------------------------------+
| INSTALL(...)                                  |  Function that defines the install rules. Several options that         |
|                                               |  control binary, library and file installation. See the CMake          |
|                                               |  documentation.                                                        |
+-----------------------------------------------+------------------------------------------------------------------------+
| ENABLE_TESTING()                              |  Activates the ADD_TEST() to build the ``test`` target. If this        |
|                                               |  call is not made ADD_TEST does nothing.                               |
+-----------------------------------------------+------------------------------------------------------------------------+
| ADD_TEST(testname Executable [arg1 arg2 ...]) |  If ENABLE_TESTING command has been called, add test labeled           |
|                                               |  ``testname`` to the project.                                          |
+-----------------------------------------------+------------------------------------------------------------------------+
| IF(...)/ELSE()/ELSEIF()/ENDIF()               |  Conditional command control                                           |
+-----------------------------------------------+------------------------------------------------------------------------+
| FOREACH()/WHILE()                             |  Loop control commands                                                 |
+-----------------------------------------------+------------------------------------------------------------------------+
| FIND_FILE                                     |  Set of commands that find the full path to a file, a library,         | 
| FIND_LIBRARY                                  |  the directory containing a particular file and an executable.         |
| FIND_PATH                                     |  The search rules CMake follows is extensive. Consult the CMake        |
| FIND_PROGRAM                                  |  documentation for a full explanation.                                 |
+-----------------------------------------------+------------------------------------------------------------------------+
| FIND_PACKAGE(PackName [OPTIONS])              |  Search for a CMake configuration for software package PackName        |
|                                               |  or execute the CMake module file Find**PackName**.cmake. This         |
|                                               |  is the primary method to include header files and libraries that      |
|                                               |  are external to the project.                                          |
+-----------------------------------------------+------------------------------------------------------------------------+

The directory where the original ``cmake`` command is envoked is referred to as the build or binary directory.
CMake *strongly* encourages building outside the source directory. When CMake finishes processing the source
directory there are several new directories and files present. Two that any user should be aware of is 
the ``CMakeCache.txt`` file and the ``CMakeFiles`` directory.  



A Simple Library Example
++++++++++++++++++++++++

In this section, will explain a simple example creating a project that builds and installs a library.
In this project, ``Foo`` builds and installs a C++ library.

First step is to define a project name, the additional ``C++`` argument defines the project's language.
This will envoke CMake language checks.

::

# CMake Project name is Foo
# This call will define Foo_SOURCE_DIR and Foo_BINARY_DIR
PROJECT(Foo C++)

Next, this library needs source files and header files to build. To avoid updating the ``CMakeLists.txt`` 








Additional CMake Resources
++++++++++++++++++++++++++

CMake User Guide
----------------
Any CMake binary will generate a full user's guide from the command line by calling the
binary with one of the following options,

* ``cmake --help-full [file]`` for standard text output,

* ``cmake --help-html [file]`` for HTML output, or  

* ``cmake --help-man [file]`` for man page output.  

If a file ``[file]`` is specified, then the output is written to that file. This output
lists and documents all the CMake commands, available ``Find*`` modules, and target properties.

CMake Tutorials
---------------

The Cmake projects has web pages on how to get started with CMake. These are

* the basic tutorial, http://www.cmake.org/cmake/help/cmake_tutorial.html, and

* and a series of webinars http://www.cmake.org/cmake/resources/webinars.html.

Looking at what other projects have done with CMake is also an excellent way to learn CMake. 
We found the following software projects to be well organized and had clear CMake usage guidelines
that developers will recognize in the Amanzi project. They are:

* ParaView, an open-source data analysis and visualization tool, http://www.cmake.org/cmake/resources/webinars.html
  ParView builds on both UNIX and Window systems using CMake. It also builds the entire third party
  software stack using the ExternalProject CMake module. Amanzi's SuperBuild project design was based on ParaView's
  SuperBuild design.

* HDF5, http://www.hdf5.org, a well-known data file format library, is converting to CMake after using GNU
  autoconf. Good example of a project transitioning to CMake from GNU and producing CMake configuration files
  for other CMake projects.

* OpenSceneGraph, a 3D graphics toolkit that build on several platforms and operating systems
  http://www.openscenegraph.org/projects/osg

For more projects that use CMake, see the CMake Wiki page http://en.wikipedia.org/wiki/CMake.
