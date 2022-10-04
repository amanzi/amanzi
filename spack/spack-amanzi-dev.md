# Spack Package for Amanzi

The goal of this file is to document the details and more subtle parts
of the Amanzi spack package.  This document will augment github issues
with a more complete and coherent guide to details of the TPLs build,
the Amanzi and ATS builds and testing.  Simplified instructions will
go in the README.md

## Building the Tools (compiler, mpi, cmake)

This is part of our bootstrap.sh script, but is really covered in the spack-tips.md discussion. and so I think we need to provide some guideance here.  Although, just pointing to the section in spack-tips.md is probably enough.


## Building the Third Party Libraries with Spack

Here we can dig in to the build (and dependencies) of each TPL,
describe how each TPL can be built directly with spack install, and
how they can be built as a collection.

### List of TPLs, Versions and Dependencies

The file /config/SuperBuild/TPLVersions.cmake provides the list of TPLs and their versions (at least the versions we have tested with).

Dependencies can be found with spack spec, or looking at the CMake files for each TPL within our current build system.



## Building Amanzi with Spack




## Evaluating the Amanzi build

  * ldd, check amanzi libraries are linked correctly
  * ctest, run all the tests


## Building ATS with Spack




