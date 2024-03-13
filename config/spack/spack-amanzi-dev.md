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

Ideally the simplest command just works.  In other words ideally, 
```
spack install amanzi
```
would just work.  In our other notes we have some stuff on selecting variants, mpi, etc. 

But how do we see what Spack was really doing, can we keep the build around even when it succeeds.  Yes, just use the *--keep\_stage* option, for example, 
```
spack install --keep-stage amanzi +geochemistry -tests ^openmpi@4.0.4 %gcc@10.2.0
then look in /tmp/<user_name>/spack-stage/spack-stage-amanzi-<version>-<spack\_hash>.
```

To build in the current source directory (i.e., top level of the repository clone) as a developer or for other debugging reasons use

```
spack dev-build amanzi <variant_options>
```
This command will create the build directory within the current working directory, and use the current respository as the source.  Log files will be stored in the current working directory. 

## Evaluating the Amanzi build

  * ldd, check amanzi libraries are linked correctly
  * ctest, run all the tests


## Building ATS with Spack




