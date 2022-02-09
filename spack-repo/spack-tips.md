# Spack Tips for Amanzi

To build Amanzi, and for many spack workflows, you only need a few simple spack commands.  This document highlights the commands we find most useful, gives examples of how to use them, and what information they provide.

## Spack Commands

Here's a quick list of the commands we use most often:

  * Setup and Build Environment 
    * spack compiler find
    * spack external find
  * Build and Install
    * spack install <package>
    * spack spec <package>
  * ...
  
  
## Spack Setup

Here we need to develop instructions on how to setup the compilers,
mpi, cmake, and other tools that we want for our build environment.
I think there are three cases we need to discuss here.

  * Let spack build everything (including, compilers, mpi, cmake)
  * Tell spack where local cmake, gcc, mpi are located (without modules)
  * Tell spack where local cmake, gcc, mpi are located with modules
 
A condensed/simple version of these instructions would probably go in
the README.md. But here I'm hoping we can include a bit more detail
and provide instructions on how to verify that the setup is correct.

## Working with Spack

Here we look at the commands used most often, e.g., install, spec,
info, find etc.




