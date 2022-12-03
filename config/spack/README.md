# Amanzi-ATS Spack Package 

**DISCLAIMER**: The Amanzi-ATS team is striving to provide complete and concise directions on how to build Amanzi-ATS with Spack, but please be aware there are some issues due to instabilities in Spack itself, and challenges with certain configurations. We are continuosly working to improve the Amanzi-ATS Spack package and the installation instructuctions, particularly with machine specific information. Generally, the Spack build of Amanzi-ATS has been successful on all platforms listed here, although there are exceptions.  Please check this documentation [(README.md)](https://github.com/amanzi/amanzi/blob/master/config/spack/README.md)
on the master branch for the latest updates.

## Quickstart (Part 1): Setup Spack

### Clone Spack 
Download Spack on your machine
```
git clone https://github.com/spack/spack
```

Note that Spack is improving rapidly, and with our large number of TPLs we've found using a specific version/commit on the develop branch to be more robust and consistent.  We encourage you to try this commit for your builds,   
```
cd spack
git reset --hard 45043bcdf5b383af573375300f2ba6b28210e3c6 
```

### Add Spack to your environment

Add Spack in your environment (and consider adding it to your shell initialization (e.g., to .bash_profile)):
```
source ${PATH_TO_SPACK}/spack/share/spack/setup-env.sh
```

### Setup compilers and tools for Spack

Find compilers:
```
spack compiler find
```

Choose the compiler and its version for your builds.  The command
```
spack compilers
```
displays the list of compilers that spack (e.g., gcc@10.2.0). 

Find external libraries:
```
spack external find
```
**NOTE**: letting Spack find external libraries could be dangerous because it may cause Spack to select system libraries during the Amanzi build instead of the libraries required by the package installation (that would be installed as requested by our package when not found on the system). Please be aware of the implications of letting Spack find external libraries for you.


## Quickstart (Part 2): Build Amanzi or ATS with Spack

### Clone Amanzi

To get the source code clone the Amanzi repository:
```
git clone --recursive https://github.com/amanzi/amanzi
```

### Add the Amanzi spack repositories

Change to the _amanzi_ subdirectory, and then add the Amanzi repo to the list of Spack repos:

```
spack repo add amanzi/config/spack
```

This command adds repositories (spack packages) for the following nine packages:

```
alquimia	amanzi		ascemio		ccse		crunchtope	mstk		petsc		pflotran	trilinos
```

Note that packages such as alquimia, petsc, pflotran and trilinos are maintained by importing the published spack packages and then augmenting them with Amanzi-specific requirements such as patches, additional dependencies, or additional compiler flags.

### Install Amanzi with Spack:

With spack configuration complete, all you need to type is
```
spack install amanzi %compiler@version
```

Once this process completes successfully, it's useful to have spack 'load' the amanzi installation
```
spack load amanzi
```
This command will add amanzi to you PATH and setup any other environment variables that are needed.

Now you can run amanzi.  As a quick test, for example, you can run
```
amanzi --print_version
```
to show the branch and global git hash associated with the source you just built.

### Install with a specific version of openmpi

See what versions of openmpi are available on Spack and select one:

```
spack info openmpi
```
Install Amanzi with custom openmpi:

```
spack install amanzi %compiler@version ^openmpi@mpiVersion
```

## Amanzi-ATS Spack Package Variants
The Amanzi-ATS spack package supports the following options, 

```
Currently, we are striving to support the following variants:
    Name [Default]                 When    Allowed values          Description
    ===========================    ====    ====================    ============================================

    build_type [RelWithDebInfo]    --      Debug, Release,         CMake build type
                                           RelWithDebInfo,         
                                           MinSizeRel              
    data_model [epetra]            --      epetra, tpetra          Trilinos data model
    geochemistry [off]             --      on, off                 Enable geochemistry support
    hypre [on]                     --      on, off                 Enable Hypre solver support
    mesh_framework [mstk]          --      mstk, moab              Unstructured mesh framework
    mesh_type [unstructured]       --      unstructured,           Select mesh type: unstructured or structured
                                           structured              (the structured option is currently NOT supported)
    physics [amanzi]               --      amanzi, ats             Physics implementation
    shared [on]                    --      on, off                 Build Amanzi shared
                                                                   (note that this variant ONLY applies to Amanzi and NOT to its dependencies)
    silo [off]                     --      on, off                 Enable Silo reader for binary files
    tests [on]                     --      on, off                 Enable the unit test suite
                                                                   (currently working on the OFF option for tests)
```

For example, to build ATS with geochemistry you could run,
```
spack install amanzi +geochemistry physics=ats
```


## Detailed instructions with system specific directions

This expands on the _quick_ _start_ section of these instructions and provides additional details as well as system specific best practices to build Amanzi with Spack.

If compilers are available on your system, you can load them using: 
``` 
spack compiler find
```
This will update the file: ${HOME}/.spack/_system_/compilers.yaml

The _system_ could be for instance a supercomputer (NERSC), local cluster, or your desktop/laptop.

If you are on a system using modules, first load the compiler module and then the Spack command, for instance:
```
module load gcc/XXX
spack compiler find
```

### Find external libraries

If external libraries are available on your system, and if you
don’t want Spack to build them all, Spack can find them for you by doing:
```
spack external find
```
This will update the file: ${HOME}/.spack/packages.yaml.

External packages (such as openmpi or mpich) can be used to specify dependencies when installing spack packages using `^`, for instance on LANL's Darwin:

```
spack install amanzi ^openmpi@4.1.2
```
will install Amanzi with the explicit dependence on the module `openmpi/4.1.2-gcc_11.2.0` to which it is associated the Spack spec `^openmpi@4.1.2`.
Note that specifying a compiler is not necessary here since it is already taken into account in the module.

The spec is in the packages.yaml file, see below.

### Manage external packages

In ${HOME}/.spack/packages.yaml you can find the external libraries detected by spack, and also add others that might be associated with existing modules
present in your system. For instance, on LANL's Darwin, the module `openmpi/4.1.2-gcc_11.2.0` can be used to create the following spec:

```
openmpi:
  externals:
    - spec: openmpi@4.1.2
      modules:
        - openmpi/4.1.2-gcc_11.2.0
```
If your packages.yaml file was blank, then the line `packages:` has to be added above `openmpi:`. If other specs are already present in the packages.yaml file,
then you can just add more at the bottom and omit the `packages:` line before the library name.

Note that if you ran `spack external find` the above spec may have been added to the packages.yaml file and look like this:

```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        prefix: /projects/opt/centos8/x86_64/openmpi/4.1.2-gcc_11.2.0
```

While the above spec might work in some cases, we have found that it is safer to just use the module syntax.



## Building Amanzi with Spack on Specific Systems

Depending on the different systems you might be using, it is recommended to use specific specs for the Amanzi build.
Below you can find some system specific tested commands for the Amanzi installation via Spack that may help you identify the best approach for your system.

### Ubuntu Linux Compute Server (LANL, T-5)

These systems are fairly generic Ubuntu 18.04 compute servers, and use modules to manage tools (e.g., compilers, mpi, cmake, etc.). 
The following steps have been tested on the `piquillo` servers. 

First, LANL systems are behind a proxy, so the proxy variables need to be set (if you don't already have it in your environment):
```
export http_proxy=proxyout:8080
export https_proxy=$http_proxy
```
Then, load the module `gcc/11.2.0` and then do a `spack compiler find` to be able to use `gcc@11.2.0`.
Next, add the following spec to the packages.yaml file:
```
  openmpi:
    externals:
    - spec: openmpi@4.0.4
      modules:
      - openmpi/4.0.4/gcc-11.2.0
```
To build Amanzi, run from any folder:
```
spack install amanzi <desired_variant> ^openmpi@4.0.4 %gcc@11.2.0
```

### HPC Linux Clusters (LANL, Darwin)
Darwin is a collection of advanced compute nodes administered as an HPC Linux cluster.  It uses modules to manage tools, and slurm as the queuing system.

From either a login node or a compute node, load the following modules: `miniconda3/py39_4.12.0` and  `openmpi/4.1.2-gcc_11.2.0`.
Then, run a `spack compiler find` to add `gcc@11.2.0` to the list of available compilers for Spack. Next, add the following spec to the packages.yaml file:
```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        modules:
          - openmpi/4.1.2-gcc_11.2.0
```
To build Amanzi, run from any folder:
```
spack install amanzi <desired_variant>  ^openmpi@4.1.2
```
Note that tests are currently run as part of the Spack build by default (the `+tests` variant is on by default).

### MacOS (Native tools with Spack) 

Macs are _almost_ standard unix systems, but the _almost_ can make things a little tricky.  Specifically, there are a number of options for augmenting OSX with the missing tools and components, such as brew and macports.  Here we outline an approch that minimizes the role of those additional tools and focuses on using the tool chain that is provided with XCode.  

Note, however, that some Amanzi dependencies need a fortran compiler and Apple does not include one in its standard tools. Here we use spack itself to add this to our setup.  First, look at the ${HOME}/.spack/_system_/compilers.yaml file and see if fortran compilers are present. In a case when they are not present, the apple-clang spec might look something like this:

```
compilers:
- compiler:
    spec: apple-clang@12.0.0
    paths:
      cc: /usr/bin/clang
      cxx: /usr/bin/clang++
      f77: null
      fc: null
    flags: {}
    operating_system: catalina
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

If that is the case, one option is to install a different compiler such as gcc and then create a [mixed toolchain](https://spack.readthedocs.io/en/latest/getting_started.html):

```
spack install gcc@11.2.0
```

Once it is installed, do:

```
spack load gcc
```

And then:

```
spack compiler find
```

This will add the newly installed gcc compiler to ${HOME}/.spack/_system_/compilers.yaml. The installed spec might look like this:

```
- compiler:
    spec: gcc@11.2.0
    paths:
      cc: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/gcc
      cxx: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/g++
      f77: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/gfortran
      fc: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/gfortran
    flags: {}
    operating_system: catalina
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

Then, to create the mixed toolchain, just copy the fortran lines over to the apple-clang compiler as follows:

```
compilers:
- compiler:
    spec: apple-clang@12.0.0
    paths:
      cc: /usr/bin/clang
      cxx: /usr/bin/clang++
      f77: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/gfortran
      fc: /Users/$USER/Repos/spack/opt/spack/darwin-catalina-skylake/apple-clang-12.0.0/gcc-11.2.0-pgkoliksg6eooixvz37lfydfvfvnr67y/bin/gfortran
    flags: {}
    operating_system: catalina
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

Amanzi can then be installed with:

```
spack install amanzi %apple-clang@12.0.0
```

### Cori (NERSC)

From a login node, run a `spack compiler find` to have `gcc@8.3.0` added to the list of available compilers for Spack. Then, add the following specs to the packages.yaml file:
```
packages:
  mpich:
    externals:
     - spec: mpich@7.7.19
       modules:
         - cray-mpich/7.7.19
  cray-libsci:
    externals:
    - spec: cray-libsci@20.09.1
      modules:
      - cray-libsci/20.09.1
```
To build Amanzi, run from any folder:

```
spack install amanzi <desired_variant> ^mpich@7.7.19 ^cray-libsci@20.09.1 %gcc@8.3.0
```

The specific dependence on the cray scientific library is to make sure that `blas` and `lapack` are taken from there.
Note that with versions of gcc higher than 8.3.0 there could be an argument mismatch error when building some libraries (such as [pflotran](https://github.com/spack/spack/issues/30498)). In that case, users should specify the following flag in the compilers.yaml file for the desired gcc version: `fflags:-fallow-argument-mismatch`. 

