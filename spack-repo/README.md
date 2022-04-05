Amanzi/ATS Spack package 

# Spack installation 
Download spack on your machine:
```
git clone git@github.com:spack/spack 
```

Add spack in your environement (and maybe your bash_profile script.) 
```
source ${PATH_TO_SPACK}/spack/share/spack/setup-env.sh
```

If you want support for the modules installed in spack, install lmod: 
```
spack install lmod
```

# Find compiler 

If compilers are available on your system you can load them using: 
``` 
spack compiler find
```
This will update the file: ${HOME}/.spack/linux/compilers.yaml

If you are on a system using modules, first load the compiler module and then the spack command. 
```
module load gcc/XXX
spack compiler find
```

Note that you might find your compilers also on ${HOME}/.spack/cray/compilers.yaml (if you are on Cori for instance).

# Find external libraries

If external libraries are available on your system, spack can find them for you by doing:
```
spack external find
```
This will update the file: ${HOME}/.spack/packages.yaml
External packages (such as openmpi or mpich) can be used to specify dependencies when installing spack packages using ^, for instance on Darwin

```
spack install amanzi@spack ^openmpi@4.1.2
```
will install amanzi with the explicit dependence on the module openmpi/4.1.2-gcc_11.2.0 to which it is associated the spack spec ^openmpi@4.1.2
The spec is in the packages.yaml file, see below.

# Manage external packages

In ${HOME}/.spack/packages.yaml you can find the external libraries detected by spack, and also add others that might be associated with existing modules
present in your system. For instance, on Darwin, the module openmpi/4.1.2-gcc_11.2.0 can be used to create the following spec

```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        modules:
          - openmpi/4.1.2/gcc-11.2.0
```
The first line above can be omitted if other specs are already present in the packages.yaml file.
Note that if you ran spack external find the above spec may have been added to the packages.yaml file and look like this

```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        prefix: /projects/opt/centos8/x86_64/openmpi/4.1.2-gcc_11.2.0
```

While the above spec might work in some cases, we have found that it is safer to just use the module syntax.

# Amanzi Spack

This current version of Amanzi's Spack package is not (yet) available on the remote Spack repository. 
You will need to download Amanzi/ATS: 

```
git clone --recursive git@github.com:amanzi/amanzi
cd amanzi
git checkout spack
git pull
```

You will then be able to add the repository in your local spack: 

```
spack repo add amanzi/spack-repo
```

The command above will add repositories for the following four packages:

```
amanzi  ascemio  crunchtope  mstk  ccse
```

# Spack variants

Variants:
    Name [Default]                 When    Allowed values          Description
    ===========================    ====    ====================    ============================================

    build_type [RelWithDebInfo]    --      Debug, Release,         CMake build type
                                           RelWithDebInfo,         
                                           MinSizeRel              
    data_model [epetra]            --      epetra, tpetra          Trilinos data model
    geochemistry [off]             --      on, off                 Enable geochemistry support
    hypre [on]                     --      on, off                 Enable Hypre solver support
    ipo [off]                      --      on, off                 CMake interprocedural optimization
    mesh_framework [mstk]          --      mstk, moab              Unstructure mesh framework
    mesh_type [unstructured]       --      unstructured,           Select mesh type: unstructured or structured
                                           structured              
    physics [amanzi]               --      amanzi, ats             Physics implementation
    shared [on]                    --      on, off                 Build shared libraries and TPLs
    silo [off]                     --      on, off                 Enable Silo reader for binary files
    tests [on]                     --      on, off                 Enable the unit test suite

# Building and testing Amanzi with spack

Depending on the different systems you might be using, it is recommended to use specific specs for the Amanzi installation.
Below you can find some system specific tested commands for the amanzi installation via spack

# Darwin
Go on a compute node and add the folling spec to the packages.yaml file:
```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        modules:
          - openmpi/4.1.2/gcc-11.2.0
```
Then, cd into amanzi and build and test amanzi by doing:
```
spack dev-build --test all amanzi@spack <desired_variant>  ^openmpi@4.1.2
```
The dev-build command has to be run from within the amanzi folder and it allows to add the --test all command in order
to test amanzi. Spack does have a spack test run command but at the moment we do not have an amanzi test suite for spack.

# Cori
Compute nodes cannot be allocated for enough time to compute the spack installation of amanzi, so it is better to proceed from a login node.
Add the following specs to the packages.yaml file:
```
packages:
  mpich:
    externals:
     - spec: mpich@7.7.10
       modules:
         - cray-mpich/7.7.10
  cray-libsci:
    externals:
    - spec: cray-libsci@20.09.1
      modules:
      - cray-libsci/20.09.1
```
Then, cd into amanzi and build and test amanzi by doing:

```
spack dev-build --test all amanzi@spack ^mpich@7.7.19 ^cray-libsci@20.09.1 %gcc@8.3.0
```

The specific dependence on the cray scientific library is to make sure that blas and lapack are taken from there.

# Piquillo

# Laptop

# Notes

Currently spack does not propagate variants to dependencies, hence the static variant (shared = off) is currently a work in progress.

The ccse spack package is currently undergoing testing so the structured variant (mesh_type=structured) is currently unavailable.

