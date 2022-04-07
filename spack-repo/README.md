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

Note that you might find your compilers also in ${HOME}/.spack/cray/compilers.yaml (if you are on Cori for instance).

# Find external libraries

If external libraries are available on your system, spack can find them for you by doing:
```
spack external find
```
This will update the file: ${HOME}/.spack/packages.yaml

External packages (such as openmpi or mpich) can be used to specify dependencies when installing spack packages using `^`, for instance on Darwin:

```
spack install amanzi@spack ^openmpi@4.1.2
```
will install amanzi with the explicit dependence on the module `openmpi/4.1.2-gcc_11.2.0` to which it is associated the spack spec `^openmpi@4.1.2`.

The spec is in the packages.yaml file, see below.

# Manage external packages

In ${HOME}/.spack/packages.yaml you can find the external libraries detected by spack, and also add others that might be associated with existing modules
present in your system. For instance, on Darwin, the module `openmpi/4.1.2-gcc_11.2.0` can be used to create the following spec

```
openmpi:
  externals:
    - spec: openmpi@4.1.2
      modules:
        - openmpi/4.1.2/gcc-11.2.0
```
If your packages.yaml file was blank, then the line `packages:` has to be added above `openmpi:`. If other specs are already present in the packages.yaml file,
then you can just add more at the bottom and omit the `packages:` line before the library name.

Note that if you ran `spack external find` the above spec may have been added to the packages.yaml file and look like this

```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        prefix: /projects/opt/centos8/x86_64/openmpi/4.1.2-gcc_11.2.0
```

While the above spec might work in some cases, we have found that it is safer to just use the module syntax.

# Amanzi Spack

This current version of amanzi's spack package is not (yet) available on the remote spack repository. 
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
```
Variants:
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
                                           structured              
    physics [amanzi]               --      amanzi, ats             Physics implementation
    shared [on]                    --      on, off                 Build shared libraries and TPLs
    silo [off]                     --      on, off                 Enable Silo reader for binary files
    tests [on]                     --      on, off                 Enable the unit test suite
```

# Building and testing Amanzi with spack

Depending on the different systems you might be using, it is recommended to use specific specs for the amanzi build.
Below you can find some system specific tested commands for the amanzi installation via spack

# Darwin
Go on a compute node and add the following spec to the packages.yaml file:
```
packages:
  openmpi:
    externals:
      - spec: openmpi@4.1.2
        modules:
          - openmpi/4.1.2/gcc-11.2.0
```
To build amanzi, run from any folder:
```
spack install amanzi@spack <desired_variant>  ^openmpi@4.1.2
```
Currently amanzi does not have a spack testing suite.
To run tests for amanzi after installation, cd into amanzi and replace `install` in the above command string with `dev-build --test all`. Note that the testing cannot be done from any folder but it as to be done from within the amanzi folder.

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
To build amanzi, run from any folder:

```
spack install amanzi@spack <desired_variant> ^mpich@7.7.19 ^cray-libsci@20.09.1 %gcc@8.3.0
```

The specific dependence on the cray scientific library is to make sure that `blas` and `lapack` are taken from there.
To test amanzi, proceed as explained above using `dev-build --test all`.

# T-5 systems

The following steps have been carried out on `piquillo`.

Add the following spec to the packages.yaml file:
```
  openmpi:
    externals:
    - spec: openmpi@4.0.4
      modules:
      - openmpi/4.0.4/gcc-11.2.0
```
To build amanzi, run from any folder:
```
spack install amanzi@spack <desired_variant> ^openmpi@4.0.4 %gcc@11.2.0
```
Note that you might have to first load the module `gcc/11.2.0` and then do a `spack compiler find` to be able to use `%gcc@11.2.0`.
To test amanzi, proceed as explained above using `dev-build --test all`.

# Laptop

# Notes

Currently spack does not propagate variants to dependencies, hence the static variant (`-shared`) is currently a work in progress.

The ccse spack package is currently undergoing testing so the structured variant (`mesh_type=structured`) is unavailable at the moment.

