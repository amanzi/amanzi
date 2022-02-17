Amanzi/ATS Spack package 

# Spack installation 
Download spack on your machine:
```
git clone git@github.com:spack/spack 
```

Add spack in your environement (and maybe you bash_profile script.) 
```
source ${PATH_TO_SPACK}/spack/share/spack/setup-env.sh
```

If you want support for the modules installed in spack, install lmod: 
```
spack install lmod
```

## Find compiler 

If compilers are available on your system you can load them using: 
``` 
spack compiler find
```
This will update the file: ${HOME}/.spack/linux/compilers.yaml

If you are on a system using module, first load the compiler module and then the spack command. 
```
module load gcc/XXX
spack compiler find
```

## Add modules

Some modules might be already present on your system. 
You can add them to spack by editing the file: ${HOME}/.spack/linux/packages.yaml

You can then edit this file following this example for openpmi: 

```
packages:
  openmpi: 
    externals:
      - spec: openmpi@3.1.4 
        modules: 
          - openmpi/3.1.4-gcc_9.2.0 
      - spec: openmpi@1.2.3
        modules: 
          - openmpi/1.2.3-gcc_2.3.4
```

# Amanzi Spack

This current version of Amanzi's Spack package is not (yet) available on the remove Spack repository. 
You will need to download Amanzi/ATS: 

```
git clone --recursive git@github.com:amanzi/amanzi
```

You will then be able to add the repository in your local spack: 

```
spack repo add amanzi/spack-repo
```

The basic install option for amanzi will be available as soon as the remote spack will contain the Amanzi package. 
For now we will proceed using a dev-build installation, installing the Amanzi/ATS in your current directory: 

```
spack dev-build amanzi@1.0.0 +alquimia +crunchtope +ATSPhysics +AmanziPhysics 
```

If you want to provide a specific library, for example mpi, you can provide them at the end of the install line: 

```
spack dev-build amanzi@1.0.0 +alquimia +crunchtope +ATSPhysics +AmanziPhysics ^openmpi@3.1.4 
```

# Spack variants: 

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


