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
  
  
## Spack Environments 

Spack environments allow you to create specific containers of spack modules for your project. 
See the spack complete guide on the [environments](https://spack.readthedocs.io/en/latest/environments.html)

Create the environment: 
```shell
spack env create amanzi
```

Load the environment:
```shell
spack env activate amanzi 
```

Unload/quit an enviroment: 
```shell 
spack env disactivate 
# Shortcut (same as previous): 
despacktivate 
```

From inside the environment you can now install packages as you would do using the `spack install` command. 
They will be added to the current environment and loaded directly when the environment is activated. 

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

### Let Spack build and manage it all

Here we'd be looking at spack building gcc, cmake, openmpi, and openblas

ParFlow provides this guideance, it's a good start for us to think about.
Building GCC will take considerable amount of time.

```
spack install gcc@10.3.0
```

Add the compiler to the set of compilers spack can use:

```
spack load gcc@10.3.0
spack compiler find
spack compilers
```

### Configure Spack to use locally installed software

Here we'd be looking at first make sure that cmake, gcc, openmpi are in your path and then running the compiler and externals setup

```
  spack compiler find
```

```
  spack externals find
```

Finally, we'd need to pick and configure blas/lapack within packages/yaml file.  I wonder if we need an entire section on blas/lapack, particularly if we end up supporting mkl.

### Configure Spack to use modules to manage locally installed software

Here we would still use 
```
  spack compiler find
```

```
  spack externals find
```

But then we would update by hand the entries associated with the modules we want spack to be able to use.   Not sure how to make this easiest on the user.


## Working with Spack

### Find available and loaded packages

Here we look at the commands used most often. The `spack find` command can be used with different options. For instance:

```
spack find -ld
```

shows installed packages, their hashes (-l option) and dependency info (-d option). With Spack each package has its own unique hash. Note that there are some packages shown in the `spack find -d` output that may have not been installed explicitly (i.e. with a `spack install` command). These are dependencies that were installed implicitly. A few packages installed implicitly are not shown as dependencies in the `spack find -d` output. These are build dependencies. The following command:

```
spack find -px
```
shows  packages that were installed explicitly rather than pulled in as a dependency (-x option) and their path (-p option). The following command:

```
spack find ^mpich
```
returns every installed package that depends on MPICH. To show the full spec of a package (without dependencies), use:

```
spack find -v

```

To see what has been loaded using `spack load` type:

```
spack find --loaded
```

So to see the list of available builds that are installed
```
spack find -v -l <package_name>
```
will show both the variant info as well as the hash you can use to reference it. 

### Loading a package for easier use

Much like loading a module, loading a package in spack puts the binaries in your PATH and takes care of dependencies.  So for most packages you would load it with the *spack load* command,
```
spack load <package_name> <variants>
```
For example, to load the ats you could provide the variant information as you did for the install,
```
spack load amanzi@master physics=ats
```
Or the hash you found before with the *spack find -v*, 
```
spack load /<sha1_hash>
```
Either way you can see if you packages is loaded as you expected, 
```
spack find --loaded -v amanzi
```
which on the T-5 compute server generates the following output,
```
piquillo4(296)$ spack find --loaded -v amanzi
==> 1 loaded package
-- linux-ubuntu18.04-skylake_avx512 / gcc@10.2.0 ----------------
amanzi@master~geochemistry+hypre~ipo+shared~silo+tests build_type=RelWithDebInfo data_model=epetra mesh_framework=mstk mesh_type=unstructured physics=ats
```

### The spec shows the dependencies

To see the full spec of a package with dependencies use:

```
spack spec <package>
```

### Compilers

To see the list of available compilers, use:

```
spack compilers list
```

To add a compiler to the list of available compilers (without knowing its path) one can use:

```
spack compiler add $(spack location -i gcc@8.4.0)
```

### Repositories for spack packages

To see the list of available repos, use:
```
spack repo list
```

### zlib: A simple example

The following is an example of a spack installation of zlib

```
spack install zlib@1.2.8 %gcc@6.5.0
```

The above command installs `zlib`, version 1.2.8 with `gcc` version 6.5.0. Note that one can install only the dependencies of a package by doing

```
spack install --only dependencies <package>
```

It is also possible to install a certain package by specifying an explicit dependence using a hash. This is helpful in case there are different versions of a package and one wants to make sure the correct one is given as dependency. The syntax is the following:

```
spack install zlib ^/hash
```

Recall that the hash of a certain package can be found by doing `spack find -l`. Variants of packages can be turned on with `+` or off with either `-` or `~` when specifying the install command. For instance:

``` 
spack install hdf5 ~mpi +fortran
``` 

installs hdf5 without mpi support but with fortran support.

 









