= Amanzi and ATS Installation Guide =

First, a word of warning -- please be patient.  Amanzi and therefore ATS depend upon a lot of third party libraries.  This allows us to use existing, mature code to make ATS a much better software tool.  It also means installing the code and its dependencies can be quite painful and time/labor intensive.  Installation time is a very bimodal distribution -- if it "just works" this process will take 25 minutes.  If it doesn't "just work" it can take much longer.

The basics of the process are these.

0. Install pre-reqs.
1. Download amanzi-ats.
2. Use SuperBuild to install TPLs.
3. Install and test Amanzi.
4. Download and install ATS.
5. Download test problems, test ATS.

== Linux-based systems ==

All instructions assume you use bash.  Alter as needed for other shells.

0. Ensure you have cmake and an MPI installation.
  a. cmake >= 2.8.6
    ```
    which cmake  # this will test if you have any cmake installed
    ``` 
    If not, go to: http://www.cmake.org/cmake/resources/software.html
    To install cmake from command line see https://cmake.org/install/
    
    Binary installations are likely fine?
    
  b. MPI
    ```
    which mpicc # see if you have any MPI installed.
    locate mpi.h # see if you have header files installed (is this enough?)
    ```
    On ubuntu, these are available with:
    ```
    sudo apt-get install libopenmpi-dev openmpi-bin  # I believe these are the sufficient conditions?
    ```
    On Mac, How to install the openmpi:
    https://wiki.helsinki.fi/display/HUGG/Open+MPI+install+on+Mac+OS+X
   
    
  c. Lapack (linear algebra package)
    On ubuntu, this is available with:
    ```
    sudo apt-get install libblas3gf libblas-doc libblas-dev liblapack3gf liblapack-doc liblapack-dev # are there all needed?
    ```
    On Mac, Instructions on how to install lapack  
    https://pheiter.wordpress.com/2012/09/04/howto-installing-lapack-and-blas-on-mac-os/
    Follow instruction to install the BLAS. LAPACK requies to extra steps that is not in the instrcutions. 
    Before compiling lapack package go to BLAS/SRC, LAPACKE and CBLAS folders run make for each of them. This will producce three .a files required make.inc. Note, the make.inc in the lapack has to be modified accordingly. To successfully pass all test copy librefblas.a to all TESTING folder. 

1. Set up your directory structure.  Life gets much easier if you get this right first.
  a. I use the following:
    ```
    mkdir /my/path/to/ats
    export ATS_BASE=/my/path/to/ats
    export ATS_SRC_DIR=${ATS_BASE}/repos/ats
    export AMANZI_SRC_DIR=${ATS_BASE}/repos/amanzi
    export ATS_BUILD_DIR=${ATS_BASE}/ats-build
    export AMANZI_BUILD_DIR=${ATS_BASE}/amanzi-build    
    export AMANZI_TPLS_BUILD_DIR=${ATS_BASE}/amanzi-tpls-build    
    export ATS_DIR=${ATS_BASE}/ats-install
    export AMANZI_DIR=${ATS_BASE}/amanzi-install    
    export AMANZI_TPLS_DIR=${ATS_BASE}/amanzi-tpls-install    
    ```    
    You may want to put that line in your ~/.bashrc or similar files (~/.bash_profile on Mac OS X).
  b. Download and unzip the amanzi source tar at https://github.com/amanzi/amanzi/, moving the resulting downloaded directory to $AMANZI_SRC_DIR.  Note that most people should use the current release, from https://github.com/amanzi/amanzi/releases called ``ats-amanzi-X.YY.pZ``, **not** to be confused with either the current dev version or the current amanzi-only release.  Numbering of Amanzi and ATS releases are close, but not always identical.

    Currently this is ``ats-amanzi-0.86.p1``

  c. Download and unzip the ats source tar at https://github.com/amanzi/ats-dev/ .  Again, most people want the most current ATS release at https://github.com/amanzi/ats/releases which will be named ats-X.YY.pW  Note specifically that the patch number (Z in Amanzi, W in ATS) need not match, but the major and minor version numbers X and YY **must match**.

    Currently this is ``ats-0.86.p0``

2. Build the Amanzi TPLs.
  a. Set up a directory for configuration scripts and a TPL installation directory.
    ```
    cd ${ATS_BASE}
    mkdir ats-config-files  # this will hold scripts that call cmake to build things
    mkdir amanzi-tpls # this will be where TPLs get built and installed
    ```
  b. Write a configure script for SuperBuild / Amanzi TPLs.  Edit a file (I call it configure-amanzi-tpls-debug.sh).  Add something like the following text:
    ```
    #!/bin/bash
    # configure-amanzi-tpls.sh
    # script for running CMake to build TPLs for Amanzi.

    rm -rf ${AMANZI_TPLS_BUILD_DIR}
    mkdir -p ${AMANZI_TPLS_BUILD_DIR}
    cd ${AMANZI_TPLS_BUILD_DIR}

    CC=/my/path/to/mpicc  # <--- edit me
    CXX=/my/path/to/mpicxx  # <--- edit me
    FC=/my/path/to/mpif90 # <--- edit me

    cmake \
      -D CMAKE_C_COMPILER=${CC} \
      -D CMAKE_CXX_COMPILER=${CXX} \
      -D CMAKE_Fortran_COMPILER=${FC} \
      -D ENABLE_HYPRE:BOOL=ON \
      -D ENABLE_Structured:BOOL=OFF \
      -D ENABLE_STK_Mesh:BOOL=OFF \
      -D TPL_INSTALL_PREFIX=${AMANZI_TPLS_DIR} \
    ${AMANZI_SRC_DIR}/config/SuperBuild

    make -j8  # <--- edit me to be approximately your number of cores + a few for faster build time.
    make install
    ```

  c. Run the script.
  d. Go get a cup of coffee.  If this works, it will take ~15 minutes, depending upon your machine.  When you come back, it should tell you if it passed or not.


3. Build Amanzi.
  a. Edit another script to build Amanzi.  (I call this one configure-amanzi-debug.sh).  Add something like:
    ```
    #!/bin/bash
    # configure-amanzi.sh
    # script for running CMake to build Amanzi.

    rm -rf ${AMANZI_BUILD_DIR}
    mkdir -p ${AMANZI_BUILD_DIR}
    cd ${AMANZI_BUILD_DIR}

    CC=/my/path/to/mpicc  # <--- edit me
    CXX=/my/path/to/mpicxx  # <--- edit me
    FC=/my/path/to/mpif90 # <--- edit me

    cmake \
       -C  ${AMANZI_TPLS_DIR}/share/cmake/amanzi-tpl-config.cmake \
       -D CMAKE_INSTALL_PREFIX=${AMANZI_DIR} \
       -D ENABLE_Physics:BOOL=OFF \
       -D ENABLE_Structured:BOOL=OFF \
    ${AMANZI_SRC_DIR}

    make -j8
    make install
    ```
  c. Enter your T&E for the week.  This isn't as long (~3-5 minutes?)
  d. Run Amanzi tests.  Some may not pass, but having them run means you should be in good shape.
    ```
    cd $AMANZI_BUILD_DIR
    make test
    ```

4. Install ATS.
  a. Edit another script to build ATS.  (I call this one configure-ats-debug.sh).  Add something like:
    ```
    #!/bin/bash
    # configure-ats.sh
    # script for running CMake to build ATS.

    rm -rf ${ATS_BUILD_DIR}
    mkdir -p ${ATS_BUILD_DIR}
    cd ${ATS_BUILD_DIR}

    cmake \
      -D Amanzi_DIR=${AMANZI_DIR}/lib \
      -D CMAKE_INSTALL_PREFIX=${ATS_DIR} \
    ${ATS_SRC_DIR}

    make -j8
    make install
    ```
  b. Twiddle your thumbs.  ~1-5 minutes?

5. Download ats-demos to get at some examples, and run one!
  a. Get the example repo from http://github.com/amanzi/ats-demos
  b. Run a test problem.
    ```
    cd ${ATS_BASE}/ats-demos/richards-steadystate
    mkdir run0
    cd run0
    ${ATS_DIR}/bin/ats --xml_file=../richards-steadystate.xml
    ```
  c. Run as a real test, visualized the results.
    ```
    cd ${ATS_BASE}/ats-demos
    ipython regression-tests.py -n richards-steadystate/richards-steadystate.cfg
    cd richards-steadystate
    ipython -notebook richards-steadystate.ipynb
    ```
   
