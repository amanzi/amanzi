= Amanzi and ATS Installation Guide =

First, a word of warning -- please be patient.  Amanzi and therefore ATS depend upon a lot of third party libraries.  This allows us to use existing, mature code to make ATS a much better software tool.  It also means installing the code and its dependencies can be quite painful and time/labor intensive.  Installation time is a very bimodal distribution -- if it "just works" this process will take 25 minutes.  If it doesn't "just work" it can take much longer.

The basics of the process are these.

0. Install pre-reqs.
1. Download amanzi-ats.
2. Use SuperBuild to install TPLs.
3. Install and test Amanzi.
4. Download and install ATS.
5. Download test problems, test ATS.

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
    On Mac, install mpi through ```homebrew``` or other package management systems.
    
  c. Lapack (linear algebra package)
    On ubuntu, this is available with:
    ```
    sudo apt-get install libblas3gf libblas-doc libblas-dev liblapack3gf liblapack-doc liblapack-dev # are there all needed?
    ```
    On Mac, install lapack through ```homebrew``` or other package management systems, i.e.
    ```
    brew install lapack
    ```


1. Set up your directory structure.  Life gets much easier if you get this right first.
  a. I use the following:
    ```
    mkdir /my/path/to/atsl
    export ATS_BASE=/my/path/to/ats
    export ATS_SRC_DIR=${ATS_BASE}/repos/ats
    export AMANZI_SRC_DIR=${ATS_BASE}/repos/amanzi
    export ATS_BUILD_DIR=${ATS_BASE}/ats-build
    export AMANZI_BUILD_DIR=${ATS_BASE}/amanzi-build    
    export AMANZI_TPLS_BUILD_DIR=${ATS_BASE}/amanzi-tpls-build    
    export ATS_DIR=${ATS_BASE}/ats-install
    export AMANZI_DIR=${ATS_BASE}/amanzi-install    
    export AMANZI_TPLS_DIR=${ATS_BASE}/amanzi-tpls-install
    export PATH=${ATS_DIR}/bin:${AMANZI_TPLS_DIR}/bin:${PATH}
    ```    
    You may want to put that line in your ~/.bashrc or similar files (~/.bash_profile on Mac OS X).
  b. Clone the Amanzi source for the latest release.  Currently this is ``0.86``
    ```git clone -b ats-amanzi-0.86 http://github.com/amanzi/amanzi $AMANZI_SRC_DIR```

  c. Clone the ATS source for the latest release.
    ```git clone -b ats-0.86 http://github.com/amanzi/ats $ATS_SRC_DIR```

  d. Set up a directory for configuration scripts and a TPL installation directory.
    ```
    cd ${ATS_BASE}
    mkdir ats-config-files  # this will hold scripts that call cmake to build things
    ```


2. Build the Amanzi TPLs.
  a. Write a configure script for SuperBuild / Amanzi TPLs.  Edit a file (I call it $ATS_BASE/ats-config-files/configure-amanzi-tpls-debug.sh).  Add something like the following text:

    ```
    #!/bin/bash
    # configure-amanzi-tpls.sh
    # script for running CMake to build TPLs for Amanzi.

    rm -rf ${AMANZI_TPLS_BUILD_DIR}
    rm -rf ${AMANZI_TPLS_DIR}
    mkdir -p ${AMANZI_TPLS_BUILD_DIR}
    mkdir -p ${AMANZI_TPLS_DIR}
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

  b. Run the script.
    ```. $ATS_BASE/ats-config-files/configure-amanzi-tpls-debug.sh```
  c. Go get a cup of coffee.  If this works, it will take ~15 minutes, depending upon your machine.  When you come back, it should tell you if it passed or not.


3. Build Amanzi.
  a. Edit another script to build Amanzi.  (I call this one $ATS_BASE/ats-config-files/configure-amanzi-debug.sh).  Add something like:
    ```
    #!/bin/bash
    # configure-amanzi.sh
    # script for running CMake to build Amanzi.

    rm -rf ${AMANZI_BUILD_DIR}
    mkdir -p ${AMANZI_BUILD_DIR}
    rm -rf ${AMANZI_DIR}
    mkdir -p ${AMANZI_DIR}
    cd ${AMANZI_BUILD_DIR}

    cmake \
       -C  ${AMANZI_TPLS_DIR}/share/cmake/amanzi-tpl-config.cmake \
       -D CMAKE_INSTALL_PREFIX=${AMANZI_DIR} \
       -D ENABLE_Physics:BOOL=OFF \
       -D ENABLE_Structured:BOOL=OFF \
       -D CMAKE_BUILD_TYPE=Debug \
    ${AMANZI_SRC_DIR}

    make -j8
    make install
    ```
  b. Run the script.
    ```. $ATS_BASE/ats-config-files/configure-amanzi-debug.sh```
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
    rm -rf ${ATS_DIR}
    mkdir -p ${ATS_DIR}
    cd ${ATS_BUILD_DIR}

    cmake \
      -D Amanzi_DIR=${AMANZI_DIR}/lib \
      -D CMAKE_INSTALL_PREFIX=${ATS_DIR} \
       -D CMAKE_BUILD_TYPE=Debug \
    ${ATS_SRC_DIR}

    make -j8
    make install
    ```
  b. Twiddle your thumbs.  ~1-5 minutes?

5. Download ats-demos to get at some examples, and run one!
  a. Get the example repo:
  ```
  cd $ATS_BASE
  mkdir testing
  cd testing
  git clone -b ats-demos-0.86 http://github.com/amanzi/ats-demos
  ```
  b. Run a test problem.
    ```
    cd ats-demos
    python regression_tests.py -n richards-steadystate
    ```
  c. Visualize the results
    ```
    cd richards-steadystate
    jupyter notebook richards-steadystate.ipynb
    ```
   
