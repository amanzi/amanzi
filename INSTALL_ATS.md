ATS Installation Guide
==================================

NOTE: This has changed significantly since ATS version 0.88.  If you are building versions of ATS including version 0.88 or older, do not follow these instructions, but follow those at: https://github.com/amanzi/ats/blob/ats-0.88/INSTALL.md .  This file describes installation ATS version 1.0 and newer.

First, a word of warning -- please be patient.  Amanzi and therefore ATS depend upon a lot of third party libraries.  This allows us to use existing, mature code to make ATS a much better software tool.  It also means installing the code and its dependencies can be quite painful and time/labor intensive.  Installation time is a very bimodal distribution -- if it "just works" this process will take 10=20 minutes.  If it doesn't "just work" it can take much longer.

The basics of the process are these.

0. Install pre-reqs.
1. Download Amanzi source.
2. Use bootstrap to build TPLs, Amanzi, and ATS in one go.
3. Download test problems, test ATS.

All instructions assume you use bash.  Change as needed for other shells.

0. Ensure you have cmake and an MPI installation.

  * cmake >= 3.10
    ```
    which cmake  # this will test if you have any cmake installed
    ``` 
    If not, go to: http://www.cmake.org/cmake/resources/software.html
    To install cmake from command line see https://cmake.org/install/
    
    Alternative, most binary installations are new enough, including Homebrew, Ubuntu, etc.
    
  * MPI
    ```
    which mpicc # see if you have any MPI installed.
    locate mpi.h # see if you have header files installed (is this enough?)
    ```
    On ubuntu, these are available with:
    ```
    sudo apt-get install libopenmpi-dev openmpi-bin  # I believe these are the sufficient conditions?
    ```
    On Mac, install mpi through ```homebrew``` or other package management systems.
    
  * Lapack (linear algebra package)
    On ubuntu, this is available with:
    ```
    sudo apt-get install libblas3gf libblas-doc libblas-dev liblapack3gf liblapack-doc liblapack-dev # are there all needed?
    ```
    On Mac, install lapack through ```homebrew``` or other package management systems, i.e.
    ```
    brew install lapack
    ```


1. Set up your directory structure.  Life gets much easier if you get this right first.

  * I use the following:
    ```
    # EDIT THESE!
    export ATS_BASE=/my/path/to/all/things/ats
    export ATS_BUILD_TYPE=Release
    export ATS_VERSION=master
    export OPENMPI_DIR=/my/path/to/mpi
    # END EDIT THESE!

    export AMANZI_TPLS_BUILD_DIR=${ATS_BASE}/amanzi_tpls-build-${ATS_VERSION}-${ATS_BUILD_TYPE}
    export AMANZI_TPLS_DIR=${ATS_BASE}/amanzi_tpls-install-${ATS_VERSION}-${ATS_BUILD_TYPE}

    export AMANZI_SRC_DIR=${ATS_BASE}/repos/amanzi
    export AMANZI_BUILD_DIR=${ATS_BASE}/amanzi-build-${ATS_VERSION}-${ATS_BUILD_TYPE}
    export AMANZI_DIR=${ATS_BASE}/amanzi-install-${ATS_VERSION}-${ATS_BUILD_TYPE}

    export ATS_SRC_DIR=${AMANZI_SRC_DIR}/src/physics/ats
    export ATS_DIR=${AMANZI_DIR}

    export PATH=${ATS_DIR}/bin:${PATH}
    export PATH=${AMANZI_TPLS_DIR}/bin:${PATH}
    export PYTHONPATH=${ATS_SRC_DIR}/tools/utils:${PYTHONPATH}
    export PYTHONPATH=${AMANZI_TPLS_DIR}/SEACAS/lib:${PYTHONPATH}
    ```    

    Note two things here -- first, you may want to build both a ATS_BUILD_TYPE=Debug and ATS_BUILD_TYPE=Release builds.  Debug is extremely useful for catching errors in your input files, while Release is significantly faster.  I HIGHLY RECOMMEND BUILDING BOTH, then using Debug until your input files work and give reasonable results, then swapping to Release to do production runs.

    You may want to put these lines in your `~/.bashrc` or similar files (`~/.bash_profile` on Mac OS X), or better yet use Environment modules.

    Then make your base directory and go there:
    ```
    mkdir -p ${ATS_BASE}
    cd ${ATS_BASE}
    ```

  * Clone the Amanzi source for the latest commit.
    ```
    git clone -b master http://github.com/amanzi/amanzi $AMANZI_SRC_DIR
    ```

2. Configure and build the Amanzi TPLs, Amanzi, and ATS.

  * Run build_ATS_generic.sh. This script simply invokes Amanzi's bootstrap.sh with appropriate options for an ATS build. If you follow the above instructions exactly, you can probably use that file as is, but you should look through it to make sure everything makes sense. More information about the build process options are available within bootstrap. Go ahead and look through it if you're interested. Any options available in bootstrap can be added to build_ATS_generic.
    ```
    . ${AMANZI_SRC_DIR}/build_ATS_generic.sh
    ```

  * Go get a cup of coffee.  This will take a bit, from 10-15 minutes on a fast workstation to much longer on a cluster.

3. Download ats-demos to get at some examples, and run one!

  * Get the example repo:
    ```
    cd $ATS_BASE
    mkdir testing
    cd testing
    git clone -b master http://github.com/amanzi/ats-demos
    ```

  * Run a test problem.
    ```
    cd ats-demos
    cd 01_richards_steadystate
    ats --xml_file=richards_steadystate.xml &> out.log
    ```

  * Visualize the results
    ```  
    jupyter notebook richards_steadystate.ipynb
    ```
   
