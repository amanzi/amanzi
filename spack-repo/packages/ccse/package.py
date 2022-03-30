# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class ccse(CMakePackage):

    """ structured AMR package """

    homepage = "https://github.com/amanzi/amanzi-tpls"
    url      = "https://github.com/amanzi/amanzi-tpls/"\
        "raw/master/src/ccse-1.3.4.tar.gz"

    version('1.3.4', sha256='faa52bb553cea8ca9ea436c1a7135b12')

    variant('shared', default=True, description='Build shared library')

    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-1.3.4-dependency.patch', sha256='4949e9d20031c83663db06af2a82932c')
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-1.3.4-tools-compilers.patch', sha256='7953619207e0ec5c217ff6f63338f95b' )
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-1.3.4-tools-plot1d.patch', sha256='69e412a72adcd8fc97f5bdf843f28de9')
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-1.3.5-cmake.patch', sha256='3ec6e137675ca6b4afff62b5394124b5')
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-1.3.5-rvalue.patch', sha256='dbcf5a13f7cbf2b96c10eb61eeafc87c')
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-16.10-f90.patch', sha256='0cc9013ebe2547b39a4e1179365d6d6d')
    patch('https://github.com/amanzi/amanzi/tree/master/config/SuperBuild/templates/ccse-mpi4.patch', sha256='ff7ea30004e91421cf1b7f051b02d774')

    depends_on('mpi')

    def cmake_args(self):
        options = ['-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc]
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)
        options.append('-DENABLE_Config_Report:BOOL=TRUE')
        options.append('-DENABLE_MPI:INT=1')
        options.append('-DENABLE_TESTS:BOOL=FALSE')
        options.append('-DVERBOSE:BOOL=ON')
        # setting spacedim to 2 below 
        options.append('-DBL_SPACEDIM:INT=2')
        if '+shared' in self.spec:
            options.append('-DBUILD_SHARED_LIBS:BOOL=TRUE')
        else :
            options.append('-DBUILD_SHARED_LIBS:BOOL=FALSE')
        return options
