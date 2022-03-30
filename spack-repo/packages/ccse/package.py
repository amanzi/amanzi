# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Ccse(CMakePackage):

    """ structured AMR package """

    homepage = "https://github.com/amanzi/amanzi-tpls"
    url      = "https://github.com/amanzi/amanzi-tpls/"\
        "raw/master/src/ccse-1.3.4.tar.gz"

    version('1.3.4', sha256='faa52bb553cea8ca9ea436c1a7135b12')

    variant('shared', default=True, description='Build shared library')

    patch('ccse-1.3.4-dependency.patch')
    patch('ccse-1.3.4-tools-compilers.patch')
    patch('ccse-1.3.4-tools-plot1d.patch')
    patch('ccse-1.3.5-cmake.patch')
    patch('ccse-1.3.5-rvalue.patch')
    patch('ccse-16.10-f90.patch')
    patch('ccse-mpi4.patch')

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
