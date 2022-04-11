# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Ccse(CMakePackage):

    """ structured AMR post processing tools """

    homepage = "https://ccse.lbl.gov/BoxLib/"
    url      = "https://github.com/BoxLib-Codes/BoxLib/archive/17.05.1.tar.gz"

    version('17.05.1', sha256='97c3e1615cc649e2748fa9a7291724fa')

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
        # setting spacedim to 2 below 
        options.append('-DBL_SPACEDIM:INT=2')

        return options
