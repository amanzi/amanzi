# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class Ccse(CMakePackage):

    """ structured AMR post processing tools """

    homepage = "https://ccse.lbl.gov/BoxLib/"
    url      = "https://github.com/BoxLib-Codes/BoxLib/archive/17.05.1.tar.gz"

    version('17.05.1', sha256='97c3e1615cc649e2748fa9a7291724fa')
    variant('dims', default='2', values=('2', '3'), multi=False, description='Number of spatial dimensions')

    patch('../../../config/SuperBuild/templates/ccse-1.3.4-dependency.patch')
    patch('../../../config/SuperBuild/templates/ccse-1.3.4-tools-compilers.patch')
    patch('../../../config/SuperBuild/templates/ccse-1.3.4-tools-plot1d.patch')
    patch('../../../config/SuperBuild/templates/ccse-1.3.5-cmake.patch')
    patch('../../../config/SuperBuild/templates/ccse-1.3.5-rvalue.patch')
    patch('../../../config/SuperBuild/templates/ccse-16.10-f90.patch')
    patch('../../../config/SuperBuild/templates/ccse-mpi4.patch')

    depends_on('mpi')

    def cmake_args(self):
        spec = self.spec
        options = ['-DCMAKE_C_COMPILER=%s' % spec['mpi'].mpicc]
        options.append('-DCMAKE_CXX_COMPILER=%s' % spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=%s' % spec['mpi'].mpifc)
        options.append('-DBL_SPACEDIM=%d' % int(spec.variants['dims'].value))
        options.append('-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON')
        options.append('-DENABLE_FBASELIB=ON')

        return options
