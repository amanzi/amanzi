# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class Crunchtope(CMakePackage):

    """CrunchTope """

    homepage = "https://bitbucket.org/crunchflow/crunchtope-dev/wiki/Home"
    git      = "https://bitbucket.org/crunchflow/crunchtope-dev.git"

    maintainers = ['julienloiseau']

    # Submodule is ON for ATS 
    version('spack', branch='spack', default=True, submodules=True)

    variant('shared', default=True, description='Build sharedd library')

    depends_on('git', type='build')
    depends_on('cmake@3.15:',  type='build')

    # Mandatory 
    depends_on('mpi')
    depends_on('petsc')
    depends_on('hdf5 +mpi+fortran+hl+shared',when='+shared')
    depends_on('hdf5 +mpi+fortran+hl~shared',when='-shared')
    #depends_on('alquimia@1.0.9')
    depends_on('pflotran@3.0.2')

    def cmake_args(self):
        options = ['-DCMAKE_BUILD_TYPE=%s' % self.spec.variants["build_type"].value]
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)
        options.append('-DCMAKE_Fortran_FLAGS=-w -DALQUIMIA -Wall -Wno-unused0variable -ffree-line-length-0 -O3')
        options.append('-DBUILD_SHARED_LIBS:BOOL=TRUE')
        options.append('-DTPL_PETSC_LIBRARIES:PATH='+self.spec['petsc'].prefix+'/lib')
        options.append('-DTPL_PETSC_INCLUDE_DIRS:PATH='+self.spec['petsc'].prefix+'/include')
        options.append('-DPETSC_ARCH=.')
        return options

