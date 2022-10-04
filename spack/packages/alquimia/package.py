# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Alquimia(CMakePackage):
    """Alquimia is an interface that exposes the capabilities
    of mature geochemistry codes such as CrunchFlow and PFLOTRAN"""

    homepage = "https://github.com/LBL-EESA/alquimia-dev"
    git      = "https://github.com/LBL-EESA/alquimia-dev.git"

    maintainers = ['smolins', 'balay']

    version('1.0.9', commit='2ee3bcfacc63f685864bcac2b6868b48ad235225')  # tag v.1.0.9

    variant('shared', default=True,
            description='Enables the build of shared libraries')

    depends_on('crunchtope')
    depends_on('mpi')
    depends_on('hdf5')
    depends_on('pflotran@3.0.2', when='@1.0.9')
    depends_on('petsc')

    patch('alquimia-cmake.patch')
    patch('alquimia-FindPETSc.patch')

    def cmake_args(self):
        spec = self.spec

        options = ['-DCMAKE_C_COMPILER=%s' % spec['mpi'].mpicc,
                   '-DCMAKE_BUILD_TYPE=%s' % spec.variants['build_type'].value,
                   '-DCMAKE_Fortran_COMPILER=%s' % spec['mpi'].mpifc,
                   '-DUSE_XSDK_DEFAULTS=YES',
                   self.define_from_variant('BUILD_SHARED_LIBS', 'shared'),
                   '-DTPL_ENABLE_MPI:BOOL=ON',
                   '-DMPI_BASE_DIR:PATH=%s' % spec['mpi'].prefix,
                   '-DTPL_ENABLE_HDF5:BOOL=ON',
                   '-DXSDK_WITH_PFLOTRAN:BOOL=ON',
                   # This is not good.
                   # It assumes that the .a file exists and is not a .so
                   '-DTPL_PFLOTRAN_LIBRARIES=%s' % (
                       spec['pflotran'].prefix.lib + "/libpflotranchem.a"),
                   '-DTPL_PFLOTRAN_INCLUDE_DIRS=%s' % (
                       spec['pflotran'].prefix.include),
                   '-DPETSC_DIR=' + spec['petsc'].prefix, 
                   '-DCMAKE_INSTALL_NAME_DIR:PATH=%s/lib' % self.prefix,
                   '-DXSDK_WITH_CRUNCHFLOW:BOOL=ON',
                   '-DTPL_CRUNCHFLOW_LIBRARIES:FILEPATH=' + spec['crunchtope'].prefix + '/lib/libcrunchchem.so',
                   '-DTPL_CRUNCHFLOW_INCLUDE_DIRS:FILEPATH=' + spec['crunchtope'].prefix + '/lib/']
        return options

