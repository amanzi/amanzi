# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Amanzi(CMakePackage):

    """Amanzi, the Multi-Process HPC Simulator is a highly modular
    and extensible computational engine for simulations of flow and
    reactive transport. It is capable of simulating
    transient saturated and variably saturated (Richards) flows,
    transport with non-grid-aligned dispersion and a variety of
    reactions. In the future it will include non-isothermal,
    multi-phase multi-component, geo-mechanical, and surface water
    models. To achive this ambitious goal we are building Amanzi
    as a grass-roots collaboration of an emerging broader community
    of geoscienists, computational and computer scientists, and
    applied mathematicians. This broader community is leveraging
    its extensive experience, existing capabilities, and untapped
    advances from their areas of expertise to develop Amanzi."""

    homepage = "http://www.amanzi.github.io"
    git      = "https://github.com/amanzi/amanzi"

    maintainers = ['julienloiseau']

    #version('master', branch='master', submodules=True)
    version('tpetra', branch='tpetra')
    version('0.89', tag='amanzi-0.89-dev', default=True)

    variant('structured', default=False,
            description='Build structured mesh capability')
    # Need to have at least one mesh enable
    variant('unstructured', default=False,
            description='Build unstructured mesh capability')
    variant('mstk', default=True, description='Enable MSTK mesh support for '
            'unstructured mesh')
    variant('moab', default=False, description='Enable MOAB mesh support for '
            'unstructured mesh')
    variant('tests', default=True, description='Enable the unit test suite')
    # Should always be off
    variant('silo', default=False, description='Enable Silo reader for binary '
            'files')
    variant('petsc', default=False, description='Enable PETsC support')
    variant('alquimia', default=False, description='Enable alquimia support')
    variant('hypre', default=True, description='Enable Hypre solver support')
    variant('ats', default=False, description='Enable ATS support')
    variant('gpu_support', default=False, description='Enable GPU support')

    depends_on('cuda', when='@tpetra +gpu_support')
    depends_on('kokkos', when='@tpetra -gpu_support')
    depends_on('kokkos +debug +cuda +cuda_uvm +cuda_relocatable_device_code ' 
        '+cuda_lambda +power9 +volta70 +serial', when='@tpetra +gpu_support')

    depends_on('git', type='build')
    depends_on('cmake@3.15:',  type='build')
    depends_on('silo', when='+silo')
    depends_on('mpi')
    depends_on('zlib')
    depends_on('hdf5@1.10.5 +hl+mpi')
    depends_on('superlu')
    #depends_on('superlu-dist', when='+hypre')
    #depends_on('hypre +superlu-dist', when='+hypre')
    depends_on('metis')
    depends_on('parmetis')
    depends_on('moab', when='+moab')
    depends_on('seacas')
    depends_on('mstk partitioner=all +exodusii +parallel', when='+mstk')
    depends_on('nanoflann', when='+mstk')
    depends_on('boost@1.59.0: cxxstd=11 +program_options')
    depends_on('cgns@develop')
    depends_on('netcdf-c +parallel-netcdf')
    depends_on('xerces-c')
    depends_on('unittest-cpp', when='+tests')
    depends_on('ascemio')
    depends_on('petsc', when='+petsc')
    depends_on('alquimia', when='+alquimia')

    depends_on('trilinos +pnetcdf +boost +cgns +hdf5 +hypre +metis '
               '+superlu-dist +zlib +amesos2 +anasazi +epetra +ml +kokkos '
               '+teuchos +tpetra +zoltan +nox +ifpack ', when='@tpetra')
    depends_on('trilinos +pnetcdf +boost +cgns +hdf5 +hypre +metis '
               '+zlib +amesos2 +anasazi +epetra +ml +teuchos +superlu-dist '
               '+zoltan +nox +ifpack +muelu', when='@0.89')

    def cmake_args(self):
        options = ['-DCMAKE_BUILD_TYPE=debug']
        options.append('-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc)
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)
        # not supported or always off/on options
        options.append('-DENABLE_OpenMP=OFF')
        options.append('-DENABLE_STK_Mesh=OFF')
        options.append('-DSPACKAGE_BUILD=ON')
        # options based on variants
        if '+tests' in self.spec:
            options.append('-DENABLE_UnitTest=ON')
        else:
            options.append('-DENABLE_UnitTest=OFF')
        if '+petsc' in self.spec:
            options.append('-DENABLE_PETSC=ON')
        else:
            options.append('-DENABLE_PETSC=OFF')
        if '+mstk' in self.spec:
            options.append('-DENABLE_MSTK_Mesh=ON')
        else:
            options.append('-DENABLE_MSTK_Mesh=OFF')
        if '+moab' in self.spec:
            options.append('-DENABLE_MOAB_Mesh=ON')
        else:
            options.append('-DENABLE_MOAB_Mesh=OFF')
        if '+silo' in self.spec:
            options.append('-DENABLE_Silo=ON')
        else:
            options.append('-DENABLE_Silo=OFF')
        if '+unstructured' in self.spec:
            options.append('-DENABLE_Unstructured=ON')
        else:
            options.append('-DENABLE_Unstructured=OFF')
        if '+structured' in self.spec:
            options.append('-DENABLE_Structured=ON')
        else:
            options.append('-DENABLE_Structured=OFF')
        if '+ascemio' in self.spec:
            options.append('-DENABLE_ASCEMIO=ON')
        else:
            options.append('-DENABLE_ASCEMIO=OFF')
        #if '+gpu' in self.spec:
        #    options.append('-DAMANZI_ARCH="Summit"')
        #else: 
        #    options.append('-DAMANZI_ARCH=\'\'')
        return options
