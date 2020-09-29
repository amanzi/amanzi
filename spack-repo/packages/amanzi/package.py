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
    version('master', branch='master', default=True)
    version('1.1-dev', tag='amanzi-1.1-dev')
    version('1.0.0', tag='amanzi-1.0.0')

    variant('structured', default=False,
            description='Build structured mesh capability')
    # Need to have at least one mesh enable
    variant('unstructured', default=True,
            description='Build unstructured mesh capability')
    variant('mstk', default=True, description='Enable MSTK mesh support for '
            'unstructured mesh')
#    variant('moab', default=False, description='Enable MOAB mesh support for '
#            'unstructured mesh')
    variant('tests', default=False, description='Enable the unit test suite')
#    variant('silo', default=False, description='Enable Silo reader for binary '
#            'files')
    #variant('petsc', default=True, description='Enable PETsC support')
    variant('alquimia', default=False, description='Enable alquimia support')
    variant('hypre', default=True, description='Enable Hypre solver support')
    variant('ats', default=False, description='Enable ATS support')

    patch('exprtk.patch', when = '@1.0.0')

#    depends_on('moab', when='+moab')
#    depends_on('silo', when='+silo')
    #depends_on('hdf5@1.10.6 +hl+mpi', when='-alquimia')
    #depends_on('superlu-dist@6.0.0', when='+hypre')
    #depends_on('hypre +superlu-dist +mpi', when='+hypre')


    depends_on('git', type='build')
    depends_on('cmake@3.15:',  type='build')

    # Mandatory 
    depends_on('mpi')
    depends_on('zlib')
    depends_on('metis')
    depends_on('parmetis')
    depends_on('seacas')
    depends_on('boost@1.59.0: cxxstd=11 +program_options')
    depends_on('xerces-c')
    depends_on('cgns@develop +mpi +parallel')
    depends_on('ascemio')
    depends_on('netcdf-c@4.7.3 +parallel-netcdf')
    # Alquimia
    depends_on('petsc@3.10.2')
    depends_on('hdf5@1.10.6 +mpi+fortran+hl', when='+alquimia')
    depends_on('alquimia@xsdk-0.4.0', when='+alquimia')
    depends_on('pflotran@xsdk-0.4.0')
    # Hypre
    depends_on('superlu', when='+hypre')
    depends_on('superlu-dist@5.4.0', when='+hypre')
    depends_on('hypre@2.18.2 +mpi', when='+hypre')
    # MSTK 
    depends_on('mstk@3.3.5 partitioner=all +exodusii +parallel', when='+mstk')
    depends_on('nanoflann', when='+mstk')
    # Other
    depends_on('unittest-cpp', when='+tests')

    #depends_on('trilinos@12.14.1 +pnetcdf +boost +cgns +hdf5 +metis '
    #           '+zlib +anasazi +amesos2 +epetra +ml +teuchos +superlu-dist '
    #           '+zoltan +nox +ifpack +muelu')
    depends_on('trilinos@12.18.1 +pnetcdf +boost +cgns +hdf5 +metis '
               '+zlib +anasazi +amesos2 +epetra +ml +teuchos '
               '+zoltan +nox +ifpack +muelu')

    def cmake_args(self):
        options = ['-DCMAKE_BUILD_TYPE=debug']
        options.append('-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc)
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)

        # not supported or always off/on options
        options.append('-DENABLE_OpenMP=OFF')
        options.append('-DSPACKAGE_BUILD=ON')
        
        options.append('-DXERCES_LIBRARY_DIR=' + self.spec['xerces-c'].prefix + '/lib')
        #options.append('-DSuperLU_DIR=' + self.spec['superlu'].prefix)
        options.append('-DTrilinos_INSTALL_PREFIX:PATH=' + self.spec['trilinos'].prefix)

        if '+alquimia' in self.spec:
            options.append('-DENABLE_ALQUIMIA=ON')
            options.append('-DENABLE_PETSC=ON')
            options.append('-DENABLE_PFLOTRAN=ON')
            options.append('-DENABLE_CRUNCHTOPE=ON')
        else:
            options.append('-DENABLE_ALQUIMIA=OFF')
            options.append('-DENABLE_PETSC=OFF')
            options.append('-DENABLE_PFLOTRAN=OFF')
            options.append('-DENABLE_CRUNCHTOPE=OFF')

        # options based on variants
        if '+tests' in self.spec:
            options.append('-DENABLE_TESTS=ON')
            options.append('-DENABLE_UnitTest=ON')
        else:
            options.append('-DENABLE_TESTS=OFF')
            options.append('-DENABLE_UnitTest=OFF')

        if '+mstk' in self.spec:
            options.append('-DMSTK_VERSION=3.3.5')
            options.append('-DENABLE_MSTK_Mesh=ON')
        else:
            options.append('-DENABLE_MSTK_Mesh=OFF')

        if '+unstructured' in self.spec:
            options.append('-DENABLE_Unstructured=ON')
            options.append('-DENABLE_STK_Mesh=OFF')
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

        if '+hypre' in self.spec: 
            options.append('-DENABLE_SUPERLU=ON')
            options.append('-DENABLE_HYPRE=ON')
        else: 
            options.append('-DENABLE_SUPERLU=OFF')
            options.append('-DENABLE_HYPRE=OFF')

        #options.append('-DENABLE_ATSPhysicsModule=ON')
        options.append('-DENABLE_AmanziPhysicsModule=OFF')
        options.append('-DENABLE_CLM=OFF')
        options.append('-DENABLE_DBC=ON')

        # Change to the type of kokkos backend
        # Need trilinos to support CUDA
        #if '+gpu' in self.spec:
        #    options.append('-DAMANZI_ARCH="Summit"')
        #else: 
        #    options.append('-DAMANZI_ARCH=\'\'')

        # unused 
        #if '+moab' in self.spec:
        #    options.append('-DENABLE_MOAB_Mesh=ON')
        #else:
        #    options.append('-DENABLE_MOAB_Mesh=OFF')
        #if '+silo' in self.spec:
        #    options.append('-DENABLE_Silo=ON')
        #else:
        #    options.append('-DENABLE_Silo=OFF')
        
        return options
