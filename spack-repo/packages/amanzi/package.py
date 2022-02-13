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

    # Submodule is ON for ATS 
    version('spack', branch='spack', default=True, submodules=True)
    version('master', branch='master', default=True, submodules=True)
    version('1.1-dev', tag='amanzi-1.1-dev', submodules=True)
    version('1.0.0', tag='amanzi-1.0.0', submodules=True)

    # Trilinos framework: epetra|tpetra
    variant('eptra', default=True, description='Enable Epetra support')
    variant('tpetra', default=False, description='Enable Tpetra support')

    # Mesh Type: unstructured|structured
    #            (could be both, but currently not support for structured) 
    variant('mesh_type', default='unstructured', 
        values=('unstructured', 'structured'),
        description='Select mesh type: unstructured or structured', 
        multi=False)
    variant('mstk', default=True, description='Enable MSTK mesh support for '
            'unstructured mesh') 
    variant('moab', default=False, description='Enable MOAB mesh support for '
            'unstructured mesh')
    variant('silo', default=False, description='Enable Silo reader for binary '
            'files')
    
    variant('alquimia', default=False, description='Enable alquimia support')
    variant('hypre', default=True, description='Enable Hypre solver support')
    variant('ats', default=False, description='Enable ATS support')
    variant('AmanziPhysics', default=True, description='Enable Amanzi Physics support')
    variant('ATSPhysics', default=False, description='Enable Amanzi Physics support')
    variant('crunchtope', default=False, description='Enable CrunchTope support')
    #variant('petsc', default=True, description='Enable PETsC support')
#    variant('tests', default=False, description='Enable the unit test suite')

    patch('exprtk.patch', when = '@1.0.0')

    #depends_on('moab', when='+moab')
    #depends_on('silo', when='+silo')
    #depends_on('hdf5@1.10.6 +hl+mpi', when='-alquimia')
    #depends_on('superlu-dist@6.0.0', when='+hypre')
    #depends_on('hypre +superlu-dist +mpi', when='+hypre')

    depends_on('git', type='build')
    depends_on('cmake@3.15:',  type='build')

    # Mandatory 
    depends_on('mpi')
    depends_on('zlib +shared')
    depends_on('metis +shared')
    depends_on('parmetis +shared')
    depends_on('seacas +shared')
    depends_on('boost@1.67.0 cxxstd=11 +program_options +shared')
    depends_on('xerces-c')
    # depends_on('cgns@develop +mpi')
    depends_on('ascemio')
    depends_on('netcdf-c +parallel-netcdf +shared')
    depends_on('unittest-cpp')
    # Alquimia
    # xsdk-0.7.0 -> alquimia v1.0.9, pflotran v3.0.2
    depends_on('alquimia@1.0.9', when='+alquimia')
    depends_on('pflotran@3.0.2', when='+alquimia')
    # pflotran depends on petsc 3.10 (or newer). 
    depends_on('petsc@3.13.3', when="+alquimia")
    # HDF5 dependency is embeded in ascemio.
    # depends_on('hdf5@1.10.6 +mpi+fortran+hl', when='+alquimia')
    # Hypre
    depends_on('superlu@5.2.2', when='+hypre')
    depends_on('superlu-dist@6.2.0', when='+hypre')
    depends_on('hypre@2.22.1 +mpi +shared', when='+hypre')
    # MSTK 
    depends_on('mstk@3.3.6 partitioner=all +seacas +parallel +shared', when='+mstk')
    depends_on('nanoflann', when='+mstk')
    # Silo
    # There finally is a new version 4.11 (but we haven't tested it yet).
    depends_on('silo@4.10.2', when='+silo')
    # Moab
    # There is a newer version 5.3.0 (but we haven't tested it yet).
    depends_on('moab@5.2.0', when='+moab')
    # Other
    depends_on('crunchtope', when='+crunchtope')
    #depends_on('trilinos@12.14.1 +pnetcdf +boost +cgns +hdf5 +metis '
    #           '+zlib +anasazi +amesos2 +epetra +ml +teuchos +superlu-dist '
    #           '+zoltan +nox +ifpack +muelu')
    depends_on('trilinos@13.0.0 +shared +boost +hdf5 '
               '+anasazi +amesos2 +epetra +ml '
               '+zoltan +nox +ifpack +muelu')

    # Conflicts 
    conflicts('+crunchtope', when="-alquimia", msg="+crunchtope needs +alquimia") 

    def cmake_args(self):
        options = ['-DCMAKE_BUILD_TYPE=debug']
        options.append('-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc)
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)

        # Provide information normally in the cache?
        options.append('-DXERCES_LIBRARY_DIR=' + self.spec['xerces-c'].prefix + '/lib')
        #options.append('-DSuperLU_DIR=' + self.spec['superlu'].prefix)
        options.append('-DTrilinos_INSTALL_PREFIX:PATH=' + self.spec['trilinos'].prefix)

        # not supported or always off/on options
        options.append('-DENABLE_OpenMP=OFF')
        options.append('-DENABLE_SPACK_BUILD=ON')

        #if '+eptra' in self.spec: 
        #    options.append('-DENABLE_EPETRA=ON')
        #else: 
        #    options.append('-DENABLE_EPETRA=OFF')

        #if '+tpetra' in self.spec: 
        #    options.append('-DENABLE_TPETRA=ON')
        #else: 
        #    options.append('-DENABLE_TPETRA=OFF')

        if '+alquimia' in self.spec:
            options.append('-DENABLE_ALQUIMIA=ON')
            options.append('-DENABLE_PETSC=ON')
            options.append('-DENABLE_PFLOTRAN=ON')
            options.append('-DPFLOTRAN_LIBRARY_DIR=' + self.spec['pflotran'].prefix + '/lib')
            options.append('-DALQUIMIA_DIR=' + self.spec['alquimia'].prefix)
        else:
            options.append('-DENABLE_ALQUIMIA=OFF')
            options.append('-DENABLE_PETSC=OFF')
            options.append('-DENABLE_PFLOTRAN=OFF')

        if '+crunchtope' in self.spec: 
            options.append('-DENABLE_CRUNCHTOPE=ON')
            options.append('-DCRUNCHTOPE_DIR=' + self.spec['crunchtope'].prefix)
        else: 
            options.append('-DENABLE_CRUNCHTOPE=OFF')

        if '+AmanziPhysics' in self.spec: 
            options.append('-DENABLE_AmanziPhysicsModule=ON')
        else: 
            options.append('-DENABLE_AmanziPhysicsModule=OFF')

        if '+ATSPhysics' in self.spec: 
            options.append('-DENABLE_ATSPhysicsModule=ON')
        else: 
            options.append('-DENABLE_ATSPhysicsModule=OFF')


        # options based on variants
        #if '+tests' in self.spec:
        #    options.append('-DENABLE_TESTS=ON')
        #    options.append('-DENABLE_UnitTest=ON')
        #else:
        options.append('-DENABLE_TESTS=ON')
        options.append('-DENABLE_UnitTest=ON')

        if '+mstk' in self.spec:
            options.append('-DMSTK_VERSION=3.3.5')
            options.append('-DENABLE_MESH_MSTK=ON')
        else:
            options.append('-DENABLE_MESH_MSTK=OFF')

        if self.spec.variants['mesh_type'].value == 'unstructured':
            options.append('-DENABLE_Unstructured=ON')
            options.append('-DENABLE_MESH_STK=OFF')
            options.append('-DENABLE_MESH_MOAB=OFF')
        else:
            options.append('-DENABLE_Unstructured=OFF')

        if self.spec.variants['mesh_type'].value == 'structured':
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

        options.append('-DENABLE_CLM=OFF')
        options.append('-DENABLE_DBC=ON')

        # Change to the type of kokkos backend
        # Need trilinos to support CUDA
        #if '+gpu' in self.spec:
        #    options.append('-DAMANZI_ARCH="Summit"')
        #else: 
        #    options.append('-DAMANZI_ARCH=\'\'')

        # unused 
        if '+moab' in self.spec:
            options.append('-DENABLE_MESH_MOAB=ON')
        else:
            options.append('-DENABLE_MESH_MOAB=OFF')
        if '+silo' in self.spec:
            options.append('-DENABLE_Silo=ON')
        else:
            options.append('-DENABLE_Silo=OFF')
        
        return options
