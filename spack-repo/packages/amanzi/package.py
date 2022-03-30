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

    maintainers = ['julienloiseau','jd-moulton','gcapodag']

    # Submodule is ON for ATS 
    version('spack', branch='spack', default=True, submodules=True)
    version('master', branch='master', default=False, submodules=True)
    version('1.1-dev', tag='amanzi-1.1-dev', default=False, submodules=True)
    version('1.0.0', tag='amanzi-1.0.0', default=False, submodules=True)

    # Trilinos data model: epetra|tpetra
    variant('data_model', default='epetra', values=('epetra','tpetra'),
        description='Trilinos data model', multi=False)

    # Mesh Type: unstructured|structured
    #            (could be both, but currently not support for structured) 
    variant('mesh_type', default='unstructured', 
        values=('unstructured', 'structured'),
        description='Select mesh type: unstructured or structured')
    variant('shared', default=True, description='Build shared libraries and TPLs')
    variant('mesh_framework', default='mstk', values=('mstk','moab'), 
        description='Unstructure mesh framework')

    # Solvers 
    variant('hypre', default=True, description='Enable Hypre solver support')

    # Physics 
    variant('physics', default='amanzi', values=('amanzi','ats'), description='Physics implementation')

    # I/O
    variant('silo', default=False, description='Enable Silo reader for binary '
            'files')

    # Geochemistry 
    variant('geochemistry', default=False, description='Enable geochemistry support')


    variant('tests', default=True, description='Enable the unit test suite')

    patch('exprtk.patch')

    conflicts('physics=ats', when='mesh_type=structured', msg='ERROR: ats physics is not supported in the structured mesh framework.')

    ##### Build dependencies #####

    depends_on('git', type='build')
    depends_on('cmake@3.15:',  type='build')

    ##### CORE DEPENDENCIES ##### 

    depends_on('mpi')

    core_dependencies = {
        'zlib','metis', 'parmetis', 'seacas -x11',
        'boost@1.67.0 cxxstd=11 +program_options',
        'netcdf-c +parallel-netcdf', 'hdf5 +mpi+fortran+hl api=default', 
        'ascemio'
    }

    for dep in core_dependencies: 
        depends_on(dep+' +shared',when='+shared'); 
        depends_on(dep+' -shared',when='-shared');

    #double check the version of cgns  below
    depends_on('cgns@develop -shared', when='-shared')

    depends_on('ccse', when='mesh_type=structured')

    # The following core dependencies do not support +shared/~shared
    depends_on('xerces-c')
    depends_on('unittest-cpp')

    #### Geochemistry ####
    geochemistry = {
        'alquimia@1.0.9','petsc@3.16:'
    }
    for dep in geochemistry: 
        depends_on(dep+' +shared',when='+geochemistry +shared'); 
        depends_on(dep+' -shared',when='+geochemistry -shared'); 

    depends_on('pflotran@3.0.2', when='+geochemistry')

    ##### Hypre #####
    depends_on('superlu@5.2.2', when='+hypre')
    depends_on('superlu-dist@6.2.0 +shared', when='+hypre +shared')
    depends_on('superlu-dist@6.2.0 -shared', when='+hypre -shared')
    depends_on('hypre@2.22.0 +mpi +shared', when='+hypre +shared')
    depends_on('hypre@2.22.0 +mpi -shared', when='+hypre -shared')


    ##### MSTK #####
    depends_on('mstk@3.3.6 partitioner=all +seacas +parallel +shared', when='mesh_framework=mstk +shared')
    depends_on('mstk@3.3.6 partitioner=all +seacas +parallel -shared', when='mesh_framework=mstk -shared')

    depends_on('nanoflann', when='mesh_framework=mstk')

    ##### Silo #####
    # There finally is a new version 4.11 (but we haven't tested it yet).
    depends_on('silo@4.10.2 +shared', when='+silo +shared')
    depends_on('silo@4.10.2 -shared', when='+silo -shared')

    ##### Moab #####
    # There is a newer version 5.3.0 (but we haven't tested it yet).
    depends_on('moab@5.2.0 +shared', when='mesh_framework=moab +shared')
    depends_on('moab@5.2.0 -shared', when='mesh_framework=moab -shared')


    ##### Other #####
    depends_on('trilinos@13.0.0 +shared +boost +hdf5 +hypre '
               '+anasazi +amesos2 +epetra +ml '
               '+zoltan +nox +ifpack +muelu -ifpack2 cxxstd=11',when='+shared')

    depends_on('trilinos@13.0.0 -shared +boost +hdf5 +hypre '
               '+anasazi +amesos2 +epetra +ml '
               '+zoltan +nox +ifpack +muelu -ifpack2 cxxstd=11',when='-shared')

    def cmake_args(self):
        options = ['-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc]
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)

        # Tags 
        if self.spec.satisfies('@spack'):
            options.append('-DSPACK_AMANZI_VERSION_MAJOR=1')
            options.append('-DSPACK_AMANZI_VERSION_MINOR=0')
            options.append('-DSPACK_AMANZI_VERSION=1.0')
            options.append('-DSPACK_AMANZI_VERSION_PATCH=0')
            options.append('-DSPACK_AMANZI_VERSION_HASH=0')
        

        # Provide information normally in the cache?
        options.append('-DUnitTest_DIR=' + self.spec['unittest-cpp'].prefix)
        options.append('-DZLIB_DIR=' + self.spec['zlib'].prefix)
        options.append('-DMETIS_DIR=' + self.spec['metis'].prefix)
        options.append('-DBOOST_ROOT=' + self.spec['boost'].prefix)
        options.append('-DHDF5_ROOT=' + self.spec['hdf5'].prefix)
        options.append('-DHDF5_HL_LIBRARIES=' + self.spec['hdf5'].prefix + '/lib')
        options.append('-DHDF5_hdf5_hl_LIBRARY=' + self.spec['hdf5'].prefix + '/lib')
        options.append('-DASCEMIO_DIR=' + self.spec['ascemio'].prefix)
        options.append('-DNetCDF_DIR=' + self.spec['netcdf-c'].prefix)
        options.append('-DMSTK_LIBRARY_DIR=' + self.spec['mstk'].prefix + '/lib')
        options.append('-DMSTK_INCLUDE_DIR=' + self.spec['mstk'].prefix + '/include')
        options.append('-DXERCES_DIR=' + self.spec['xerces-c'].prefix)
        options.append('-DSEACAS_DIR=' + self.spec['seacas'].prefix)
        options.append('-DSuperLU_DIR=' + self.spec['superlu'].prefix)
        options.append('-DHYPRE_DIR=' + self.spec['hypre'].prefix)
        options.append('-DTrilinos_INSTALL_PREFIX:PATH=' + self.spec['trilinos'].prefix)

        # not supported or always off/on options
        options.append('-DENABLE_OpenMP=OFF')
        options.append('-DENABLE_SPACK_BUILD=ON')
        options.append('-DENABLE_ASCEMIO=ON')
        options.append('-DENABLE_CLM=OFF')
        options.append('-DENABLE_DBC=ON')

        if '+shared' in self.spec: 
            options.append('BUILD_SHARED_LIBS=ON')
        else: 
            options.append('BUILD_SHARED_LIBS=OFF')

        if 'data_model=epetra' in self.spec: 
            options.append('-DENABLE_EPETRA=ON')
        else: 
            options.append('-DENABLE_EPETRA=OFF')

        if 'data_model=tpetra' in self.spec: 
            options.append('-DENABLE_TPETRA=ON')
        else: 
            options.append('-DENABLE_TPETRA=OFF')

        if '+geochemistry' in self.spec:
            options.append('-DENABLE_ALQUIMIA=ON')
            options.append('-DENABLE_PETSC=ON')
            options.append('-DENABLE_PFLOTRAN=ON')
            options.append('-DPFLOTRAN_LIBRARY_DIR=' + self.spec['pflotran'].prefix + '/lib')
            options.append('-DAlquimia_DIR=' + self.spec['alquimia'].prefix) 
            options.append('-DAlquimia_INCLUDE_DIR=' + self.spec['alquimia'].prefix + 'include/alquimia') 
        else:
            options.append('-DENABLE_ALQUIMIA=OFF')
            options.append('-DENABLE_PETSC=OFF')
            options.append('-DENABLE_PFLOTRAN=OFF')

        if '+crunchtope' in self.spec: 
            options.append('-DENABLE_CRUNCHTOPE=ON')
            options.append('-DCRUNCHTOPE_DIR=' + self.spec['crunchtope'].prefix)
        else: 
            options.append('-DENABLE_CRUNCHTOPE=OFF')

        if 'physics=amanzi' in self.spec: 
            options.append('-DENABLE_AmanziPhysicsModule=ON')
        else: 
            options.append('-DENABLE_AmanziPhysicsModule=OFF')

        if 'physics=ats' in self.spec: 
            options.append('-DENABLE_ATSPhysicsModule=ON')
        else: 
            options.append('-DENABLE_ATSPhysicsModule=OFF')

        if '+tests' in self.spec:
            options.append('-DENABLE_TESTS=ON')
            options.append('-DENABLE_UnitTest=ON')
        else:
            options.append('-DENABLE_TESTS=OFF')
            options.append('-DENABLE_UnitTest=OFF')

        if 'mesh_framework=mstk' in self.spec:
            options.append('-DMSTK_VERSION=3.3.6')
            options.append('-DENABLE_MESH_MSTK=ON')
        else:
            options.append('-DENABLE_MESH_MSTK=OFF')

        if 'mesh_framework=moab' in self.spec:
            options.append('-DENABLE_MESH_MOAB=ON')
        else:
            options.append('-DENABLE_MESH_MOAB=OFF')
       
        if 'mesh_type=unstructured' in self.spec:
            options.append('-DENABLE_Unstructured=ON')
            options.append('-DENABLE_MESH_STK=OFF')
            options.append('-DENABLE_MESH_MOAB=OFF')
        else:
            options.append('-DENABLE_Unstructured=OFF')

        if 'mesh_type=structured'in self.spec: 
            options.append('-DENABLE_Structured=ON')
            options.append('-DENABLE_CCSE_TOOLS=ON')
            options.append('-DCCSE_DIR=' + self.spec['ccse'].prefix)
        else:
            options.append('-DENABLE_Structured=OFF')
        
        if '+hypre' in self.spec: 
            options.append('-DENABLE_SUPERLU=ON')
            options.append('-DENABLE_HYPRE=ON')
        else: 
            options.append('-DENABLE_SUPERLU=OFF')
            options.append('-DENABLE_HYPRE=OFF')

        if '+silo' in self.spec:
            options.append('-DENABLE_Silo=ON')
        else:
            options.append('-DENABLE_Silo=OFF')
        
        # Change to the type of kokkos backend
        # Need trilinos to support CUDA
        #if '+gpu' in self.spec:
        #    options.append('-DAMANZI_ARCH="Summit"')
        #else: 
        #    options.append('-DAMANZI_ARCH=\'\'')

        # unused 
        
        
        return options
