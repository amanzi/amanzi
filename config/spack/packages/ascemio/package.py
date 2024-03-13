# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Ascemio(CMakePackage):

    """ ASCEM-IO Parallel I/O module for Environmental Management Applications
    This library provides software that read/write data sets from/to parallel
    file systems in an efficient and scalable manner. """

    homepage = "https://github.com/amanzi/amanzi-tpls"
    url      = "https://github.com/amanzi/amanzi-tpls/"\
        "raw/master/src/ascem-io-2.3.tar.gz"

    version('2.3', sha256='de556a774b4ef7dc223f4611a39c978c')

    patch('../../../SuperBuild/templates/ascemio-2.2-hdf5.patch')

    variant("shared", default=True, description="Builds a shared version of the library")

    depends_on('mpi')
    depends_on('hdf5 +hl+mpi+fortran+shared',when="+shared")
    depends_on('hdf5 +hl+mpi+fortran~shared',when="-shared")


    def cmake_args(self):

        options = ['-DCMAKE_C_COMPILER=' + self.spec['mpi'].mpicc]
        options.append('-DCMAKE_CXX_COMPILER=' + self.spec['mpi'].mpicxx)
        options.append('-DCMAKE_Fortran_COMPILER=' + self.spec['mpi'].mpifc)

        return options
