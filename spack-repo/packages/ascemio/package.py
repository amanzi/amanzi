# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Ascemio(CMakePackage):

    """ ASCEM-IO Parallel I/O module for Environmental Management Applications
    This library provides software that read/write data sets from/to parallel
    file systems in an efficient and scalable manner. """

    homepage = "https://github.com/amanzi/amanzi-tpls"
    url      = "https://github.com/amanzi/amanzi-tpls/"\
        "raw/master/src/ascem-io-2.3.tar.gz"

    version('2.3', sha256='de556a774b4ef7dc223f4611a39c978c')

    patch('https://raw.githubusercontent.com/amanzi/amanzi/master/config/'\
        'SuperBuild/templates/ascemio-2.2-hdf5.patch',
        sha256='f6c14e802604b56c723a69296a77835912e8282b86b7cef84fb4f2923079bbc5',
        when='@2.3')

    depends_on('mpi')
    depends_on('hdf5@1.10.6 +hl+mpi')

    def cmake_args(self):
        options = []
        return options
