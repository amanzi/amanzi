# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class Pflotran(AutotoolsPackage):
    """PFLOTRAN is an open source, state-of-the-art massively parallel
       subsurface flow and reactive transport code.
    """

    homepage = "https://www.pflotran.org"
    git      = "https://bitbucket.org/pflotran/pflotran.git"

    version('0.0.1',commit='1ba735f1a64a12315c52a64107a75571c53492b3')

    depends_on('mpi')
    depends_on('hdf5@1.8.12:+mpi+fortran+hl')
    depends_on('petsc')

    @property
    def parallel(self):
        return (self.spec.satisfies('@xsdk-0.4.0:'))
