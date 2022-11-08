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

    version("3.0.2", commit="9e07f416a66b0ad304c720b61aa41cba9a0929d5")  # tag v3.0.2

    depends_on('mpi')
    depends_on('hdf5@1.8.12:+mpi+fortran+hl')
    depends_on('petsc')

