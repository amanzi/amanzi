# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
from spack.pkg.builtin.petsc import Petsc


class Petsc(Petsc):

    patch('../../../SuperBuild/templates/petsc-cmake.patch')

    # modified hack from https://github.com/spack/spack/issues/27508
    @run_before('configure')
    def fixup_bug(self):
        spack.pkg.builtin.petsc.python = python
