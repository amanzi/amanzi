# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *
from spack.pkg.builtin.trilinos import Trilinos


class Trilinos(Trilinos):
    version('13.0.0.afc4',commit='afc4e52595ab82f449f8a4676febbcfbf8223afc')

    patch('../../../SuperBuild/templates/trilinos-superludist.patch',when="@13.0.0.afc4")
    patch('../../../SuperBuild/templates/trilinos-ifpack.patch',when="@13.0.0.afc4")
    patch('../../../SuperBuild/templates/trilinos-duplicate-parameters.patch',when="@13.0.0.afc4")
