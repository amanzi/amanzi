# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Exprtk(Package):

    """ C++ Mathematical Expression Toolkit Library """

    homepage = "https://www.partow.net/programming/exprtk/index.html"
    git      = "https://github.com/ArashPartow/exprtk"

    maintainers = ['julienloiseau','jd-moulton','arashpartow']
    version('release', commit='a4b17d543f072d2e3ba564e4bc5c3a0d2b05c338', default=True)

    def install(self, prec, prefix):
        mkdirp(prefix.include)
        install("exprtk.hpp", prefix.include)
