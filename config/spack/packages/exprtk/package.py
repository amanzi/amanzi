# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Exprtk(Package):

    """ C++ Mathematical Expression Toolkit Library """

    homepage = "https://github.com/ArashPartow/exprtk"
    git      = "https://github.com/ArashPartow/exprtk"

    maintainers = ['julienloiseau','jd-moulton','gcapodag']
    version('master', commit='806c519c91fd08ba4fa19380dbf3f6e42de9e2d1', default=True)

    def install(self, prec, prefix): 
        mkdirp(prefix.include)
        install("exprtk.hpp", prefix.include)