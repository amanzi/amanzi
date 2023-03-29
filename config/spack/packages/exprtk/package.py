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
    version('release', commit='f46bffcd6966d38a09023fb37ba9335214c9b959', default=True)

    def install(self, prec, prefix):
        mkdirp(prefix.include)
        install("exprtk.hpp", prefix.include)
