# Copyright 2010-202x held jointly by LANL, ORNL, LBNL, and PNNL.
# Amanzi is released under the three-clause BSD License.
# The terms of use and "as is" disclaimer for this license are
# provided in the top-level COPYRIGHT file.
#
# Author: Ethan Coon
#

from setuptools import setup, find_packages

setup(
    name='amanzi_xml',
    version='1.0.0',
    url='https://github.com/amanzi/amanzi',
    author='Amanzi-ATS Development Team',
    author_email='ats-users@googlegroups.com',
    description='Input spec manipulation tools for Amanzi and ATS input xml files.',
    packages=find_packages(),
    )

