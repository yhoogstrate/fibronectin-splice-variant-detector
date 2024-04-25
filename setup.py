#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]

    fibronectin-splice-variant-detector: counts Fibronectin (FN1) alt. splicing in BAM files
    Copyright (C) 2024  Youri Hoogstrate, Tobias Weiss and Pim French

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/fibronectin-splice-var-determiner>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""


import fn1splicevardeterminer
from setuptools import setup


def get_requirements():
    with open('requirements.txt', 'r') as fh:
        content = fh.read().strip().split()
    
    return content


setup(name="fibronectin-splice-var-determiner",
      scripts=['bin/fibronectin-splice-var-determiner'],
      packages=["fn1splicevardeterminer"],
      test_suite="tests",
      tests_require=['nose', 'pytest', 'pytest-cov'],
      #setup_requires=['scipy', 'numpy'],
      install_requires=[get_requirements()],
      version=fn1splicevardeterminer.__version__,
      description="Free open-source tool to extract counts of splice variants of FN1",
      author=fn1splicevardeterminer.__author__,
      url=fn1splicevardeterminer.__homepage__,
      keywords=["rna-seq", "fibronectin", "FN1"],
      classifiers=[
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ])
