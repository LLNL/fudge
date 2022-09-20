# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import numpy
from distutils.core import setup, Extension

# find numpy include path:
numpyPath = os.path.split(numpy.__file__)[0]
numpyPath = os.path.join(numpyPath, 'core/include/numpy')

getBreitWignerSums = Extension('_getBreitWignerSums',
        sources=['getBreitWignerSums.c'],
        include_dirs=['./', numpyPath])

getScatteringMatrices = Extension('_getScatteringMatrices',
        sources=['getScatteringMatrices.c'],
        include_dirs=['./', numpyPath])

getCoulombWavefunctions = Extension('_getCoulombWavefunctions',
        sources = ['getCoulombWavefunctions.c','coulfg2.c'],
        include_dirs = ['./', numpyPath])

if __name__ == '__main__':
    setup(name='extensions',
            version='1.0',
            description='Extensions (written in c) for better performance in reconstructing resonances',
            ext_modules=[getBreitWignerSums, getScatteringMatrices, getCoulombWavefunctions])
