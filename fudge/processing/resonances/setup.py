# <<BEGIN-copyright>>
# <<END-copyright>>

import os
import numpy
from distutils.core import setup, Extension

# find numpy include path:
numpyPath = os.path.split( numpy.__file__ )[0]
numpyPath = os.path.join( numpyPath, 'core/include/numpy' )

getScatteringMatrices = Extension( '_getScatteringMatrices',
        sources = ['getScatteringMatrices.c'],
        include_dirs = [ './', numpyPath ] )

getCoulombWavefunctions = Extension( '_getCoulombWavefunctions',
        sources = ['getCoulombWavefunctions.c','coulfg2.c'],
        include_dirs = [ './', numpyPath ] )

setup(name='extensions',
        version='1.0',
        description = 'Extensions (written in c) for better performance in reconstructing resonances',
        ext_modules=[ getScatteringMatrices, getCoulombWavefunctions ] )

