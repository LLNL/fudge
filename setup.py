# <<BEGIN-copyright>>
# <<END-copyright>>

import glob, os
from distutils.core import setup, Extension, run_setup

try:
    import numpy
    NUMPYSETUPOK = True
    # find numpy include path:
    numpyPath = os.path.split( numpy.__file__ )[0]
    numpyPath = os.path.join( numpyPath, 'core/include/numpy' )

except ImportError:
    print "Warning: numpy was not found! Extensions requiring numpy will not be built"
    NUMPYSETUPOK = False

import fudge
setup( 
    name = 'fudge',\
    version = fudge.__version__,\
    author = 'Computational Nuclear Physics Group, LLNL',\
    author_email  = 'mattoon1@llnl.gov',\
    maintainer_email = 'mattoon1@llnl.gov',\
    packages = [ \
        'fudge', \
        'fudge.core', \
        'fudge.core.math', \
        'fudge.core.math.test', \
        'fudge.core.math.xData', \
        'fudge.core.utilities', \
        'fudge.core.utilities.test', \
        'fudge.gnd', \
        'fudge.gnd.channelData', \
        'fudge.gnd.productData', \
        'fudge.gnd.productData.distributions', \
        'fudge.gnd.productData.distributions.test', \
        'fudge.gnd.reactionData', \
        'fudge.gnd.reactionData.test', \
        'fudge.gnd.reactions', \
        'fudge.gnd.test', \
        'fudge.legacy', \
        'fudge.legacy.converting', \
        'fudge.legacy.endl', \
        'fudge.legacy.endl.test', \
        'fudge.particles', \
        'fudge.processing', \
        'fudge.processing.deterministic', \
        'fudge.processing.montecarlo', \
        'fudge.processing.resonances', \
        'fudge.structure', \
        'fudge.vis', \
        'fudge.vis.gnuplot', \
        'fudge.vis.matplotlib', \
        'pqu', \
        'pqu.test' \
    ],\
    package_data = { \
        'fudge.legacy.endl.test': [ 'testdb/ascii/yi01/za001001/y*','testdb/ascii/yi01/za001001/*.txt', 'testdb/ascii/yi01/za001001/*xml' ], \
    }, \
    ext_modules=[
        Extension( 'fudge.processing.resonances._getScatteringMatrices',
            sources = ['fudge/processing/resonances/getScatteringMatrices.c'],
            include_dirs = [ numpyPath ], ),
        Extension( 'fudge.processing.resonances._getCoulombWavefunctions',
            sources = ['fudge/processing/resonances/getCoulombWavefunctions.c', 'fudge/processing/resonances/coulfg2.c'],
            include_dirs = [ numpyPath ], ),
    ], \
    url = 'http://nuclear.llnl.gov/fudge',\
    license = open( 'LICENSE.txt' ).read(),\
    description = '',\
    long_description = open( 'README.txt' ).read()\
)

# Also call the setup.py in externals packages
curdir = os.path.realpath( os.curdir )
for extension in ('numericalFunctions', 'crossSectionAdjustForHeatedTarget'):
    os.chdir( extension )
    run_setup('setup.py', ['--quiet','build'])
    os.chdir( curdir )
