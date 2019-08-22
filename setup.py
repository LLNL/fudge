# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
        'fudge.evaluationBuilder',\
        'fudge.gnd', \
        'fudge.gnd.channelData', \
        'fudge.gnd.covariances',\
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
