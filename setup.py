# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, glob, shutil
import setuptools
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from setuptools import setup, Extension
import subprocess

minimumNumpy = 1.15
try:
    import numpy
    assert float('.'.join(numpy.__version__.split('.')[:2])) >= minimumNumpy, numpyErrorMessage

except (ImportError, ModuleNotFoundError):
    sys.exit('Install numpy>=%s before installing FUDGE' % minimumNumpy)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        # install submodules
        #subprocess.call('python ./crossSectionAdjustForHeatedTarget/setup.py install', shell=True)
        #subprocess.call('python ./numericalFunctions/setup.py install', shell=True)
        #subprocess.call('python ./pqu/setup.py install', shell=True)
        #subprocess.call('python ./xData/setup.py install', shell=True)
        #subprocess.call('python ./PoPs/setup.py install', shell=True)
        #subprocess.call('python ./brownies/setup.py install', shell=True)

        # copy C executables Merced/bin/merced and upscatter/bin/calcUpscatterKernel to Python environment bin folder
        workingFolder = os.getcwd()
        binFolder = os.path.join(sys.prefix, 'bin')
        os.chdir('Merced')
        subprocess.check_call('make -j', shell=True)
        shutil.copy('bin/merced', binFolder)
        os.chdir(workingFolder)

        os.chdir('upscatter')
        subprocess.check_call('make -j', shell=True)
        shutil.copy('bin/calcUpscatterKernel', binFolder)
        os.chdir(workingFolder)

        super().run()


class CustomBuildExt(build_ext):
    def run(self):
        # find numpy include path:
        numpyPath = os.path.split( numpy.__file__ )[0]
        numpyPath = os.path.join( numpyPath, 'core/include/numpy' )
        assert os.path.isdir(numpyPath), 'Numpy path "%s" NOT FOUND' % numpyPath

        for ext in self.extensions:
            ext.include_dirs.append(numpyPath)

        super().run()


import getFudgeVersion as versionModule

setup(
    name='fudge',
    #    version = fudge.__version__,\
    version=versionModule.getVersionNumber(),
    author='Computational Nuclear Physics Group, LLNL',
    author_email='mattoon1@llnl.gov',
    maintainer_email='mattoon1@llnl.gov',
    packages=[
        'fudge',
        'fudge.core',
        'fudge.core.math',
        'fudge.core.math.test',
        'fudge.core.utilities',
        'fudge',
        'fudge.channelData',
        'fudge.channelData.fissionFragmentData',
        'fudge.covariances',
        'fudge.covariances.test',
        'fudge.productData',
        'fudge.productData.distributions',
        'fudge.productData.distributions.test',
        'fudge.reactionData',
        'fudge.reactionData.test',
        'fudge.reactionData.doubleDifferentialCrossSection',
        'fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic',
        'fudge.reactionData.doubleDifferentialCrossSection.photonScattering',
        'fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw',
        'fudge.reactions',
        'fudge.resonances',
        'fudge.processing',
        'fudge.processing.deterministic',
        'fudge.processing.montecarlo',
        'fudge.processing.resonances',
        'fudge.processing.resonances.test',
        'fudge.vis',
        'fudge.vis.gnuplot',
        'fudge.vis.matplotlib',
        'LUPY',        
        #'crossSectionAdjustForHeatedTarget',
        #'numericalFunctions',
        #'pqu',
        #'xData',
        #'PoPs',
        #'brownies'
    ],
    package_dir = {'': '.'},
    scripts = glob.glob('bin/*.py'),
    package_data = {
        'fudge.legacy.endl.test': [ 'testdb/ascii/yi01/za001001/y*', 'testdb/ascii/yi01/za001001/*.txt', 'testdb/ascii/yi01/za001001/*xml' ],
        'fudge.processing.resonances.test': ['*.py'],
        'fudge.legacy.endl': ['bdfls'],
        'fudge': ['gnds.xsd'],
        'fudge.covariances.test': ['*.py'],
    },
    ext_modules=[
        Extension( 'fudge.processing.resonances._getBreitWignerSums',
            sources = ['fudge/processing/resonances/getBreitWignerSums.c'], ),
        Extension( 'fudge.processing.resonances._getScatteringMatrices',
            sources = ['fudge/processing/resonances/getScatteringMatrices.c'], ),
        Extension( 'fudge.processing.resonances._getCoulombWavefunctions',
            sources = ['fudge/processing/resonances/getCoulombWavefunctions.c', 'fudge/processing/resonances/coulfg2.c'], ),
    ],
    url = 'https://github.com/llnl/fudge',
        install_requires=[
        'numpy', 
        'crossSectionAdjustForHeatedTarget @ file://localhost/./crossSectionAdjustForHeatedTarget'
    ],
    license = open( 'LICENSE' ).read(),
    description = '',
    long_description = open( 'README.md' ).read(), requires=['numpy'],
    cmdclass={'install': CustomInstall, 'build_ext': CustomBuildExt}
)

# Also call the setup.py in externals packages
# curdir = os.path.realpath( os.curdir )
# for extension in ('numericalFunctions', 'crossSectionAdjustForHeatedTarget'): #'Merced', 'statusMessageReporting'):
#     os.chdir( extension )
#     run_setup('setup.py', ['--quiet','build'])
#     os.chdir( curdir )
