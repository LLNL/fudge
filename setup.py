# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, sys, glob, shutil
import setuptools
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from setuptools import setup, Extension
import subprocess

minimumNumpy = 1.15
cwd = 'file://localhost%s/' % os.getcwd()
numpyErrorMessage = f'Install numpy>={minimumNumpy} before installing FUDGE'
try:
    import numpy
    assert float('.'.join(numpy.__version__.split('.')[:2])) >= minimumNumpy, numpyErrorMessage

except (ImportError, ModuleNotFoundError):
    sys.exit(numpyErrorMessage)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        # copy C executables Merced/bin/merced and fudge/processing/deterministic/upscatter/bin/calcUpscatterKernel to Python environment bin folder
        workingFolder = os.getcwd()
        binFolder = os.path.join(sys.prefix, 'bin')
        os.chdir('Merced')
        subprocess.check_call('make -j', shell=True)
        shutil.copy('bin/merced', binFolder)
        os.chdir(workingFolder)

        os.chdir('fudge/processing/deterministic/upscatter')
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
    author='Nuclear Data and Theory Group, LLNL',
    author_email='mattoon1@llnl.gov',
    maintainer_email='mattoon1@llnl.gov',
    packages=[
        'fudge',
        'fudge.gnds',
        'fudge.core',
        'fudge.core.math',
        'fudge.core.math.test',
        'fudge.core.utilities',
        'fudge',
        'fudge.outputChannelData',
        'fudge.outputChannelData.fissionFragmentData',
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
        'isotopicAbundances',
        'isotopicAbundances.bin'
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
        'crossSectionAdjustForHeatedTarget @ %s/crossSectionAdjustForHeatedTarget#egg=crossSectionAdjustForHeatedTarget' % cwd,
        'numericalFunctions @ %s/numericalFunctions#egg=numericalFunctions' % cwd,
        'pqu @ %s/pqu#egg=pqu' % cwd,
        'xData @ %s/xData#egg=xData' % cwd,
        'PoPs @ %s/PoPs#egg=PoPs' % cwd,
        'brownies @ %s/brownies#egg=brownies' % cwd
    ],
    license = open( 'LICENSE' ).read(),
    description = '',
    long_description = open( 'README.md' ).read(), requires=['numpy'],
    cmdclass={'install': CustomInstall, 'build_ext': CustomBuildExt}
)
