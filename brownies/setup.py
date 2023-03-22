# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
import glob

def setup():
    from setuptools import setup

    setup(
        name='brownies',
        packages = [
            'brownies.bin',
            'brownies.BNL',
            'brownies.BNL.utilities', 'brownies.BNL.utilities.test',
            'brownies.BNL.restools', 'brownies.BNL.restools.test',
            'brownies.BNL.plot_evaluation', 'brownies.BNL.plot_evaluation.test',
            'brownies.BNL.inter', 'brownies.BNL.inter.spectra', 'brownies.BNL.inter.test',
            'brownies.BNL.reaclib',
            'brownies.LANL', 'brownies.LANL.toACE',
            'brownies.LLNL', 'brownies.LLNL.fetePy',
            'brownies.legacy',
            'brownies.legacy.converting',
            'brownies.legacy.converting.ENDFToGNDS',
            'brownies.legacy.endl',
            'brownies.legacy.endl.structure',
            'brownies.legacy.endl.test',
            'brownies.legacy.toENDL', 'brownies.legacy.toENDL.productData', 'brownies.legacy.toENDL.productData.distributions',
            'brownies.legacy.toENDF6', 'brownies.legacy.toENDF6.PoPs_toENDF6', 'brownies.legacy.toENDF6.PoPs_toENDF6.atomic',
            'brownies.legacy.toENDF6.PoPs_toENDF6.decays', 'brownies.legacy.toENDF6.PoPs_toENDF6.fissionFragmentData',
            'brownies.legacy.toENDF6.outputChannelData', 'brownies.legacy.toENDF6.covariances', 'brownies.legacy.toENDF6.differentialCrossSection',
            'brownies.legacy.toENDF6.productData', 'brownies.legacy.toENDF6.productData.distributions', 'brownies.legacy.toENDF6.reactions',
            'brownies.legacy.toENDF6.reactionData', 'brownies.legacy.toENDF6.reactionData.chargedParticleElastic',
            'brownies.legacy.toENDF6.reactionData.photonScattering', 'brownies.legacy.toENDF6.resonances'],
        package_dir={'brownies': '.'},
        scripts=['legacy/bin/prepro.py'] + glob.glob('bin/*.py'),
        package_data = {
            'brownies.BNL.inter': ['README.txt', '*.json'],
            'brownies.BNL.inter.test': ['*.endf'],
            'brownies.BNL.inter.spectra': ['*.001', 'README.txt', '*.dat', '*.json', '*.endf'],
            'brownies.BNL.reaclib': ['SkyNetReaclib/*'],
            'brownies.BNL.plot_evaluation': ['*.json', '*.DAT'],
            'brownies.LANL': ['build_xsdir/*.py', 'build_xsdir/*.dat', 'dismemberACE/*.py']
        },
        version='0.9.1',
        description='Fudge brownies',
        author='David Brown',
        author_email='dbrown@bnl.gov',
        url='http://www.nndc.bnl.gov',
        download_url='',
        keywords=['ENDF', 'Fudge'],
        classifiers=[
            'Development Status :: ???',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: ???',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Fortran',
            'Programming Language :: Python :: 2',
            'Topic :: Text Processing :: General',
            'Topic :: Software Development :: Interpreters',
            'Topic :: Scientific/Engineering',
        ],
        long_description=''''''
    )


if __name__ == '__main__':
    setup()
