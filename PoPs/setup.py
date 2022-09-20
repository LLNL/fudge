# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


def setup():
    from setuptools import setup
    setup(
        name='PoPs',
        version='1.0.0',
        maintainer='mattoon1@llnl.gov',
        packages=[
            'PoPs',
            'PoPs.atomic',
            'PoPs.decays',
            'PoPs.families',
            'PoPs.fissionFragmentData',
            'PoPs.chemicalElements',
            'PoPs.quantities',
        ],
        package_dir={'PoPs': '.'},
        install_requires=['pqu', 'xData',],
        description='',
        license=''
    )


if __name__ == '__main__':
    setup()
