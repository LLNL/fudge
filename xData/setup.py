#!/usr/bin/env python
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import glob

"""
This module is used to build and install the xData module.
"""

def setup():
    from setuptools import setup

    setup(
        name='xData',
        version='1.0.0',
        maintainer='mattoon1@llnl.gov',
        packages=[
            'xData',
            'xData.Documentation',
            'xData.uncertainty',
            'xData.uncertainty.physicalQuantity',
            'xData.interactivePlot'
        ],
        package_dir={'xData': '.'},
        scripts = glob.glob('bin/*.py'),
        install_requires=[
            'numericalFunctions',
            'pqu',
            'numpy>=1.15'
        ],
        description='',
        license=''
    )

if __name__ == '__main__':
    setup()
