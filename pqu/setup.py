#!/usr/bin/env python
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


def setup():
    from setuptools import setup

    setup(
        name='pqu',
        version='1.1.0',
        maintainer='mattoon1@llnl.gov',
        packages=['pqu'],
        package_dir={'pqu': '.'},
        description='',
        license=''
    )

if __name__ == '__main__':
    setup()
