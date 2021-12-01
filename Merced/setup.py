# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
setup.py to support installing Merced via pip.

This setup.py was based on the example at
https://stackoverflow.com/questions/33168482/compiling-installing-c-executable-using-pythons-setuptools-setup-py

Run as 'python setup.py install'
"""

from setuptools import setup
from setuptools.command.install import install
import subprocess


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        cmd = f'make -j'
        subprocess.check_call(cmd, shell=True)
        super().run()


setup(
    name='merced',
    version='1.0.0',
    maintainer='mattoon1@llnl.gov',
    packages=['merced'],
    package_dir={'merced': 'bin'},
    data_files=[('bin', ['bin/merced'])],
    cmdclass={'install': CustomInstall}
)

