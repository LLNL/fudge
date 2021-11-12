# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
setup.py to support installing statusMessageReporting via pip.

This setup.py was based on the example at
https://stackoverflow.com/questions/33168482/compiling-installing-c-executable-using-pythons-setuptools-setup-py

Run as 'python setup.py install'
"""

import os
import sys
from setuptools import setup
from setuptools.command.install import install
import subprocess

    
def compile_with_make():
    cmd = "make -j"
    subprocess.check_call(cmd, shell=True)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_with_make()
        super().run()


setup(
    name='statusMessageReporting',
    version='1.0.0',
    maintainer='mattoon1@llnl.gov',
    packages=['statusMessageReporting'],
    package_dir={'statusMessageReporting': ''},
    package_data = {'statusMessageReporting': ['Src/*.c', 'Src/*.h']},
    cmdclass={'install': CustomInstall}
)

