# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from setuptools import setup, Extension
import os, glob, shutil

libs = glob.glob(os.path.join('build', 'lib*', 'crossSectionAdjustForHeatedTarget*'))
for lib in libs: os.remove(lib)

sources = glob.glob(os.path.join('Src', '*.c'))
sources += glob.glob(os.path.join('Python', 'Src', '*.c'))

crossSectionAdjustForHeatedTarget = Extension('crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget', sources=sources, include_dirs=['Src'])

setup(
    name='crossSectionAdjustForHeatedTarget',
    version='1.0.0',
    maintainer='mattoon1@llnl.gov',
    packages=['crossSectionAdjustForHeatedTarget'],
    package_dir={'crossSectionAdjustForHeatedTarget': '.'},
    ext_modules=[crossSectionAdjustForHeatedTarget],
    description='',
    license=''
)
