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

from distutils.core import setup, Extension
import os, glob, shutil

sep = os.path.sep
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'crossSectionAdjustForHeatedTarget.so' )
for lib in libs : os.remove( lib )

sources =  glob.glob( 'Src' + sep + '*.c' ) 
sources += glob.glob( os.path.join( 'Python', 'Src' ) + sep + '*.c' )

crossSectionAdjustForHeatedTarget = Extension( 'crossSectionAdjustForHeatedTarget', sources = sources, include_dirs = [ 'Src' ] )

setup(  name = 'crossSectionAdjustForHeatedTarget', \
        version = '1.0', \
        description = 'This is the crossSectionAdjustForHeatedTarget package', \
        ext_modules = [ crossSectionAdjustForHeatedTarget ] \
)
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'crossSectionAdjustForHeatedTarget.so' )
if len(libs)>0: shutil.copyfile( libs[0], os.path.join( 'lib', 'crossSectionAdjustForHeatedTarget.so' ) )
