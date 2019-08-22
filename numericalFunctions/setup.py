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

import os, glob, shutil
from distutils.core import setup, Extension

#
# Stuff to build listOfDoubles.so.
#
sep = os.path.sep
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'listOfDoubles.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

listOfDoubles_C = Extension( 'listOfDoubles_C', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 'ptwX' + sep + 'Src', 'Src' ] )

setup( name = 'listOfDoubles_C', 
    version = '1.0', \
    description = 'This module contains the listOfDoubles_C class and support routines.', \
    ext_modules = [ listOfDoubles_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'listOfDoubles_C.so' )
if len(libs)>0: shutil.copyfile( libs[0], os.path.join( 'lib', 'listOfDoubles_C.so' ) )

#
# Stuff to build pointwiseXY_C.so.
#
sep = os.path.sep
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

pointwiseXY_C = Extension( 'pointwiseXY_C', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 'ptwX' + sep + 'Src', 'ptwXY' + sep + 'Src', 'ptwXY' + sep + 'Python' + sep + 'Src' ] )

setup( name = 'pointwiseXY_C', 
    version = '1.0', \
    description = 'This module contains the pointwiseXY_C class and support routines.', \
    ext_modules = [ pointwiseXY_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
if len(libs)>0: shutil.copyfile( libs[0], os.path.join( 'lib', 'pointwiseXY_C.so' ) )

#
# Stuff to build Legendre.so.
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'Legendre.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_Legendre' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_Legendre' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

nf_Legendre_C = Extension( 'Legendre', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src',                   'ptwX' + sep + 'Src',       'ptwXY' + sep + 'Src', 
                     'ptwXY' + sep + 'Python' + sep + 'Src', 'nf_Legendre' + sep + 'Src' ] )

setup( name = 'Legendre', 
    version = '1.0', \
    description = 'This module contains the Legendre class and support routines.', \
    ext_modules = [ nf_Legendre_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'Legendre.so' )
if len(libs)>0: shutil.copyfile( libs[0], os.path.join( 'lib', 'Legendre.so' ) )

#
# Stuff to build specialFunctions.so
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'specialFunctions.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_specialFunctions' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_specialFunctions' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

specialFunctions = Extension( 'specialFunctions', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 'nf_specialFunctions' + sep + 'Python' + sep + 'Src', 'nf_specialFunctions' + sep + 'Src' ] )

setup( name = 'specialFunctions', 
    version = '1.0', \
    description = 'This module contains some special math functions not in the python math module.', \
    ext_modules = [ specialFunctions ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'specialFunctions.so' )
if len(libs)>0: shutil.copyfile( libs[0], os.path.join( 'lib', 'specialFunctions.so' ) )
