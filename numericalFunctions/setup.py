# <<BEGIN-copyright>>
# <<END-copyright>>

import os, glob, shutil
from distutils.core import setup, Extension

sep = os.path.sep
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

pointwiseXY_C = Extension( 'pointwiseXY_C', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 'ptwX' + sep + 'Src', 'ptwXY' + sep + 'Src',
                     'ptwXY' + sep + 'Python' + sep + 'Src' ] )

setup( name = 'pointwiseXY_C', 
    version = '1.0', \
    description = 'This module contains the pointwiseXY_C class and support routines.', \
    ext_modules = [ pointwiseXY_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
shutil.copyfile( libs[0], os.path.join( 'lib', 'pointwiseXY_C.so' ) )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'nf_Legendre_C.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'ptwXY' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_Legendre' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_Legendre' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

nf_Legendre_C = Extension( 'nf_Legendre_C', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 
                     'ptwX' + sep + 'Src', 
                     'ptwXY' + sep + 'Src', 
                     'ptwXY' + sep + 'Python' + sep + 'Src',
                     'nf_Legendre' + sep + 'Src' ] )

setup( name = 'nf_Legendre_C', 
    version = '1.0', \
    description = 'This module contains the nf_Legendre_C class and support routines.', \
    ext_modules = [ nf_Legendre_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'nf_Legendre_C.so' )
shutil.copyfile( libs[0], os.path.join( 'lib', 'nf_Legendre_C.so' ) )






libs = glob.glob( 'build' + sep + 'lib*' + sep + 'nf_specialFunctions_C.so' )
for lib in libs : os.remove( lib )

sources =   glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_specialFunctions' + sep + 'Src' + sep + '*.c' ) + \
            glob.glob( 'nf_specialFunctions' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

nf_specialFunctions_C = Extension( 'nf_specialFunctions_C', \
    sources = sources, \
    include_dirs = [ 'ptwC' + sep + 'Src', 
                     'nf_specialFunctions' + sep + 'Python' + sep + 'Src',
                     'nf_specialFunctions' + sep + 'Src' ] )

setup( name = 'nf_specialFunctions_C', 
    version = '1.0', \
    description = 'This module contains some special math functions not in the python math module.', \
    ext_modules = [ nf_specialFunctions_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'nf_specialFunctions_C.so' )
shutil.copyfile( libs[0], os.path.join( 'lib', 'nf_specialFunctions_C.so' ) )
