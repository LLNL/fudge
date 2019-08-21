# <<BEGIN-copyright>>
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
shutil.copyfile( libs[0], os.path.join( 'lib', 'crossSectionAdjustForHeatedTarget.so' ) )
