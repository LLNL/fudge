from __future__ import print_function
import glob, os, shutil, filecmp

files = sorted( glob.glob( 't*.py' ) )

if( os.path.exists( 'Out' ) ) : shutil.rmtree( 'Out' )
os.mkdir( 'Out' )

for file in files :
    base = file[:-3]
    status = os.system( 'python %s > Out/%s.out' % ( file, base ) )
    if( status ) : print( '=========== %s ===========' % file )

outs = sorted( glob.glob( 'Out/t*.out' ) )
for out in outs :
    file = os.path.basename( out )
    if( not( filecmp.cmp( os.path.join( 'Out.checked', file ), out ) ) ) : print( 'ERROR: %s' % out )
