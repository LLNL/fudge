# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

fIn = open( 'check3.out' )
lines = fIn.readlines( )
fIn.close( )

splitString = '    At line '

fOut = open( 'check3.out', 'w' )
for line in lines :
    if( splitString in line ) :
        start, end = line.split( splitString )
        end = ' '.join( end.split( )[1:] )
        line = "%s%snn %s\n" % ( start, splitString, end )
    fOut.write( line )
fOut.close( )
