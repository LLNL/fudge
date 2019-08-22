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

import os
from fudge import fudgeParameters
import endlZA
import endlFile

EOD = '                                                                       1\n'
tmpFileName = None

def endlSplit( file, outDir, endlretHeaders = True, verbose = 0, format = None ) :
    """The function undoes what endlret does. That is, it takes an ascii file with a collect of
    isotopes in ENDL ascii format, and splits the isotopes up into separate 'za' directories
    with the appropriate 'yo' files. The za directories are placed in outDir."""

    def getData( fIn, l, ZA ) :

        data = [ l ]
        while( l != EOD ) :
            l = fIn.readline( )
            data.append( l )
        l = fIn.readline( )
        return( l, data )

    def splitZA( fIn, l, yi, workDir, format ) :

        ZA = int( l[:6] )
        if( verbose != 0 ) : print 'splitting ZA %5d' % ZA
        ZA_ = ZA
        target = endlZA.endlZA( ZA, yi, workDir = workDir )
        while( ZA == ZA_ ) :
            l, data = getData( fIn, l, ZA )
            yo = int( data[0][9:12] )
            C = int( data[1][:2] )
            I = int( data[1][2:5] )
            S = int( data[1][5:8] )
            tmpFile = os.path.join( outDir, 'yo%2.2dc%2.2di%3.3ds%3.3d' % ( yo, C, I, S ) )
            fTmp = open( tmpFile, 'w' )
            fTmp.write( ''.join( data ) )
            fTmp.close( )
            endlData = endlFile.endlFile( 0, os.path.basename( tmpFile ), tmpFile, ZA, yi )
            endlData.read( )
            file = target.findFile( yo = yo, C = C, I = I, S = S )
            if( file == None ) : file = target.addFile( yo = yo, C = C, I = I, S = S, printWarnings = False )
            if( format != None ) : endlData[0].setFormat( format )
            file.addEndlData( endlData[0] )
            os.remove( tmpFile )
            if( l == '' ) : break
            ZA_ = int( l[:6] )
        target.save( )
        return( l )

    fudgeParameters.VerboseMode = verbose > 0
    fIn = open( file )
    if( endlretHeaders ) :
        l = fIn.readline( )
        if( 'date of processing' not in l ) : raise Exception( 'Invalid first line for endlret type file %s' % file )
        l = fIn.readline( )
        if( 'date of last physics change' not in l ) : raise Exception( 'Invalid second line for endlret type file %s' % file )

    l = fIn.readline( )
    yi = int( l[6:9] )
    workDir = os.path.join( outDir, 'yi%2.2d' % yi )
    tmpFileName = os.path.join( workDir, 'tmpFile' )
    while( l != '' ) : l = splitZA( fIn, l, yi, workDir, format = format )
    fIn.close( )
