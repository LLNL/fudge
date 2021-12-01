# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
from fudge import fudgeParameters
from brownies.legacy.endl import endlZA
from brownies.legacy.endl import endlFile

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
        if( verbose != 0 ) : print('splitting ZA %5d' % ZA)
        ZA_ = ZA
        target = endlZA( ZA, yi, workDir = workDir )
        while( ZA == ZA_ ) :
            l, data = getData( fIn, l, ZA )
            yo = int( data[0][9:12] )
            C = int( data[1][:2] )
            I = int( data[1][2:5] )
            S = int( data[1][5:8] )
            if( C == 93 ) :             # Happens in EPICS.
                if( I == 941 ) :
                    C = 71
                elif( I == 942 ) :
                    C = 72
                elif( I == 942 ) :
                    C = 72
                elif( I == 942 ) :
                    C = 72
                else :
                    pass
                data[1] = ( '%2s' % C ) + data[1][2:]
            if( C != 93 ) :
                tmpFile = os.path.join( outDir, 'yo%2.2dc%2.2di%3.3ds%3.3d' % ( yo, C, I, S ) )
                fTmp = open( tmpFile, 'w' )
                fTmp.write( ''.join( data ) )
                fTmp.close( )
                endlData = endlFile.endlFile(0, os.path.basename(tmpFile), tmpFile, ZA, yi)
                endlData.read( )
                file = target.findFile( yo = yo, C = C, I = I, S = S )
                if( file is None ) : file = target.addFile( yo = yo, C = C, I = I, S = S, printWarnings = False )
                if( format is not None ) : endlData[0].setFormat( format )
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
