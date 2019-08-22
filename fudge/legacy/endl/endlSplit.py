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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
