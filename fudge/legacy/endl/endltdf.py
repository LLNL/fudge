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
import sys
from fudge.core.utilities import fudgeFileMisc
import bdfls
import endl2
import endl_Z

def processTDF_Reaction( target, C, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None, outputFile = 'tdfgen.out', workDir = None,
        bdflsFile = None ) :

    if( bdflsFile == None ) : bdflsFile = bdfls.getDefaultBdfls( )
    AMUToMeV = bdflsFile.constant( 4 )
    xsec = target.findData( C = C, I = 0, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )

    residualZA, yos, Q = endl2.residualZA_yos_Q( target.yi, target.ZA, C, bdflsFile = bdflsFile )
    projectileMass = bdflsFile.mass( target.yi )
    targetMass = bdflsFile.mass( target.ZA )

    yi = endl2.ZAToYo( target.yi )
    yiZA = endl2.yoToZA( yi )
    yiZ, yiA = endl2.ZandAFromZA( yiZA )

    ZA = target.ZA
#    ZAZA = endl2.yoToZA( ZA )
    ZAZA = ZA
    ZAZ, ZAA = endl2.ZandAFromZA( ZAZA )

    if( projectileMass > targetMass ) :
        reaction = '%s%d__%s%d_' % ( endl_Z.endl_ZSymbol( yiZ ), yiA, endl_Z.endl_ZSymbol( ZAZ ), ZAA )
    else :
        reaction = '%s%d__%s%d_' % ( endl_Z.endl_ZSymbol( ZAZ ), ZAA, endl_Z.endl_ZSymbol( yiZ ), yiA )

    outGoing = []
    print yos
    for yo in yos : 
        if( yo < 1000 ) :
            outGoing.append( [ endl2.yoToZA( yo ) ] )
        else :
            outGoing.append( [ yo ] )
    outGoing.append( [ residualZA ] )

    for i in outGoing :
        iZA = i[0]
        i.insert( 0, bdflsFile.mass( iZA ) )
        Z, A = endl2.ZandAFromZA( iZA )
        i.append( Z )
        i.append( A )
    outGoing.sort( )
    s = ''
    for mass, iZA, Z, A in outGoing :
        reaction += '%s%s%d' % ( s, endl_Z.endl_ZSymbol( Z ), A )
        s = '__'

    outputStr  = [ '## Fudge generated data for tdfgen version:0.9.9' ]
    outputStr.append( '## Data generated from:fudge' )
    outputStr.append( '' )
    outputStr.append( '## Reaction:%s' % reaction )
    outputStr.append( '' )
    outputStr.append( '# Masses of particles in MeV.' )
    outputStr.append( '## Mass of projectile:%.12e' % ( projectileMass * AMUToMeV ) )
    outputStr.append( '## Mass of target:%.12e' % ( targetMass * AMUToMeV ) )

    outputStr.append( '' )
    outputStr.append( '## Number of final particles:%d' % len( outGoing ) )
    for mass, iZA, Z, A in outGoing : outputStr.append( '## %.12e' % ( mass * AMUToMeV ) )

    outputStr.append( '' )
    outputStr.append( '## Lab of CM frame:Lab' )
    outputStr.append( '## Number of data points:%d' % len( xsec ) )
    outputStr.append( '# E(MeV)    Sigma( barn )' )
    outputStr.append( '#------------------------' )
    outputStr.append( xsec.toString( ) )

    outputStr = '\n'.join( outputStr )

    inputFile = fudgeFileMisc.fudgeTempFile( dir = workDir )
    inputFile.write( outputStr )

    inputName = inputFile.getName( )
    print inputName
    os.system( './tdfgen -i %s -o %s' % ( inputFile.getName( ), outputFile ) )
