# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
from fudge.core.utilities import fudgeFileMisc
from brownies.legacy.endl import bdfls, endl2, endl_Z


def processTDF_Reaction( target, C, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None, outputFile = 'tdfgen.out', workDir = None,
        bdflsFile = None ) :

    if( bdflsFile is None ) : bdflsFile = bdfls.getDefaultBdfls()
    AMUToMeV = bdflsFile.constant( 4 )
    xsec = target.findData( C = C, I = 0, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )

    residualZA, yos, Q = endl2.residualZA_yos_Q(target.yi, target.ZA, C, bdflsFile = bdflsFile)
    projectileMass = bdflsFile.mass( target.yi )
    targetMass = bdflsFile.mass( target.ZA )

    yi = endl2.ZAToYo(target.yi)
    yiZA = endl2.yoToZA(yi)
    yiZ, yiA = endl2.ZandAFromZA(yiZA)

    ZA = target.ZA
#    ZAZA = endl2.yoToZA( ZA )
    ZAZA = ZA
    ZAZ, ZAA = endl2.ZandAFromZA(ZAZA)

    if( projectileMass > targetMass ) :
        reaction = '%s%d__%s%d_' % (endl_Z.endl_ZSymbol(yiZ), yiA, endl_Z.endl_ZSymbol(ZAZ), ZAA)
    else :
        reaction = '%s%d__%s%d_' % (endl_Z.endl_ZSymbol(ZAZ), ZAA, endl_Z.endl_ZSymbol(yiZ), yiA)

    outGoing = []
    print( yos )
    for yo in yos : 
        if( yo < 1000 ) :
            outGoing.append([endl2.yoToZA(yo)])
        else :
            outGoing.append( [ yo ] )
    outGoing.append( [ residualZA ] )

    for i in outGoing :
        iZA = i[0]
        i.insert( 0, bdflsFile.mass( iZA ) )
        Z, A = endl2.ZandAFromZA(iZA)
        i.append( Z )
        i.append( A )
    outGoing.sort( )
    s = ''
    for mass, iZA, Z, A in outGoing :
        reaction += '%s%s%d' % (s, endl_Z.endl_ZSymbol(Z), A)
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
    print( inputName )
    os.system( './tdfgen -i %s -o %s' % ( inputFile.getName( ), outputFile ) )
