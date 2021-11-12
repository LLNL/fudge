#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
from argparse import ArgumentParser

from brownies.legacy.endl import fudgeParameters as fudgeParametersModule

fudgeParametersModule.VerboseMode = 0
from brownies.legacy.endl import endlProject as endlProjectClass

usage = \
"""
Replaces all gamma data (i.e., yo07* and *c55* data) in an ENDL neutron database with the gamma data produced by fetePy.py.
"""

epsilonDefault = 2e-6

parser = ArgumentParser( description = usage )
parser.add_argument( 'fetePath', type = str,
        help = 'Path to fete ENDL data.' )
parser.add_argument( 'fetePyPath', type = str,
        help = 'Path to fetePy ENDL data.' )
parser.add_argument( 'path', type = str,
        help = 'Path to put ENDL data with replaced fetePy ENDL data.' )
parser.add_argument( 'isotope', type = str,
        help = 'name of isotope directory (e.g., "za008016").' )
parser.add_argument( '--epsilon', default = epsilonDefault, type = float, action = 'store',
        help = 'Enable verbose output. Default = %s.' % epsilonDefault )
parser.add_argument( '-v', '--verbose', default = 0, action = 'count',
        help = 'Enable verbose output.' )

args = parser.parse_args( )

def replace( isotope, isotopePy, gammas, gammasPy, ignoreC55 = False ) :

    CSX1s = []
    gammasCs = {}
    for gamma in gammas : 
        if( ignoreC55 and ( gamma.C == 55 ) ) : continue
        C, S, X1 = gamma.C, gamma.S, gamma.X1
        CSX1 = [ C, S, X1 ]
        if( CSX1 not in CSX1s ) : CSX1s.append( CSX1 )
        if( C not in gammasCs )        : gammasCs[C] = {}
        if( S not in gammasCs[C] )     : gammasCs[C][S] = {}
        if( X1 not in gammasCs[C][S] ) : gammasCs[C][S][X1] = {}
        gammasCs[C][S][X1][gamma.I] = gamma
    CSX1s = sorted( CSX1s )

    CSX1sPy = []
    gammasCsPy = {}
    for gamma in gammasPy : 
        if( ignoreC55 and ( gamma.C == 55 ) ) : continue
        C, S, X1 = gamma.C, gamma.S, gamma.X1
        CSX1 = [ C, S, X1 ]
        if( CSX1 not in CSX1sPy ) : CSX1sPy.append( CSX1 )
        if( C not in gammasCsPy )        : gammasCsPy[C] = {}
        if( S not in gammasCsPy[C] )     : gammasCsPy[C][S] = {}
        if( X1 not in gammasCsPy[C][S] ) : gammasCsPy[C][S][X1] = {}
        gammasCsPy[C][S][X1][gamma.I] = gamma
    CSX1sPy = sorted( CSX1sPy )

    Cs = sorted( set( list( gammasCs.keys( ) ) + list( gammasCsPy.keys( ) ) ) )
    for C in Cs :
        if( C not in gammasCs ) :
            print( '    WARNING: fete missing gamma data for C=%s' %  C )
            for S in gammasCsPy[C] :
                for X1 in gammasCsPy[C][S] :
                    print( '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 ) )
            isotopePy.removeFile( C = C, S = S, yo = 7 )
            if( C == 55 ) : isotopePy.removeFile( C = C, S = S, yo = 0 )
            continue
        if( C not in gammasCsPy ) :
            print( '    WARNING: fetePy missing gamma data for C=%s' % C )
            continue
        Ss = sorted( set( list( gammasCs[C].keys( ) ) + list( gammasCsPy[C].keys( ) ) ) )
        for S in sorted( set( list( gammasCs[C].keys( ) ) + list( gammasCsPy[C].keys( ) ) ) ) :
            if( S not in gammasCs[C] ) :
                print( '    WARNING: fete missing gamma data for C=%s S=%s' % ( C, S ) )
                for X1 in gammasCsPy[C][S] :
                    print( '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 ) )
                isotopePy.removeFile( C = C, S = S, yo = 7 )
                if( C == 55 ) : isotopePy.removeFile( C = C, S = S, yo = 0 )
                continue
            if( S not in gammasCsPy[C] ) :
                print( '    WARNING: fetePy missing gamma data for C=%s S=%s' % ( C, S ) )
                continue
            gammasCSs   = gammasCs[C][S]
            gammasCSsPy = gammasCsPy[C][S]
            X1s   = sorted( gammasCSs.keys( ) )
            X1sPy = sorted( gammasCSsPy.keys( ) )
            commonX1s = []
            for X1 in X1s :
                if( X1 in X1sPy ) :
                    X1sPy.remove( X1 )
                    commonX1s.append( [ X1, X1 ] )
            for X1, X1p in commonX1s : X1s.remove( X1 )

            matches = [ [] for i2, X2 in enumerate( X1sPy ) ]
            for i1, X1 in enumerate( X1s ) :
                diffs = []
                for i2, X2 in enumerate( X1sPy ) :
                    diff = abs( X1 - X2 ) / max( abs( X1) , abs( X2 ) )
                    if( diff > args.epsilon ) : continue
                    diffs.append( [ abs( X1 - X2 ), i1, i2, X1, X2 ] )
                if( len( diffs ) == 1 ) : matches[diffs[0][2]].append( diffs[0] )       # This X1 has 1 and only 1 match.
            uniqueMatches = []
            for match in matches :
                if( len( match ) == 1 ) : uniqueMatches.append( match[0] )              # This X2 has 1 and only 1 match.
            for diff, i1, i2, X1, X2 in uniqueMatches :
                commonX1s.append( [ X1, X2 ] )
                X1s.remove( X1 )
                X1sPy.remove( X2 )
            for X1 in X1s   : print( '    WARNING: fetePy missing gamma data for C=%s S=%s X1=%s' % ( C, S, X1 ) )
            for X1 in X1sPy :
                print( '    WARNING: fete missing gamma data for C=%s S=%s X1=%s' % ( C, S, X1 ) )
                print( '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 ) )
                files = isotopePy.findFiles( C = C, S = S, yo = 7 )
                for file in files : print( file )
                raise Exception( 'Get Bret to fix this' )

# unreadData

            if( args.verbose > 1 ) : print('    For C=%2s, S=%s' % ( C, S ) )
            commonX1s.sort( )
            for X1, X2 in commonX1s :
                if( args.verbose > 1 ) : print('        replacing fete X1=%10s with fetePy X1=%10s' % ( X1, X2 ) )
                gamma = gammasCSs[X1][list( gammasCSs[X1].keys( ) )[0]]
                for I in gammasCSsPy[X2].keys( ) :
                    gammaPy = gammasCSsPy[X2][I]
                    gammaPy.setTemperature( gamma.getTemperature( ) )
                    gammaPy.setQ( gamma.getQ( ) )
                    if( S != 3 ) : gammaPy.setX1( gamma.getX1( ) )

            isotope.removeFile( yo = 7, C = C, S = S )
            if( C == 55 ) : isotope.removeFile( I = 0, C = C, S = S )

fete = endlProjectClass( database = args.fetePath, workDir = args.path )
output = os.path.join( fete.workDir, args.isotope )
if( os.path.exists( output ) ) : os.system( 'rm -rf %s' % output )
isotope = fete.readZA( args.isotope )
if( args.verbose > 0 ) : print( 'Fete source = %s' % isotope.source )

isotope.read( yo = 7 )
isotope.read( C = 55, I = 0 )

fetePy = endlProjectClass( database = args.fetePyPath, readOnly = True )
isotopePy = fetePy.readZA( args.isotope )
if( args.verbose > 0 ) : print( 'FetePy source = %s' % isotopePy.source )
if( args.verbose > 0 ) : print( 'Epsilon = %s' % args.epsilon )

isotopePy.read( yo = 7 )
isotopePy.read( C = 55, I = 0 )

C55s = isotope.findDatas( C = 55 )
C55sPy = isotopePy.findDatas( C = 55 )
replace( isotope, isotopePy, C55s, C55sPy )

yo07s = isotope.findDatas( yo = 7 )
yo07sPy = isotopePy.findDatas( yo = 7 )
replace( isotope, isotopePy, yo07s, yo07sPy, ignoreC55 = True )

isotope.read( I = 0 )
Cs = isotope.CList( )
for C in isotopePy.CList( ) :
    if( C in [ 1, 55 ] ) : continue
    if( C not in Cs ) :
        print( '    WARNING: removing C = %d from fetePy as reaction not in ENDL' % C )
        isotopePy.removeFile( C = C )

for file in isotopePy.findFiles( yo = 7 ) + isotopePy.findFiles( C = 55, I = 0 ) :
    isotope.addEndlFile( file, checkHeaders = False )
fete.save( )
