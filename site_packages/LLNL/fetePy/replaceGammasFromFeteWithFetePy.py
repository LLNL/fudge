#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import os
from argparse import ArgumentParser

from fudge.core.utilities import brb

from fudge import fudgeParameters as fudgeParametersModule
fudgeParametersModule.VerboseMode = 0
from fudge.legacy.endl import endlProject as endlProjectClass

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

    Cs = sorted( set( gammasCs.keys( ) + gammasCsPy.keys( ) ) )
    for C in Cs :
        if( C not in gammasCs ) :
            print '    WARNING: fete missing gamma data for C=%s' %  C
            for S in gammasCsPy[C] :
                for X1 in gammasCsPy[C][S] :
                    print '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 )
            isotopePy.removeFile( C = C, S = S, yo = 7 )
            if( C == 55 ) : isotopePy.removeFile( C = C, S = S, yo = 0 )
            continue
        if( C not in gammasCsPy ) :
            print '    WARNING: fetePy missing gamma data for C=%s' % C
            continue
        Ss = sorted( set( gammasCs[C].keys( ) + gammasCsPy[C].keys( ) ) )
        for S in sorted( set( gammasCs[C].keys( ) + gammasCsPy[C].keys( ) ) ) :
            if( S not in gammasCs[C] ) :
                print '    WARNING: fete missing gamma data for C=%s S=%s' % ( C, S )
                for X1 in gammasCsPy[C][S] :
                    print '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 )
                isotopePy.removeFile( C = C, S = S, yo = 7 )
                if( C == 55 ) : isotopePy.removeFile( C = C, S = S, yo = 0 )
                continue
            if( S not in gammasCsPy[C] ) :
                print '    WARNING: fetePy missing gamma data for C=%s S=%s' % ( C, S )
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
            for X1 in X1s   : print '    WARNING: fetePy missing gamma data for C=%s S=%s X1=%s' % ( C, S, X1 )
            for X1 in X1sPy :
                print '    WARNING: fete missing gamma data for C=%s S=%s X1=%s' % ( C, S, X1 )
                print '        WARNING: fetePy data for C=%s, S=%s and X1=%s not copied' % ( C, S, X1 )
                files = isotopePy.findFiles( C = C, S = S, yo = 7 )
                for file in files : print file
                raise Exception( 'Get Bret to fix this' )

# unreadData

            if( args.verbose > 1 ) : print'    For C=%2s, S=%s' % ( C, S )
            commonX1s.sort( )
            for X1, X2 in commonX1s :
                if( args.verbose > 1 ) : print'        replacing fete X1=%10s with fetePy X1=%10s' % ( X1, X2 )
                gamma = gammasCSs[X1][gammasCSs[X1].keys( )[0]]
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
if( args.verbose > 0 ) : print 'Fete source = %s' % isotope.source

isotope.read( yo = 7 )
isotope.read( C = 55, I = 0 )

fetePy = endlProjectClass( database = args.fetePyPath, readOnly = True )
isotopePy = fetePy.readZA( args.isotope )
if( args.verbose > 0 ) : print 'FetePy source = %s' % isotopePy.source
if( args.verbose > 0 ) : print 'Epsilon = %s' % args.epsilon

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
        print '    WARNING: removing C = %d from fetePy as reaction not in ENDL' % C
        isotopePy.removeFile( C = C )

for file in isotopePy.findFiles( yo = 7 ) + isotopePy.findFiles( C = 55, I = 0 ) :
    isotope.addEndlFile( file, checkHeaders = False )
fete.save( )
