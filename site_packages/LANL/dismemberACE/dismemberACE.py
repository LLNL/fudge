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

from __future__ import print_function
import os
import sys
from argparse import ArgumentParser

outputDefault = 'dismemberACE.out'
cites = []

description = """
This module dismembers the MT data in an ACE file into separate sub-directories, one for each MT.
"""


parser = ArgumentParser( description = description )
parser.add_argument( '-a', '--addresses', action = 'store_true',                help = 'Add addresses (locators) to output (for debbuging).' )
parser.add_argument( '-f', '--file', action = 'store_true',                     help = 'If present, the ACE file name is printed at startup.' )
parser.add_argument( '-k', '--keep', action = 'store_true',                     help = 'The contents of the output directory is removed unless this option is present.' )
parser.add_argument( '-o', '--output', type = str, default = outputDefault,     help = 'Name of directory to put the MT sub-directories.' )
parser.add_argument( '-q', '--quiet', action = 'store_true',                    help = 'Quiet mode. The citation data are not printed if all data have been cited.' )
parser.add_argument( '-O', '--Official', action = 'store_true',                 help = 'Only print information about official cites.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,          help = 'Verbose mode.' )
parser.add_argument( 'ace', type = str,                                         help = 'ACE file name to dismember.' )

args = parser.parse_args( )

__doc__ = description

def openFile( MT, name, subDir = None, MTP = None ) :

    if( MT > 0 ) :
        path = os.path.join( args.output, "MTs", "%.3d" % MT )
    else :
        path = os.path.join( args.output )

    if( subDir is not None ) : path = os.path.join( path, subDir )
    if( MTP is not None ) : path = os.path.join( path, str( MTP ) )
    if( not( os.path.exists( path ) ) ) : os.makedirs( path )
    path = os.path.join( path, name )
    return( open( path, 'w' ) )        

def lineColumn1BaseForXSS( offset ) :

    offset -= 1
    line = offset // 4 + 13
    column = offset % 4 + 1
    return( line, column  )

def get8Integers( line ) :

    values = []
    for i1 in range( 8 ) :
        values.append( int( line[:9] ) )
        line = line[9:]

    return( values )

class listBase1 :

        def __init__( self, values ) : self.values = values
        def __len__( self ) : return( len( self.values ) )
        def __getitem__( self, index ) : return( self.values[index-1] )

def get4Floats( line, number, lineNumber ) :

    initialLine = line
    count = min( 4, number )
    number -= count
    values = []
    for i1 in range( count ) :
        if( len( line ) == 0 ) : break
        try :
            values.append( float( line[:20] ) )
        except :
            print( lineNumber )
            print( initialLine[:-1] )
            print( line )
            raise
        line = line[20:]

    return( number, values )

def getData( label, XSS, start, numberOfPoints, offset = 0, cite = True ) :

    global cites

    start += offset - 1
    if( cite ) :
        total = sum( cites[start:start+numberOfPoints] )
        if( total > 0 ) : print( '    citings already in range %9d to %9d of length %9d: %s' % ( start, start + numberOfPoints, numberOfPoints, label ) )
        for i1 in range( start, start + numberOfPoints ) : cites[i1] += 1
    return( XSS.values[start:start+numberOfPoints] )

def toIntegers( values ) :

    return( [ int( value ) for value in values ] )

def outputXYs1d( MT, name, xs, ys, offset = 0, length = -1, subDir = None, MTP = None ) :

    if( length == -1 ) : length = len( xs )
    xs = xs[offset:offset+length]

    fOut = openFile( MT, name, subDir = subDir, MTP = MTP )
    for i1, x in enumerate( xs ) : fOut.write( "%20.12e %20.12e\n" % ( x, ys[i1] ) )
    fOut.close( )

def interpolationData( fOut, offset, XSS, distType ) :

    NR   = toIntegers( getData( 'NR (%s)' % distType,  XSS, offset, 1 ) )[0]
    if( NR != 0 ) :
        NBT = toIntegers( getData( 'NBT (%s)\n' % distType, XSS, offset + 1,      NR ) )
        INT = toIntegers( getData( 'INT (%s)\n' % distType, XSS, offset + 1 + NR, NR ) )
        for i1 in range( NR ) : fOut.write( '### NBT = %3d   INT = %2d\n' % ( NBT[i1], INT[i1] ) )

    return( 1 + 2 * NR )

def nu_bar( NU, XSS ) :

    if( NU == 0 ) : return

    KNU = toIntegers( getData( 'KNU', XSS, NU, 1 ) )[0]

    if( KNU > 0 ) :
        nu_bar2( 'prompt_or_total_nubar', NU, XSS, LNU = KNU )
    else :
        nu_bar2( 'prompt_nubar', NU + 1, XSS )
        nu_bar2( 'total_nubar', NU + abs( KNU ) + 1, XSS )

def nu_bar2( label, NU, XSS, LNU = 0 ) :

    fOut = openFile( 18, label )

    if( LNU == 0 ) : LNU = toIntegers( getData( 'LNU (%s)' % label,  XSS, NU, 1 ) )[0]
    fOut.write( '### LNU = %d\n' % LNU )
    offset = NU + 1

    if( LNU == 1 ) :
        NC = toIntegers( getData( 'NC (%s)' % label, XSS, offset, 1 ) )[0]
        fOut.write( '### NC = %d\n' % NC )
        Cs =  getData( 'Cs (%s)' % label, XSS, offset + 1,      NC )
        for C in Cs : fOut.write( '    %20.12e\n' % C )
    elif( LNU == 2 ) :
        offset += interpolationData( fOut, offset, XSS, 'energy' )
        NE = toIntegers( getData( 'NE (%s)' % label, XSS, offset, 1 ) )[0]
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (%s)' % label, XSS, offset + 1,      NE )
        nu_bar   = getData( 'nu_bar (%s)' % label,   XSS, offset + 1 + NE, NE )
        for i1 in range( NE ) : fOut.write( '    %20.12e %20.12e\n' % ( energies[i1], nu_bar[i1] ) )
    else :
        raise Exception( 'Invalid LNU = %d' % LNU )

    fOut.close( )

def delayedNeutronData( NXS, JXS, XSS ) :

    NPCR = NXS[8]
    if( NPCR == 0 ) : return

    subDir = 'delayedNeutron'
    BDD    = JXS[25]
    LDNEDL = JXS[26]
    DNEDL = JXS[27]
    energyLocators = toIntegers( getData( 'energyLocators (%s)' % subDir, XSS, LDNEDL, NPCR ) )
    offset = BDD
    for i1 in range( NPCR ) :
        label = 'DEC[%d]' % i1
        fOut = openFile( 18, 'probability_%d' % i1, subDir = subDir )

        DEC = getData( label, XSS, offset, 1 )[0]
        fOut.write( '### fission delayed neutron %20.12e\n' % DEC )
        offset += 1

        offset += interpolationData( fOut, offset, XSS, 'energy' )

        NE = toIntegers( getData( 'NE (%s)' % label, XSS, offset, 1 ) )[0]
        offset += 1
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (%s)' % label, XSS, offset,      NE )
        Ps       = getData( 'nu_bar (%s)' % label,   XSS, offset + NE, NE )
        for i2 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i2], Ps[i2] ) )
        offset += 2 * NE

        dismemberEnergy( 0, 18, 1, DNEDL, energyLocators[i1], XSS, subDir = subDir, MTP = i1 )

    fOut.close( )

def table_12( MT, AND, locator, XSS, subDir ) :
    """MCNP angular data."""

    fOut = openFile( MT, 'angular', subDir )

    NE = toIntegers( getData( 'NE (angular)', XSS, locator, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    energies = getData( 'energies (angular)', XSS, locator + 1, NE )
    LCs = toIntegers( getData( 'energies (angular)', XSS, locator + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( 'energy = %20.12e\n' % energies[i1] )
        locator2 = AND + abs( LCs[i1] ) - 1
        if( args.addresses ) : fOut.write( '    locator = %d (%d)\n' % ( LCs[i1], locator2 ) )
        if( LCs[i1] < 0 ) :
            interpolation, NP = toIntegers( getData( 'interpolation, NP (angular)', XSS, locator2, 2 ) )
            locator2 += 2
            interpolationStr = { 1 : 'flat', 2 : 'lin-lin' }
            fOut.write( '    interpolation = %s (%s)\n' % ( interpolation, interpolationStr[interpolation] ) )
            fOut.write( '    NP = %s\n' % NP )

            mus = getData( 'mus', XSS, locator2, NP )
            locator2 += NP
            pdf = getData( 'mus', XSS, locator2, NP )
            locator2 += NP
            cdf = getData( 'mus', XSS, locator2, NP )
            for i2 in range( NP ) : fOut.write( '        %20.12e %20.12e %20.12e\n' % ( mus[i2], pdf[i2], cdf[i2] ) )

        elif( LCs[i1] == 0 ) :                  # Isotopic
            fOut.write( '        isotropic' )

        else :                                  # 32 equal probable bins.
            epbs = getData( 'epbs (angular)', XSS, AND + LCs[i1] - 1, 33 )
            for epb in epbs : fOut.write( '        %20.12e\n' % epb )

    fOut.close( )

def dismemberAngular( index, MT, Type, LAND, AND, XSS, subDir, MTP = None ) :

    frame = 'lab'
    if( Type < 0 ) : frame = 'center of mass'
    locator = 0
    if( LAND[index] > 0 ) : locator = AND + LAND[index] - 1

    fOut = openFile( MT, 'info', subDir = subDir, MTP = MTP )
    fOut.write( 'Type = %d\n' % Type )          # abs( Type ) is multiplicity. If 19 fission. If > 100 is non-fission energy dependent.
    fOut.write( '    multiplicity = %d\n' % abs( Type ) )
    fOut.write( '    frame = %s\n' % frame )
    if( args.addresses ) : fOut.write( 'LAND[%d] = %s\n' % ( index, LAND[index] ) )
    if( args.addresses ) : fOut.write( 'locator (1 based) = %s\n' % locator )
    fOut.close( )

    if( LAND[index] < 0 ) :                         # Kalback/Mann data given in the energy section.
        pass
    elif( LAND[index] == 0 ) :
        fOut = openFile( MT, 'angular', subDir = subDir, MTP = MTP )
        fOut.write( '### isotropic\n' )
        fOut.close( )
    else :
        table_12( MT, AND, locator, XSS, 'neutron' )

def dismemberEnergy( level, MT, Type, DLW, offset, XSS, subDir, MTP = None ) :

    fOut = openFile( MT, 'energy_%d' % level, subDir = subDir, MTP = MTP )

    if( args.verbose > 2 ) : print( '            DLW = %d  offset = %d  start = %d: %s' % ( DLW, offset, DLW + offset - 1, subDir ) )
    offset += DLW - 1

    fOut.write( 'Type = %d\n' % Type )
    fOut.write( 'DLW = %d\n' % DLW )

    LNW = toIntegers( getData( 'LNW (energy)',   XSS, offset,     1 ) )[0]
    fOut.write( 'LNW = %d\n' % LNW )
    LAW  = toIntegers( getData( 'LAW (energy)',  XSS, offset + 1, 1 ) )[0]
    fOut.write( 'LAW = %d\n' % LAW )
    IDAT = toIntegers( getData( 'IDAT (energy)', XSS, offset + 2, 1 ) )[0]
    fOut.write( 'IDAT = %d\n' % IDAT )
    offset += 3

    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( 'NE = %d\n' % NE )
    energies =  getData( 'energies (energy)', XSS, offset + 1,      NE )
    Ps       =  getData( 'P(E) (energy)',     XSS, offset + 1 + NE, NE )
    for i1 in range( NE ) : fOut.write( 'energy = %20.12e %7d\n' % ( energies[i1], Ps[i1] ) )

    if( args.verbose > 2 ) : print( '            %sLAW = %d' % ( level * '    ', LAW ) )
# Missing LAWs 1, 22 (UK), 24 (UK)
    if( LAW == 2 ) :
        LAW2( fOut, DLW, IDAT, XSS )
    elif( LAW == 3 ) :
        LAW3( fOut, DLW, IDAT, XSS )
    elif( LAW == 4 ) :
        LAW4( fOut, DLW, IDAT, XSS )
    elif( LAW == 5 ) :
        LAW5( fOut, DLW, IDAT, XSS )
    elif( LAW == 7 ) :
        LAW7( fOut, DLW, IDAT, XSS )
    elif( LAW == 9 ) :
        LAW9( fOut, DLW, IDAT, XSS )
    elif( LAW == 11 ) :
        LAW11( fOut, DLW, IDAT, XSS )
    elif( LAW == 44 ) :
        LAW44_KalbachMann( fOut, DLW, IDAT, XSS )
    elif( LAW == 61 ) :
        LAW61( fOut, DLW, IDAT, XSS )
    elif( LAW == 66 ) :
        LAW66( fOut, DLW, IDAT, XSS )
    elif( LAW == 67 ) :
        LAW67( fOut, DLW, IDAT, XSS )
    else :
        print( "Unsupported LAW = %s" % LAW )
    fOut.close( )

    if( LNW != 0 ) :
        dismemberEnergy( level + 1, MT, Type, DLW, LNW, XSS, subDir, MTP = MTP )

def LAW2( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW2\n' )

    offset += DLW - 1
    LDAT1, LDAT2 = getData( 'LAW2 (energy)', XSS, offset, 2 )
    fOut.write( 'LDAT[1] = %d\n' % int( LDAT1 ) )
    fOut.write( 'LDAT[2] = %20.12e\n' % LDAT2 )

def LAW3( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW3\n' )

    offset += DLW - 1
    LDAT1, LDAT2 = getData( 'LAW3 (energy)', XSS, offset, 2 )
    fOut.write( 'LDAT[1] = %20.12e\n' % LDAT1 )
    fOut.write( 'LDAT[2] = %20.12e\n' % LDAT2 )

def LAW4( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW4\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )
    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '### INTT = %s\n' % INTT )
        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '### NP = %s\n' % NP )
        EPs = getData( 'EPS (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e\n' % ( EPs[i2], pdf[i2], cdf[i2] ) )

def LAW5( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW5\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )
    offset += 2 * NE

    NET = toIntegers( getData( 'NET (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NET = %d\n' % NET )
    offset += 1

    energies = getData( 'energies (energy)', XSS, offset,     NET )
    pdf = getData( 'pdf (energy)', XSS, offset + NET,     NET )
    cdf = getData( 'cdf (energy)', XSS, offset + 2 * NET, NET )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e %20.12e\n' % ( energy, pdf[i1], cdf[i1] ) )

def LAW7( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW7\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NE, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )

def LAW9( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW9\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NE, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )

def LAW11( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW11\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NEa = toIntegers( getData( 'NEa (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NEa = %d\n' % NEa )
    offset += 1
    energies_a = getData( 'energies_a (energy)', XSS, offset,       NEa )
    As =         getData( 'As (energy)',         XSS, offset + NEa, NEa )
    offset += 2 * NEa

    offset += interpolationData( fOut, offset, XSS, 'energy' )
    NEb = toIntegers( getData( 'NEb (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NEb = %d\n' % NEb )
    offset += 1
    energies_b = getData( 'energies_b (energy)', XSS, offset,       NEb )
    Bs =         getData( 'Bs (energy)',         XSS, offset + NEb, NEb )

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NEb, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    fOut.write( '### as\n' )
    for i1, energy in enumerate( energies_a ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, As[i1] ) )

    fOut.write( '\n\n### bs\n' )
    for i1, energy in enumerate( energies_b ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Bs[i1] ) )

def LAW44_KalbachMann( fOut, DLW, offset, XSS ) :

    if( args.verbose > 2 ) : print( '            offset = %d  start = %d' % ( offset, DLW + offset - 1 ) )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '###    INTT = %s\n' % INTT )

        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '###    NP = %s\n' % NP )
        EPs = getData( 'EPS (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        Rs  = getData( 'Rs (energy)',  XSS, offset + 2 + 3 * NP, NP )
        As  = getData( 'As (energy)',  XSS, offset + 2 + 4 * NP, NP )
        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e %20.12e %20.12e\n' % ( EPs[i2], pdf[i2], cdf[i2], Rs[i2], As[i2] ) )

def LAW61( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '###    INTT = %s\n' % INTT )

        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '###    NP = %s\n' % NP )

        EPs = getData( 'EPs (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        LC2s  = toIntegers( getData( 'LC2s (energy)',  XSS, offset + 2 + 3 * NP, NP ) )

        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e %8d\n' % ( EPs[i2], pdf[i2], cdf[i2], LC2s[i2] ) )
        for i2 in range( NP ) :
            fOut.write( '\n\n' )
            fOut.write( '###     energy_out = %20.12e\n' % EPs[i2] )

            if( LC2s[i2] == 0 ) :
                fOut.write( '###         isotopic\n' )
            else :
                offset = DLW + LC2s[i2] - 1
                JJ = toIntegers( getData( 'JJ (energy)', XSS, offset, 1 ) )[0]
                fOut.write( '###        JJ = %s\n' % INTT )

                NP2 = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
                mus = getData( 'mus (energy)', XSS, offset + 2,           NP2 )
                pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP2,     NP2 )
                cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP2, NP2 )
                for i3 in range( NP2 ) :
                    fOut.write( '        %20.12e %20.12e %8d\n' % ( mus[i3], pdf[i3], cdf[i3] ) )

def LAW66( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    NPSX = toIntegers( getData( 'NPSX (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NPSX = %d\n' % NPSX )
    Ap = getData( 'Ap (energy)', XSS, offset + 1, 1 )[0]
    fOut.write( '### Ap = %20.12e\n' % Ap )

def LAW67( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTMU, NMU = toIntegers( getData( 'INTMU (energy)', XSS, offset, 2 ) )
        offset += 2
        fOut.write( '###   INTMU = %d\n' % INTMU )
        fOut.write( '###   NMU = %d\n' % NMU )

        mus  = getData( 'mus (energy)', XSS, offset, NMU )
        LMU = toIntegers( getData( 'LMU (energy)', XSS, offset + NMU, NMU ) )

        for i2 in range( NMU ) :
            fOut.write( '\n###   mu = %20.12e\n' % mus[i2] )

            offset = DLW + LMU[i2] - 1
            INTEP, NPEP = toIntegers( getData( 'INTEP (energy)', XSS, offset, 2 ) )
            offset += 2
            fOut.write( '###     INTEP = %d\n' % INTEP )
            fOut.write( '###     NPEP = %d\n' % NPEP )

            eps = getData( 'eps (energy)', XSS, offset,            NPEP )
            pdf = getData( 'pdf (energy)', XSS, offset +     NPEP, NPEP )
            cdf = getData( 'cdf (energy)', XSS, offset + 2 * NPEP, NPEP )
            for i3 in range( NPEP ) :
                fOut.write( '        %20.12e %20.12e %8d\n' % ( eps[i3], pdf[i3], cdf[i3] ) )
    

def photonProductionCrossSection( MT, MTPN, SIGP, offset, XSS, energyGrid ) :

    fOut = openFile( MT, 'crossSection', subDir = 'photon', MTP = MTPN )
    MFTYPE = toIntegers( getData( 'MFTYPE', XSS, offset, 1 ) )[0]
    fOut.write( '### MFTYPE = %d\n' % MFTYPE )
    offset += 1

    if( MFTYPE in [ 12, 16 ] ) :
        MTMULT = toIntegers( getData( 'MFTYPE', XSS, offset, 1 ) )[0]
        offset += 1
        fOut.write( '### MTMULT = %d\n' % MTMULT )

        offset += interpolationData( fOut, offset, XSS, 'photonProductionCrossSection' )

        NE = toIntegers( getData( 'NE (photonProductionCrossSection)', XSS, offset, 1 ) )[0]
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (photonProductionCrossSection)', XSS, offset + 1,      NE )
        yields   = getData( 'P(E) (photonProductionCrossSection)',     XSS, offset + 1 + NE, NE )
        for i1 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i1], yields[i1] ) )

    elif( MFTYPE == 13 ) :
        IE, NE = toIntegers( getData( 'cross section info for MT = %d' % MT, XSS, offset, 2 ) )
        energyGrid = energyGrid[IE:IE+NE]
        crossSection = getData( 'cross section for MT = %d' % MT, XSS, offset + 2, NE )
        for i1, energy in enumerate( energyGrid ) : fOut.write( "%20.12e %20.12e\n" % ( energy, crossSection[i1] ) )

    else :
        raise Exception( 'Invalid MFTYPE = %d at location %d' % ( MFTYPE, offset ) )


def dismemberAceFile( fileName ) :

    global cites

    if( args.file ) : print( fileName )

    fIn = open( fileName )
    lines = fIn.readlines( )
    fIn.close( )

    fOut = openFile( 0, 'header' )
    fOut.write( ''.join( lines[:12] ) )
    fOut.write( '\n' )

    NXS  = get8Integers( lines[6] )
    NXS += get8Integers( lines[7] )
    NXS = listBase1( NXS )
    for i1 in range( 1, 16 ) : fOut.write( 'NXS[%2d] = %9d\n' % ( i1, NXS[i1] ) )

    fOut.write( '\n' )
    JXS  = get8Integers( lines[8] )
    JXS += get8Integers( lines[9] )
    JXS += get8Integers( lines[10] )
    JXS += get8Integers( lines[11] )
    JXS = listBase1( JXS )
    for i1 in range( 1, 32 ) : fOut.write( 'JXS[%2d] = %9d\n' % ( i1, JXS[i1] ) )

    fOut.close( )

    XSS = []
    count = NXS[1]
    for i1, line in enumerate( lines[12:] ) :
        while( count == 0 ) : break
        count, lineXSS = get4Floats( line, count, i1 )
        XSS += lineXSS
    fOut = openFile( 0, 'XSS' )
    count = 0
    for XX in XSS :
        if( ( abs( XX ) < NXS[1] ) and ( int( XX ) == XX ) ) :
            count += 1
            fOut.write( '%9d\n' % XX )
        else :
            count += 1
            fOut.write( '%20.12e\n' % XX )
    fOut.close( )
    XSS = listBase1( XSS )
    cites = len( XSS ) * [ 0 ]

    NES = NXS[3]                    # Number of points of energy grid.
    NTR = NXS[4]                    # Number of reactions excluding elastic.
    NR =  NXS[5]                    # Number of reactions having secondary neutrons excluding elastic.

    energyGrid = getData( 'energyGrid', XSS, JXS[1], NES )
    outputXYs1d( 0, 'totalCrossSection',        energyGrid, getData( 'totalCrossSection',       XSS, JXS[1], NES,     NES ) )
    outputXYs1d( 0, 'absorptionCrossSection',   energyGrid, getData( 'absorptionCrossSection',  XSS, JXS[1], NES, 2 * NES ) )
    outputXYs1d( 2, 'crossSection',             energyGrid, getData( 'crossSection',            XSS, JXS[1], NES, 3 * NES ) )
    outputXYs1d( 0, 'averageHeatingNumbers',    energyGrid, getData( 'averageHeatingNumbers',   XSS, JXS[1], NES, 4 * NES ) )

    nu_bar( JXS[2], XSS )
    if( JXS[24] > 0 ) : nu_bar2( 'delayed_nubar', JXS[24], XSS )
    delayedNeutronData( NXS, JXS, XSS )

    MTR = JXS[3]                    # Location of MT array.
    MTs = toIntegers( getData( 'MTs', XSS, MTR, NTR ) )
    LQR = JXS[4]                    # Location of Q-value array.
    Qs = getData( 'Qs', XSS, LQR ,NTR )
    TYR = JXS[5]                    # Location of reaction type array.
    Types = toIntegers( getData( 'Types', XSS, TYR, NTR ) )

    LSIG = JXS[6]                   # Location of table of cross section locators.
    crossSectionLocators = toIntegers( getData( 'crossSectionLocators', XSS, LSIG, NTR ) )
    SIG =  JXS[7]                   # Location of cross sections.

    LAND = JXS[8]                   # Location of table of angular distribution locators.
    angularLocators = toIntegers( getData( 'crossSectionLocators', XSS, LAND, NR + 1 ) )
    AND =  JXS[9]                   # Location of angular distribution.

    LDLW = JXS[10]                   # Location of table of energy distribution locators.
    energyLocators = toIntegers( getData( 'crossSectionLocators', XSS, LDLW, NR ) )
    DLW =  JXS[11]                   # Location of energy distribution.

    MTPsNs = {}
    NTRP  = NXS[6]                  # Number of photon production reactions.
    if( NTRP > 0 ) :
        GPD   = JXS[12]

        if( GPD > 0 ) :
            offset = GPD
            outputXYs1d( 0, 'totalPhotonProductionCrossSection', energyGrid, getData( 'totalPhotonProductionCrossSection', XSS, offset, NES ) )
            offset += NES
            if( False ) :           # Does not seem to exists for data_endfb7.1/O_016_300K.ace
                fOut = openFile( 0, 'outgoingPhotonEnergies' )
                for i1 in range( 30 ) :
                    fOut.write( "\n\n### energy index %d\n" % i1 )
                    outgoingPhotonEnergies = getData( 'outgoingPhotonEnergies', XSS, offset, 20 )
                    for energy in outgoingPhotonEnergies : fOut.write( "    %20.12e\n" % energy )
                    offset += 20
                fOut.close( )

        MTRP  = JXS[13]
        MTPs  = toIntegers( getData( 'MTPs', XSS, MTRP, NTRP ) )
        for MTP in MTPs :
            MT = MTP / 1000
            if( MT not in MTPsNs ) : MTPsNs[MT] = []
            MTPsNs[MT].append( MTP )

        LSIGP = JXS[14]
        photonProductionCrossSectionLocators = toIntegers( getData( 'crossSectionLocators', XSS, LSIGP, NTRP ) )
        SIGP  = JXS[15]

        LANDP = JXS[16]
        photonAngularLocators = toIntegers( getData( 'crossSectionLocators', XSS, LANDP, NTRP ) )
        ANDP  = JXS[17]

        LDLWP = JXS[18]
        photonEnergyLocators = toIntegers( getData( 'crossSectionLocators', XSS, LDLWP, NTRP ) )
        DLWP  = JXS[19]

    if( args.addresses ) : print( "LAND = %d:  AND = %d:  LDLW = %d:  DLW = %d" % ( LAND, AND, LDLW, DLW ) )
#
# Elastic reaction.
#
    fOut = openFile( 2, 'info' )
    fOut.write( 'type = %d\n' % -1 )
    fOut.write( 'Q = %20.12e\n' % 0.0 )
    fOut.close( )

    dismemberAngular( 0, 2, -1, angularLocators, AND, XSS, 'neutron' )         # Special case for elastic.

#
# Non-elastic reaction.
#
    for i1, MT in enumerate( MTs ) :
        if( args.verbose > 1 ) : print( "MT = %3d" % MT )
        locator = SIG + crossSectionLocators[i1] - 1

        fOut = openFile( MT, 'info' )
        fOut.write( 'type = %d\n' % Types[i1] )
        fOut.write( 'Q = %20.12e\n' % Qs[i1] )
        if( args.addresses ) : fOut.write( 'locator = %d\n' % locator )
        line, column = lineColumn1BaseForXSS( locator )
        if( args.addresses ) : fOut.write( 'line, column = %d, %d\n' % ( line, column ) )
        fOut.close( )

        offset, length = toIntegers( getData( 'cross section info for MT = %d' % MT, XSS, locator, 2 ) )
        crossSection = getData( 'cross section for MT = %d' % MT, XSS, locator + 2, length )
        outputXYs1d( MT, 'crossSection', energyGrid, crossSection, offset = offset, length = length )

        if( 4 < MT < 100 ) :
            dismemberAngular( i1 + 1, MT, Types[i1], angularLocators, AND, XSS, 'neutron' )
            dismemberEnergy( 0, MT, Types[i1], DLW, energyLocators[i1], XSS, 'neutron' )

        TY = abs( Types[i1] )
        if( TY > 100 ) :
            fOut = openFile( MT, 'multiplicityVsEnergy', subDir = 'neutron' )

            offset = DLW + TY - 101
            offset += interpolationData( fOut, offset, XSS, 'photonProductionCrossSection' )

            NE = toIntegers( getData( 'NE (photonProductionCrossSection)', XSS, offset, 1 ) )[0]
            fOut.write( '### NE = %d\n' % NE )
            energies = getData( 'energies (photonProductionCrossSection)', XSS, offset + 1,      NE )
            yields   = getData( 'P(E) (photonProductionCrossSection)',     XSS, offset + 1 + NE, NE )
            for i1 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i1], yields[i1] ) )

            fOut.close( )

    for MT in sorted( MTPsNs ) :
        MTPNs = MTPsNs[MT]
        for MTPN in MTPNs :
            i1 = MTPs.index( MTPN )
            locator = SIGP + photonProductionCrossSectionLocators[i1] - 1
            photonProductionCrossSection( MT, MTPN, SIGP, locator, XSS, energyGrid )
            dismemberAngular( i1, MT, 0, photonAngularLocators, ANDP, XSS, 'photon', MTP = MTPN )
            dismemberEnergy( 0, MT, 1, DLWP, photonEnergyLocators[i1], XSS, 'photon', MTP = MTPN )

    YP = JXS[20]
    if( YP > 0 ) :
        NYP = toIntegers( getData( 'YP', XSS, YP, 1 ) )[0]
        MTYs = toIntegers( getData( 'MTY', XSS, YP + 1, NYP ) )
        fOut = openFile( 0, 'YP' )
        fOut.write( 'NYP = %d\n' % NYP )
        for MTY in MTYs :fOut.write( '%d\n' % MTY )
        fOut.close( )

    LUNR = JXS[23]              # Location of probability tables
    if( LUNR > 0 ) :
        fOut = openFile( 0, 'URR_probabilityTables' )
        offset = LUNR
        N, M, INT, ILF, IOA, IFF = toIntegers( getData( 'N (URR)', XSS, LUNR, 6 ) )
        offset += 6
        fOut.write( '### N = %d\n'   % N )
        fOut.write( '### M = %d\n'   % M )
        fOut.write( '### INT = %d\n' % INT )
        fOut.write( '### ILF = %d\n' % ILF )
        fOut.write( '### IOA = %d\n' % IOA )
        fOut.write( '### IFF = %d\n' % IFF )

        energies = getData( 'energies (URR)', XSS, offset, N )
        offset += N
        fOut.write( '### Energies\n' )
        for energy in energies : fOut.write( '%20.12e\n' % energy )

        Ps = getData( 'P(i,j,k) (URR)', XSS, offset, N * 6 * M )
        fOut.write( '\n\n### Probality tables\n' )
        for pt in Ps : fOut.write( '%20.12e\n' % pt )

        fOut.close( )

    fOut = openFile( 0, 'cites' )
    uncited = 0
    multiplyCited = 0
    officialEnd = JXS[22]
    pastEndUncited = 0
    for i1, cite in enumerate( cites ) :
        if( cite == 0 ) : uncited += 1
        if( cite > 1 ) : multiplyCited += 1
        if( ( i1 >= officialEnd ) and ( cite == 0 ) ) : pastEndUncited += 1
        fOut.write( '%d\n' % cite )
    fOut.close( )

    if( ( uncited > 0 ) or not( args.quiet ) ) :
        if( not( args.Official and ( ( uncited - pastEndUncited ) == 0 ) ) ) :
            if( args.verbose > 0 ) : print( "    Official uncited data elements = %d, total uncited = %d, multiply cited = %d" % 
                ( ( uncited - pastEndUncited ), uncited, multiplyCited ) )

if( __name__ == '__main__' ) :

    if( not( args.keep ) ) : os.system( 'rm -rf %s/*' % args.output )
    dismemberAceFile( args.ace )
