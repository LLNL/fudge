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

import copy, sys

import xData.standards as standardsModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule

FUDGE_EPS = 1e-8

def parseENDFByMT_MF( fileName, stripMATMFMTCount = True, logFile = sys.stderr ) :
    """This function reads an endf file and creates a python dictionary with the keys being MT
    numbers and the values being the data for the MT key. For each MT, the data are stored in a python
    dictionary with the keys being MF numbers and the values being the data for the MF key."""

    MTDatas = {}
    f = open( fileName )
    header = f.readline( )[:-1]           # Occasionally, the first line does not have MT or MF values.
    line = 0
    endlLine = f.readline( )
    secondLine = endlLine
    while( True ) :
        if( endlLine == '' ) : break
        line += 1
        try :
            MF, MT = int( endlLine[70:72] ), int( endlLine[72:75] )
        except :
            logFile.write( "%s\n" % line )
            raise
        if( MT not in MTDatas ) : MTDatas[MT] = {}
        MTData = MTDatas[MT]
        if( MF not in MTData ) : MTData[MF] = []
        if( stripMATMFMTCount ) :
            MTData[MF].append( endlLine[:66].rstrip( ) )
        else :
            MTData[MF].append( endlLine.rstrip( ) )
        endlLine = f.readline( )
    f.close( )
    try :
        MAT = int( secondLine[66:70] )
    except ValueError:
        MAT = int( MTDatas[451][1][6][31:35] )
    del MTDatas[0]
    return( header, MAT, MTDatas )

def ENDFInterpolationToGNDS1d( interpolation ) :

    if( interpolation == 1 ) : return( standardsModule.interpolation.flatToken )
    if( interpolation == 2 ) : return( standardsModule.interpolation.linlinToken )
    if( interpolation == 3 ) : return( standardsModule.interpolation.linlogToken )
    if( interpolation == 4 ) : return( standardsModule.interpolation.loglinToken )
    if( interpolation == 5 ) : return( standardsModule.interpolation.loglogToken )
    if( interpolation == 6 ) : return( standardsModule.interpolation.chargedParticleToken )
    raise Exception( 'Unsupport 2d interpolation = %d' % interpolation )

def ENDFInterpolationToGNDS2plusd( interpolation ) :

    if(    0 < interpolation < 10 ) :
        interpolationQualifier = standardsModule.interpolation.noneQualifierToken
        normalInterpolation = interpolation
    elif( 10 < interpolation < 20 ) :
        interpolationQualifier = standardsModule.interpolation.correspondingPointsToken
        normalInterpolation = interpolation - 10
    elif( 20 < interpolation < 30 ) :
        interpolationQualifier = standardsModule.interpolation.unitBaseToken
        normalInterpolation = interpolation - 20
    else :
        raise Exception( 'Unsupport 2+d interpolation = %d' % interpolation )
    return( interpolationQualifier, ENDFInterpolationToGNDS1d( normalInterpolation ) )

def funkyFloatStringToFloat( index, s, logFile = sys.stderr ) :

    index6 = index % 6
    value = s[11 * index6: 11 * ( index6 + 1 )]
    try :
        f = float( value )
    except :
        value = value.lower().replace('d','')
        value = value.replace( ' ', '' )
        i1 = value.find( '.' )
        if( i1 != -1 ) :
            j = value[i1:].find( '+' )
            if( j == -1 ) : j = value[i1:].find( '-' )
            if( j == -1 ) : raise Exception( "Float value at index = %d (%d) is not funky enough <%s>" % ( index, index6, s ) )
            value = value[:i1+j] + 'e' + value[i1+j:]
        try :
            f = float( value )
        except :
            if( value == len( value ) * ' ' ) : return( 0 )
            if( logFile is not None ) : logFile.write( 'Could not convert value "%s" at index %s (%s) of string "%s"\n' % ( value, index6, index, s ) )
            raise
    return( f )

def sixFunkyFloatStringsToFloats( s, logFile = sys.stderr ) :

    return( [ funkyFloatStringToFloat( i1, s, logFile = logFile ) for i1 in xrange( 6 ) ] )

def sixFunkyFloatStringsToIntsAndFloats( s, intIndices = [], logFile = sys.stderr ) :

    values = sixFunkyFloatStringsToFloats( s, logFile = logFile )
    for i1 in intIndices : values[i1] = int( values[i1] )
    return( values )

def nFunkyFloatStringsToFloats( n1, startLine, lines, dimension = 1, logFile = sys.stderr ) :

    floats = []
    nd = dimension * n1
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i1 in xrange( n6 ) : floats.append( funkyFloatStringToFloat( i1, l, logFile = logFile ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1
    if( dimension > 1 ) :
        floats1D, floats = floats, []
        for i1 in xrange( 0, dimension * n1, dimension ) : floats.append( floats1D[i1:i1+dimension] )
    return( floats )

def readEndfINTG( string, ndigit ) :
    """Special case in ENDF: line of all integers, used to represent sparse correlation matrix."""

    ix = int( string[:5] )
    iy = int( string[5:10] )

    try :
        i1, i2 = { 2 : [ 11, 65 ], 3 : [ 11, 63 ], 4 : [ 11, 66 ], 5 : [ 11, 65 ], 6 : [ 10, 66 ] }[ndigit]
    except :
        raise ValueError, ("Encountered illegal ndigit (%d) for INTG format" % ndigit)
    sub = string[i1:i2]
    n1 = ndigit + 1
    vals = [ sub[n1*i1:n1*(i1+1)] for i1 in range( len( sub ) // n1 ) ]

    def parseInt( string ) :

        return( 0 if not string.strip() else int( string ) )

    arr = [ parseInt( a ) for a in vals ]
    return( ix, iy, arr )

def readDiscreteAndLegendre( numDiscrete, numContinua, startLine, lines, dimension, logFile = sys.stderr ) :

    floats = []
    nd = dimension * numDiscrete + dimension * numContinua  #  total amount of data
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i1 in xrange( n6 ) : floats.append( funkyFloatStringToFloat( i1, l, logFile = logFile ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1

    discretes = []                                  # Data for discrete gammas is [ gammaEnergy, multiplicity ]
    for i1 in xrange( 0, dimension * numDiscrete, dimension ) :
        discretes.append( [ floats[i1], floats[i1+1:i1+dimension] ] )
    continua = []
    offset = dimension * numDiscrete
    for i1 in xrange( offset, offset + dimension * numContinua, dimension ) :
        continua.append( floats[i1:i1+dimension] )

    return( discretes, continua )

def nStringsToInts( n_, startLine, lines, dimension = 1 ) :

    ints = []
    nd = dimension * n_
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i1 in xrange( n6 ) : ints.append( int( l[11*i1: 11*(i1+1)] ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1
    if( dimension > 1 ) :
        ints1D, ints = ints, []
        for i1 in xrange( 0, dimension * n_, dimension ) : ints.append( ints1D[i1:i1+dimension] )
    return( ints )

def getENDFDate( date ) :

    try :
        tmp = int(date[4:]) # can it be converted to int?
        y, m, d = date[4:6], date[6:8], date[8:]
        if( int( y ) < 50 ) :
            y = '20' + y
        else :
            y = '19' + y
        date = '%s-%s-%s' % ( y, m, d )
    except :
        sMonth = date[5:8]
        try :
            month = { "JAN" : '01', "FEB" : '02', "MAR" : '03', "APR" : '04', "MAY" : '05', "JUN" : '06',
                "JUL" : '07', "AUG" : '08', "SEP" : '09', "OCT" : '10', "NOV" : '11', "DEC" : '12' }[sMonth.upper( )]
        except :
            raise Exception( 'Bad evaluation date: could not convert month = "%s" to integer for "%s". See option --ignoreBadDate' % ( sMonth, date ) )
        try :
            year = 1900 + int( date[8:10] )
            if( year < 1950 ) : year += 100
        except :
            raise Exception( 'Bad evaluation date: could not convert year = "%s" to integer for "%s". See option --ignoreBadDate' % ( date[8:10], date ) )
        day = '01'
        date = '%s-%s-%s' % ( year, month, day )
    return( date )

def getTAB1( startLine, dataLines, logFile = sys.stderr ) :

    lineNumber = startLine
    C1, C2, L1, L2, NR, NP = sixFunkyFloatStringsToFloats( dataLines[startLine], logFile = logFile )
    NR = int( NR )                  # number of interpolation flags
    NP = int( NP )                  # number of pairs of data
    lineNumber += 1
    interpolationInfo = nStringsToInts( NR, lineNumber, dataLines, dimension = 2 )
    if( interpolationInfo[0][0] == NP ) : interpolationInfo = interpolationInfo[:1] # There is realy only one region, so fix it.

    lineNumber += ( NR + 2 ) / 3
    data = nFunkyFloatStringsToFloats( NP, lineNumber, dataLines, dimension = 2 )

    lineNumber += ( NP + 2 ) / 3
    return( lineNumber, { 'C1' : C1, 'C2' : C2, 'L1' : L1, 'L2' : L2, 'NR' : NR, 'NP' : NP, 'interpolationInfo' : interpolationInfo, 'data' : data } )

def getList( startLine, dataLines, logFile = sys.stderr ) :

    lineNumber = startLine
    C1, C2, L1, L2, NPL, N2 = sixFunkyFloatStringsToFloats( dataLines[startLine], logFile = logFile )
    NPL = int( NPL )                # number of items in list
    N2 = int( N2 )                  # 

    lineNumber += 1
    data = nFunkyFloatStringsToFloats( NPL, lineNumber, dataLines, dimension = 1 )

    lineNumber += ( NPL + 5 ) / 6
    return( lineNumber, { 'C1' : C1, 'C2' : C2, 'L1' : L1, 'L2' : L2, 'NPL' : NPL, 'N2' : N2, 'data' : data } )

def getTAB2Header( startLine, dataLines, logFile = sys.stderr ) :

    lineNumber = startLine
    C1, C2, L1, L2, NR, NZ = sixFunkyFloatStringsToFloats( dataLines[startLine], logFile = logFile )
    L1, L2 = int( L1 ), int( L2 )
    NR = int( NR )                  # number of interpolation flags
    NZ = int( NZ )                  # number of TAB1's

    lineNumber += 1
    interpolationInfo = nStringsToInts( NR, lineNumber, dataLines, dimension = 2 )

    lineNumber += ( NR + 2 ) / 3
    return( lineNumber, { 'C1' : C1, 'C2' : C2, 'L1' : L1, 'L2' : L2, 'NR' : NR, 'NZ' : NZ, 'interpolationInfo' : interpolationInfo } )

def getTAB2_TAB1s( startLine, dataLines, logFile = sys.stderr, axes = None ) :
    """
    This function currently does not support multiple interpolation regions or discontinuous functions vs. C2.
    If this is changed, some calling functions will need to be changed.
    """

    lineNumber, TAB2 = getTAB2Header( startLine, dataLines, logFile = logFile )
    TAB2s = []
    iStart, C2Prior = 0, -1e99          # C2 may be mu that starts at -1.
    piecewise = None
    for i1, ( NZSub, interpolation ) in enumerate( TAB2['interpolationInfo'] ) :
        TAB1s = []
        if( piecewise is not None ) : TAB1s.append( piecewise )
        for i2 in xrange( iStart, NZSub ) :
            lineNumber, TAB1, piecewise = getTAB1Regions( lineNumber, dataLines, logFile = logFile, dimension = 2, axes = axes )
            if( C2Prior == TAB1['C2'] ) : raise 'hell - see document'
            TAB1s.append( piecewise )
        iStart = NZSub
        TAB2s.append( [ interpolation, TAB1s ] )
    TAB2['TAB2s'] = TAB2s
    return( lineNumber, TAB2 )

def getTAB2_Lists( startLine, dataLines, logFile = sys.stderr ) :

    lineNumber, TAB2 = getTAB2Header( startLine, dataLines, logFile = logFile )
    Lists = []
    for i1 in xrange( TAB2['NZ'] ) :
        lineNumber, List = getList( lineNumber, dataLines, logFile = logFile )
        Lists.append( List )
    TAB2['Lists'] = Lists
    return( lineNumber, TAB2 )

def getTAB1Regions( startLine, dataLines, allowInterpolation6 = False, logFile = sys.stderr, dimension = 1, axes = None,
        cls = XYsModule.XYs1d ) :

    endLine, TAB1 = getTAB1( startLine, dataLines, logFile = logFile )
    n1, i1, mode, data = TAB1['NR'], 0, 0, TAB1['data']

    for regionsEnd, interpolation in TAB1['interpolationInfo'] :        # Search for 'flat' region next to linear region.
        if( interpolation not in [ 1, 2 ] ) : break
    if( interpolation in [ 1, 2 ] ) :                                   # Only consider if all regions are flat and/or linear.
        regionsStart, interpolationRegions2 = 0, []                     # We are now going to combine as many adjacent flat and linear regions as possible.
        for index, interpolationInfo in enumerate( TAB1['interpolationInfo'] ) :
            regionsEnd, interpolation = interpolationInfo               # This loop converts all flat regions to linear if possible
            interpolationRegions2.append( [ regionsEnd, interpolation ] )   # (i.e. x_i < x_{i+1} and all y values are the same).
            if( interpolation == 1 ) :
                convertToLinear = True
                x1, y1 = data[regionsStart]
                for x2, y2 in data[regionsStart+1:regionsEnd] :
                    if( ( x1 == x2 ) or ( y1 != y2 ) ) :                # x_i < x_{i+1} and all y values must be the same.
                        convertToLinear = False
                        break
                    x1 = x2
                if( convertToLinear ) : interpolationRegions2[-1][1] = 2
            regionsStart = regionsEnd
        x1, i1 = interpolationRegions2[0]
        interpolationRegions = []
        for x2, i2 in interpolationRegions2[1:] :                       # Now combine adjacent linear regions if possible
            if( ( i1 == 2 ) and ( i2 == 2 ) ) :
                if( data[x1][0] == data[x1+1][0] ) : interpolationRegions.append( [ x1, i1 ] )
            else :
                interpolationRegions.append( [ x1, i1 ] )
            x1, i1 = x2, i2
        interpolationRegions.append( [ x1, i1 ] )
    else :
        interpolationRegions = TAB1['interpolationInfo']

    if( interpolationRegions != TAB1['interpolationInfo'] ) :           # Print a warning message if number of regions were reduced.
# BRB, check this with real data.
        if( not( ( len( TAB1['interpolationInfo'] ) == 1 ) and ( TAB1['interpolationInfo'][0][1] == 1 ) and ( interpolationRegions[0][1] == 2 ) ) ) :
            logFile.write( ' reduced regions from %s to %s' % ( len( TAB1['interpolationInfo'] ), len( interpolationRegions ) ) )

    if (data[-2][0] == data[-1][0] and data[-1][1] == 0):               # discontinuity at final outgoing energy used to drop spectrum to 0
        if (interpolation == 1):
            data.pop(-2)                                                # flat interpolation: omit 2nd-to-last point
        else:
            data[-1][0] *= 1.000000001                                  # lin-lin interpolation: blur edge of discontinuity

    n1, regions = 0, []
    for interpolationInfo in interpolationRegions :
        n2, interpolation = interpolationInfo
        if( n2 == n1 ) :
            logFile.write( "encountered interpolation region with no points" )
            continue
        if( interpolation == 6 ) :
            threshold = 0.
            if( not( allowInterpolation6 ) ) : raise Exception( "charged-particle interpolation not allowed here" )
        np1 = n1
        if( n1 > 0 ) :
            if( data[n1-1][0] != data[n1][0] ) : np1 -= 1

        regionData2, regionData, x1 = data[np1:n2], [], None
        for x2, y2 in regionData2 :                                     # Loop to remove duplicate points.
            if( x1 is not None ) :
                if( ( x1 == x2 ) and ( y1 == y2 ) ) : continue
            regionData.append( [ x2, y2 ] )
            x1, y1 = x2, y2

        while( len( regionData ) > 0 ) :
            i3, n3 = 1, len( regionData )
            x1, y1 = regionData[0]
            x2, y2 = regionData[1]
            while( ( i3 < n3 ) and ( x1 != x2 ) ) :
                x1, y1 = x2, y2
                i3 += 1
                if( i3 < n3 ) : x2, y2 = regionData[i3]
            value = None
            if( dimension > 1 ) : value = TAB1['C2']
            regions.append( cls( data = regionData[:i3], value = value,
                    interpolation = ENDFInterpolationToGNDS1d( interpolation ), axes = axes ) )
            regionData = regionData[i3:]
            for i3 in xrange( 1, len( regionData ) ) :
                if( regionData[i3][0] != x1 ) : break
            i3 -= 1
            regionData = regionData[i3:]
        n1 = n2
    return( endLine, TAB1, regions )

def niceSortOfMTs( MTs, verbose = True, logFile = sys.stderr ) :

    def removeGetIfPresent( MT, MTs ) :

        if( MT not in MTs ) : return( [] )
        MTs.remove( MT )
        return( [ MT ] )

    MTs = copy.deepcopy( MTs )
    newMTs = removeGetIfPresent(   2, MTs )
    for MT in xrange( 50, 92 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent(   4, MTs )
    newMTs += removeGetIfPresent(  16, MTs )        # (z,2n) reactions
    for MT in xrange( 875, 892 ) : newMTs += removeGetIfPresent( MT, MTs )

    newMTs += removeGetIfPresent(  17, MTs )
    newMTs += removeGetIfPresent(  37, MTs )
    newMTs += removeGetIfPresent(  18, MTs )
    newMTs += removeGetIfPresent(  19, MTs )
    newMTs += removeGetIfPresent(  20, MTs )
    newMTs += removeGetIfPresent(  21, MTs )
    newMTs += removeGetIfPresent(  38, MTs )
    newMTs += removeGetIfPresent(  28, MTs )
    newMTs += removeGetIfPresent(  32, MTs )
    newMTs += removeGetIfPresent(  33, MTs )

    for MT in xrange( 600, 650 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent( 103, MTs )          # (z,p) reactions

    for MT in xrange( 650, 700 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent( 104, MTs )          # (z,d) reactions

    for MT in xrange( 700, 750 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent( 105, MTs )          # (z,t) reactions

    for MT in xrange( 750, 800 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent( 106, MTs )          # (z,He3) reactions

    for MT in xrange( 800, 850 ) : newMTs += removeGetIfPresent( MT, MTs )
    newMTs += removeGetIfPresent( 107, MTs )          # (z,a) reactions

    newMTs += removeGetIfPresent( 102, MTs )          # (z,g) reactions

    MT5 = removeGetIfPresent( 5, MTs )                 # (z,everything else)

    MTAtomics = []
    for MT in range( 500, 573 ) : MTAtomics += removeGetIfPresent( MT, MTs )

    skippingMTs = []
    for MT in [ 10, 27, 101, 151 ] : skippingMTs += removeGetIfPresent( MT, MTs )
    for MT in xrange( 201, 600 ) : skippingMTs += removeGetIfPresent( MT, MTs )
    for MT in xrange( 850, 875 ) : skippingMTs += removeGetIfPresent( MT, MTs )

    if( verbose and ( len( skippingMTs ) > 0 ) ) : logFile.write( 'Skipping MTs = %s\n' % skippingMTs )

    newMTs += MTs + MT5 + MTAtomics
    return( newMTs )

def getMFDataInMFList( MFs, MFData ) :

    Data = {}
    for MF in MFData :
        if( MF in MFs ) : Data[MF] = MFData[MF]
    return( Data )

def toEnergyFunctionalData( info, dataLine, MF5Data, LF, moniker, unit, xLabel = 'energy_in', xUnit = 'eV' ) :

    import fudge.gnds.productData.distributions.energy as energyModule

    axes = axesModule.axes( labelsUnits = { 0 : ( moniker , unit ), 1 : ( xLabel, xUnit ) } )
    dataLine, TAB1, data = getTAB1Regions( dataLine, MF5Data, logFile = info.logs, axes = axes )
    EFclass = None
    for tmp in ( energyModule.a, energyModule.b, energyModule.theta, energyModule.g, energyModule.T_M ):
        if ( tmp.moniker == moniker ):
            EFclass = tmp

    if( len( data ) == 1 ) :
        oneD = data[0]
    else :
        oneD = regionsModule.regions1d( axes = data[0].axes )
        for datum in data : oneD.append( datum )

    return( dataLine, TAB1, EFclass( oneD ) )
