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

import copy, sys
from fudge import gnd
from fudge.core.math.xData import axes, XYs

ENDF_Accuracy = 1e-3
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
    except :
        MAT = int( MTDatas[451][1][6][31:35] )
    del MTDatas[0]
    return( header, MAT, MTDatas )

def interpolationString( interpolation ) :

    return( "%s,%s" % ENDFInterpolationToGND2d( interpolation ) )

def ENDFInterpolationToGND2d( interpolation ) :

    ENDFInterpolationToGND = { 2 : axes.linearToken, 3 : axes.logToken }
    if( interpolation == 1 ) : return( ENDFInterpolationToGND[2], axes.flatToken )
    if( interpolation == 2 ) : return( ENDFInterpolationToGND[2], ENDFInterpolationToGND[2] )
    if( interpolation == 3 ) : return( ENDFInterpolationToGND[3], ENDFInterpolationToGND[2] )
    if( interpolation == 4 ) : return( ENDFInterpolationToGND[2], ENDFInterpolationToGND[3] )
    if( interpolation == 5 ) : return( ENDFInterpolationToGND[3], ENDFInterpolationToGND[3] )
    if( interpolation == 6 ) : return( ENDFInterpolationToGND[2], axes.chargedParticleToken )
    raise Exception( 'Unsupport 2d interpolation = %d' % interpolation )

def ENDFInterpolationToGND3plusd( interpolation ) :

    if(    0 < interpolation < 10 ) :
        header = ''
        normalInterpolation = interpolation
    elif( 10 < interpolation < 20 ) :
        header = 'correspondingPoints:'
        normalInterpolation = interpolation - 10
    elif( 20 < interpolation < 30 ) :
        header = 'unitBase:'
        normalInterpolation = interpolation - 20
    else :
        raise Exception( 'Unsupport 3+d interpolation = %d' % interpolation )
    sx, sy = ENDFInterpolationToGND2d( normalInterpolation )
    return( "%s%s,%s" % ( header, sx, sy ) )

def ENDFInterpolationToGNDAxes3plusd( interpolation ) :

    if(    0 < interpolation < 10 ) :
        interpolationQualifier = None
        normalInterpolation = interpolation
    elif( 10 < interpolation < 20 ) :
        interpolationQualifier = 'correspondingPoints'
        normalInterpolation = interpolation - 10
    elif( 20 < interpolation < 30 ) :
        interpolationQualifier = 'unitBase'
        normalInterpolation = interpolation - 20
    else :
        raise Exception( 'Unsupport 3+d interpolation = %d' % interpolation )
    dependentInterpolation, independentInterpolation = ENDFInterpolationToGND2d( normalInterpolation )
    return( dependentInterpolation, independentInterpolation, interpolationQualifier )

def funkyFloatStringToFloat( index, s, logFile = sys.stderr ) :

    index6 = index % 6
    value = s[11 * index6: 11 * ( index6 + 1 )]
    try :
        f = float( value )
    except :
        value = value.replace( ' ', '' )
        i = value.find( '.' )
        if( i != -1 ) :
            j = value[i:].find( '+' )
            if( j == -1 ) : j = value[i:].find( '-' )
            if( j == -1 ) : raise Exception( "Float value at index = %d (%d) is not funky enough <%s>" % ( index, index6, s ) )
            value = value[:i+j] + 'e' + value[i+j:]
        try :
            f = float( value )
        except :
            if( value == len( value ) * ' ' ) : return( 0 )
            logFile.write( 'Could not convert value "%s" at index %s (%s) of string "%s"\n' % ( value, index6, index, s ) )
            raise
    return( f )

def sixFunkyFloatStringsToFloats( s, logFile = sys.stderr ) :

    return( [ funkyFloatStringToFloat( i, s, logFile = logFile ) for i in xrange( 6 ) ] )

def sixFunkyFloatStringsToIntsAndFloats( s, intIndices = [], logFile = sys.stderr ) :

    values = sixFunkyFloatStringsToFloats( s, logFile = logFile )
    for i in intIndices : values[i] = int( values[i] )
    return( values )

def nFunkyFloatStringsToFloats( n, startLine, lines, dimension = 1, logFile = sys.stderr ) :

    floats = []
    nd = dimension * n
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i in xrange( n6 ) : floats.append( funkyFloatStringToFloat( i, l, logFile = logFile ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1
    if( dimension > 1 ) :
        floats1D, floats = floats, []
        for i in xrange( 0, dimension * n, dimension ) : floats.append( floats1D[i:i+dimension] )
    return( floats )

def readEndfINTG( string, ndigit ) :
    """Special case in ENDF: line of all integers, used to represent sparse correlation matrix."""

    x = int(string[:5])
    y = int(string[5:10])

    if ndigit==2:
        sub = string[11:65]
        vals = [sub[3*i:3*i+3] for i in range(len(sub)//3)]
    elif ndigit==3:
        sub = string[11:63]
        vals = [sub[4*i:4*i+4] for i in range(len(sub)//4)]
    elif ndigit==4:
        sub = string[11:66]
        vals = [sub[5*i:5*i+5] for i in range(len(sub)//5)]
    elif ndigit==5:
        sub = string[11:65]
        vals = [sub[6*i:6*i+6] for i in range(len(sub)//6)]
    elif ndigit==6:
        sub = string[10:66]
        vals = [sub[7*i:7*i+7] for i in range(len(sub)//7)]
    else:
        raise ValueError, ("Encountered illegal ndigit (%i) for INTG format" % ndigit)

    def parseInt( string ):
        return (0 if not string.strip() else int(string))
    arr = [parseInt(a) for a in vals]
    return x, y, arr

def readDiscreteAndLegendre( numDiscrete, numContinua, startLine, lines, dimension = 1, logFile = sys.stderr ) :

    floats = []
    nd = 2*numDiscrete + dimension * numContinua   #  total amount of data
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i in xrange( n6 ) : floats.append( funkyFloatStringToFloat( i, l, logFile = logFile ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1

    #   the data for discrete gammas is [ gammaEnergy, multiplicity ]       ???? Why is this???
    discretes = []
    for i in xrange( 0, 2*numDiscrete, 2 ) :
        discretes.append( [ floats[i], floats[i+1:i+2] ] )
    continua = []
    offset = 2*numDiscrete
    for i in xrange( offset, offset + dimension*numContinua, dimension ) :
        continua.append( floats[i:i+dimension] )

    return( discretes, continua )

def nStringsToInts( n, startLine, lines, dimension = 1 ) :

    ints = []
    nd = dimension * n
    while( True ) :
        l = lines[startLine]
        n6 = min( nd, 6 )
        for i in xrange( n6 ) : ints.append( int( l[11*i: 11*(i+1)] ) )
        nd -= n6
        if( nd == 0 ) : break
        startLine += 1
    if( dimension > 1 ) :
        ints1D, ints = ints, []
        for i in xrange( 0, dimension * n, dimension ) : ints.append( ints1D[i:i+dimension] )
    return( ints )

def getENDFDate( date ) :

    try :
        tmp = int(date[4:]) # can it be converted to int?
        y,m,d = date[4:6], date[6:8], date[8:]
        if int(y)<50: y = '20'+y
        else: y = '19'+y
        date = '%s-%s-%s' % (y,m,d)
    except :
        sMonth = date[5:8]
        try :
            month = { "JAN" : '01', "FEB" : '02', "MAR" : '03', "APR" : '04', "MAY" : '05', "JUN" : '06',
                "JUL" : '07', "AUG" : '08', "SEP" : '09', "OCT" : '10', "NOV" : '11', "DEC" : '12' }[sMonth.upper( )]
        except :
            raise Exception( 'Bad evaluation date: could not convert month = "%s to integer for %s' % ( sMonth, date ) )
        try :
            year = 1900 + int( date[8:10] )
            if( year < 1950 ) : year += 100
        except :
            raise Exception( 'Bad evaluation date: could not convert year = "%s for %s" to integer' % date )
        day = '01'
        date = '%s-%s-%s' % (year,month,day)
        #date = 100 * ( 100 * year + month ) + day
    return( date )

def getTAB1( startLine, dataLines, oneNR = True, logFile = sys.stderr ) :

    lineNumber = startLine
    C1, C2, L1, L2, NR, NP = sixFunkyFloatStringsToFloats( dataLines[startLine], logFile = logFile )
    NR = int( NR )                  # number of interpolation flags
    NP = int( NP )                  # number of pairs of data
    lineNumber += 1
    interpolationInfo = nStringsToInts( NR, lineNumber, dataLines, dimension = 2 )
    if( oneNR and ( NR > 1 ) ) :
        if( interpolationInfo[0][0] != NP ) :           # If True crash; otherwise, fix bad data below.
            logFile.write( "%s\n" % startLine )
            logFile.write( dataLines[startLine] + '\n' )
            logFile.write( "%s\n" % interpolationInfo )
            raise Exception( "Currently only one interpolation flag is supported" )
        interpolationInfo = interpolationInfo[:1]

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

def getTAB2_TAB1s( startLine, dataLines, asRegions = False, axes = None, logFile = sys.stderr ) :

    lineNumber, TAB2 = getTAB2Header( startLine, dataLines, logFile = logFile )
    TAB1s = []
    for i in xrange( TAB2['NZ'] ) :
        if( asRegions ) :
            lineNumber, TAB1, piecewise = getTAB1Regions( lineNumber, dataLines, axes, logFile = logFile )
            TAB1s.append( [ TAB1['C2'], piecewise ] )
        else :
            lineNumber, TAB1 = getTAB1( lineNumber, dataLines, oneNR = True, logFile = logFile )
            TAB1s.append( TAB1 )
    TAB2['TAB1s'] = TAB1s
    return( lineNumber, TAB2 )

def getTAB2_Lists( startLine, dataLines, logFile = sys.stderr ) :

    lineNumber, TAB2 = getTAB2Header( startLine, dataLines, logFile = logFile )
    Lists = []
    for i in xrange( TAB2['NZ'] ) :
        lineNumber, List = getList( lineNumber, dataLines, logFile = logFile )
        Lists.append( List )
    TAB2['Lists'] = Lists
    return( lineNumber, TAB2 )

def getTAB1Regions( startLine, dataLines, axes_, oneNR = False, allowInterpolation6 = False, logFile = sys.stderr ) :

    endLine, TAB1 = getTAB1( startLine, dataLines, oneNR = oneNR, logFile = logFile )
    n, i1, mode, data = TAB1['NR'], 0, 0, TAB1['data']

    for regionsEnd, interpolation in TAB1['interpolationInfo'] :        # Search for 'flat' region next to linear region.
        if( interpolation not in [ 1, 2 ] ) : break
    if( interpolation in [ 1, 2 ] ) :                                   # Only consider if all regions are flat and/or linear.
        regionsStart, interpolationRegions2 = 0, []                     # We are now going to combine as many adjacent flat and linear regions as possible.
        for index, interpolationInfo in enumerate( TAB1['interpolationInfo'] ) :
            regionsEnd, interpolation = interpolationInfo               # This loop converts all flat regions to linear if possible
            interpolationRegions2.append( [ regionsEnd, interpolation ] )   # (i.e. x_i < x_{i+1} and all y values are the same.
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
        if( not( ( len( TAB1['interpolationInfo'] ) == 1 ) and ( TAB1['interpolationInfo'][0][1] == 1 ) and ( interpolationRegions[0][1] == 2 ) ) ) :
            logFile.write( 'reduced regions from %s to %s' % ( len( TAB1['interpolationInfo'] ), len( interpolationRegions ) ) )

    n1, regions, index = 0, [], 0
    for interpolationInfo in interpolationRegions :
        n2, interpolation = interpolationInfo
        if n2==n1:
            logFile.write( "encountered interpolation region with no points" )
            continue
        if( ( interpolation == 6 ) and not( allowInterpolation6 ) ) : raise Exception( "charged-particle interpolation not allowed here" )
        np1 = n1
        if( n1 > 0 ) :
            if( data[n1-1][0] != data[n1][0] ) : np1 -= 1
        xi, yi = ENDFInterpolationToGND2d( interpolation )
        axes_[0].setInterpolation( axes.interpolationXY( xi, yi ) )

        regionData2, regionData, x1 = data[np1:n2], [], None
        for x2, y2 in regionData2 :                                     # Loop to remove duplicate points.
            if( x1 is not None ) :
                if( ( x1 == x2 ) and ( y1 == y2 ) ) : continue
            regionData.append( [ x2, y2 ] )
            x1, y1 = x2, y2

        while( len( regionData ) > 0 ) :
            i, n = 1, len( regionData )
            x1, y1 = regionData[0]
            x2, y2 = regionData[1]
            while( ( i < n ) and ( x1 != x2 ) ) :
                x1, y1 = x2, y2
                i += 1
                if( i < n ) : x2, y2 = regionData[i]
            changeInterpolationSubFunction = None
            if( interpolation == 6 ) : changeInterpolationSubFunction = gnd.reactionData.crossSection.chargeParticle_changeInterpolationSubFunction
            regions.append( XYs.XYs( axes_, data = regionData[:i], accuracy = ENDF_Accuracy, index = index, changeInterpolationSubFunction = changeInterpolationSubFunction ) )
            index += 1
            regionData = regionData[i:]
            for i in xrange( 1, len( regionData ) ) :
                if( regionData[i][0] != x1 ) : break
            i -= 1
            regionData = regionData[i:]
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

    skippingMTs = []
    for MT in [ 10, 27, 101, 151 ] : skippingMTs += removeGetIfPresent( MT, MTs )
    for MT in xrange( 201, 600 ) : skippingMTs += removeGetIfPresent( MT, MTs )
    for MT in xrange( 850, 875 ) : skippingMTs += removeGetIfPresent( MT, MTs )
    if( verbose ) : logFile.write( 'Skipping MTs = %s\n' % skippingMTs )

    newMTs += MTs
    newMTs += MT5
    return( newMTs )

def dull2dEdges( data_ ) :

    data = data_
    eps = FUDGE_EPS
    xy1 = data.data[0]
    i, n = 1, len( data ) - 1
    while( i < n ) :
        xy2 = data.data[i]
        if( xy1[0] == xy2[0] ) :
            xy1[0] *= ( 1. - eps )
            xy2[0] *= ( 1. + eps )
        xy1 = xy2
        i += 1
    return( data )

def getMFDataInMFList( MFs, MFData ) :

    Data = {}
    for MF in MFData :
        if( MF in MFs ) : Data[MF] = MFData[MF]
    return( Data )

def toEnergyFunctionalData( LF, moniker, unit, data ) :

    from fudge.gnd.productData import distributions

    if( data['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LF = %s, moniker = %s" % ( LF, moniker ) )

    interpolationx, interpolationy = ENDFInterpolationToGND2d( data['interpolationInfo'][0][1] )
    axes_ = axes.axes( )
    axes_[0] = axes.axis( 'energy_in', 0, 'eV', interpolation = axes.interpolationXY( interpolationx, interpolationy ) )
    axes_[1] = axes.axis( moniker, 1, unit )
    return( distributions.energy.energyFunctionalData( axes_, data['data'], accuracy = ENDF_Accuracy, moniker = moniker ) )
