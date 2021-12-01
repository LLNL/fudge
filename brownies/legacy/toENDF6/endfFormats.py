# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Routines for writing an ENDF file"""

from pqu import PQU as PQUModule

from xData import XYs as XYsModule
from xData import regions as regionsModule

useRedsFloatFormat = False

def floatToFunky( valueIn ) :

    value = float( valueIn )
    if( ( value != 0.0 ) and ( abs( value ) < 1e-99 ) ) :
        if( value < 0.0 ) :
            s = '%11.4e' % value
        else :
            s = '%12.5e' % value
        s = s.replace( 'e', '' )
    elif( useRedsFloatFormat ) :
        s, floatStr_Orig = floatToFunky2( value )
    else :
        s = '%13.6e' % value
        sValue = float( s )
        if s[ 11 ] == '0' :
            s = s[:9] + s[10] + s[-1]
        else :
            s = '%12.5e' % value
            sValue = float( s )
            s = s.replace( 'e', '' )
        if( abs( sValue - value ) > abs( 1e-11 * value ) ) :
            floatStr, floatStr_Orig = floatToFunky2( value )
            if( float( sValue ) != float( floatStr_Orig ) ) : s = floatStr
    return( s )

def floatToFunky2( value ) :

    floatStr = '%13.6e' % value
    floatStr_Orig = floatStr
    valueStr = '%.12g' % value
    valueStrLength = len( valueStr )
    eNotInStrValue = 'e' not in valueStr
    if( useRedsFloatFormat and eNotInStrValue ) : valueStr = valueStr[:10]
    length = 10
    if( value < 0.0 ) : length = 11                             # Allow for '-' sign.
    if( ( valueStrLength == length+1 ) and valueStr[:2] == '0.' ) :
        valueStr = valueStr[1:]     # some ENDF-VIII evaluations store numbers like '.700842459'
        valueStrLength -= 1
    if( ( valueStrLength <= length ) and eNotInStrValue ) :
        floatStr = valueStr.rjust( 11 )
        floatStr_Orig = floatStr
    elif( floatStr[11] == '0' ) :
        floatStr = floatStr[:9] + floatStr[10] + floatStr[-1]
    else :
        floatStr = '%12.5e' % value
        floatStr_Orig = floatStr
        floatStr = floatStr.replace( 'e', '' )
    return( floatStr, floatStr_Orig )


def endfContLine( C1, C2, L1, L2, N1, N2 ) :
    """This is the basic 'control' line in ENDF."""

    return( '%11s%11s%11d%11d%11d%11d' % ( floatToFunky( C1 ), floatToFunky( C2 ), L1, L2, N1, N2 ) )

def endfContLine2( C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS = None ) :
    """This is the basic 'control' line in ENDF."""

    line = '%11s%11s%11d%11d%11d%11d%4d%2d%3d' % ( floatToFunky( C1 ), floatToFunky( C2 ),
                                                   L1, L2, N1, N2, MAT, MF, MT )
    if NS is not None:
        line += '%5d' % NS
    return line

def endfHeadLine( ZA, AWR, L1, L2, N1, N2 ) :
    """Indicates the start of an ENDF data section."""

    return( endfContLine( ZA, AWR, L1, L2, N1, N2 ) )

def endfSENDLineNumber( ) :
    """Indicates the end of an ENDF data section for one (MF, MT) pair."""

    return( 99999 )

def endfSENDLine( MAT, MF, lineNumbers = True ) :
    """Indicates the end of an ENDF data section for one (MF, MT) pair."""

    lineNum = None
    if lineNumbers: lineNum = endfSENDLineNumber()
    return( endfContLine2( 0, 0, 0, 0, 0, 0, MAT, MF, 0, lineNum ) )

def endfFENDLine( MAT, lineNumbers = True ) :
    """Indicates the end of an ENDF data block for one MF."""

    lineNum = None
    if lineNumbers: lineNum = 0
    return( endfContLine2( 0, 0, 0, 0, 0, 0, MAT, 0, 0, lineNum ) )

def endfMENDLine( lineNumbers = True ) :
    """Indicates the end of ENDF data for one material."""

    lineNum = None
    if lineNumbers: lineNum = 0
    return endfContLine2( 0, 0, 0, 0, 0, 0,  0, 0, 0, lineNum )

def endfTENDLine( lineNumbers = True ) :
    """Indicates the end of ENDF data."""

    lineNum = None
    if lineNumbers: lineNum = 0
    return endfContLine2( 0, 0, 0, 0, 0, 0, -1, 0, 0, lineNum )

def endfComment( text, MAT, MF, MT, NS ) :
    """Used for writing the introductory comments."""

    return( '%66s%4d%2d%3d%5d' % ( text, MAT, MF, MT, NS ) )

def endfDataLine( data ) :
    """Makes one ENDF data line."""

    dData = dataListToSupportDimensionlessPQ( data )
    s = ''
    for d in dData : s += floatToFunky( d )
    return ( "%-66s" % ( s ) )

def endfDataList( data ) :
    """Writes the data in ENDF format."""

    dData = dataListToSupportDimensionlessPQ( data )
    dataOut = []
    for i1 in range( 0, len( dData ), 6 ): dataOut.append( endfDataLine( dData[i1:i1+6] ) )
    return( dataOut )

def endfNdDataList( nDdata, xUnit = 'eV', yUnit = '' ) :

    data = dataListToSupportDimensionlessPQ( nDdata )
    return( endfDataList( data ) )

def dataListToSupportDimensionlessPQ( data ) :

    dataOut = []
    for d1 in data :
        if( isinstance( d1, PQUModule.PQU ) and d1.unit.isDimensionless( ) ) :
            dataOut.append( float( d1 ) )
        elif( isinstance( d1, list ) ) :
            for dd in d1 : dataOut.append( float( dd ) )
        else :
            dataOut.append( d1 )
    return( dataOut )

def toTAB1( self, xUnitTo, yUnitTo, C1 = 0, C2 = 0, L1  = 0, L2 = 0 ) :

    from .gndsToENDF6 import gndsToENDFInterpolationFlag

    if( isinstance( self, XYsModule.XYs1d ) ) :
        ENDFDataList = [ endfContLine( C1, C2, 0, 0, 1, len( self ) ) ] + \
            endfInterpolationList( [ len( self ), gndsToENDFInterpolationFlag( self.interpolation ) ] )
        ENDFDataList += endfNdDataList( self, xUnit = xUnitTo, yUnit = yUnitTo )
    elif( isinstance( self, regionsModule.regions1d ) ) :
        interpolations, data = [], []
        for region in self :
            subData = region.copyDataToXYs( )
            if( len( data ) > 0 ) :
                if( subData[0] == data[-1] ) : subData.pop( 0 )
            data += subData
            interpolations += [ len( data ), gndsToENDFInterpolationFlag( region.interpolation ) ]
        NR = len( interpolations ) / 2
        ENDFDataList = [ endfContLine( C1, C2, 0, 0, NR, len( data ) ) ]
        ENDFDataList += endfInterpolationList( interpolations )
        ENDFDataList += endfNdDataList( data, xUnit = xUnitTo, yUnit = yUnitTo )
    else :
        raise 'hell - fix me'
    return( ENDFDataList )

def endfInterpolationLine( interpolation ) :
    """Makes one ENDF interpolation line."""

    s = ''
    for d in interpolation : s += "%11d" % d
    for i1 in range( len( interpolation ), 6 ) : s += "%11d" % 0
    return( s )

def endfInterpolationList( interpolation ) :
    """Writes the interpolation in ENDF format."""

    interpolationOut = []
    for i1 in range( 0, len( interpolation ), 6 ):
        interpolationOut.append( endfInterpolationLine( interpolation[i1:i1+6] ) )
    return( interpolationOut )

def toEndfStringList( documentation ) :

    EndfString = []
    docList = documentation.getLines( )
    for docEntry in docList : EndfString.append( '%-66s' % ( docEntry ) )
    return( EndfString )

def writeEndfINTG( row, col, vals, ndigit ):
    """Special INTG format, only used for sparse correlation matrices."""

    def toStr( val ) :
        if( val == 0 ): return ' ' * ( ndigit + 1 )
        if( abs( val ) > 10**ndigit ) :
            raise ValueError( "%d too large for INTG format with ndigit=%d!" % ( val, ndigit ) )
        return( ("%%%dd" % ( ndigit + 1 ) ) % val )

    linelength = 56 # space available for writing integers
    nints = linelength // (ndigit+1)
    if ndigit==3: nints = 13    # special case

    while len(vals) < nints: vals.append(0)

    rets = "%5d%5d" % (row,col)
    # number of padding spaces depends on ndigit:
    padleft, padright = {2: (' ',' '), 3: (' ','   '), 4: (' ',''), 5: (' ',' '), 6: ('','')}[ndigit]
    rets += padleft
    for a in vals:
        rets += toStr(a)
    rets += padright
    return rets

def endfMFListToFinalFile( endfMFList, MAT, lineNumbers = True ) :
    """From dictionary of MF/MTs, build final ENDF file as a string."""

    directory = []
    MFs = sorted(endfMFList.keys())
    for MF in MFs:
        MFData = endfMFList[MF]
        MTs = sorted( MFData.keys( ) )
        for MT in MTs :
            if( ( MF == 1 ) and ( MT == 451 ) ) : continue
            data = MFData[MT]
            directory.append( "%33d%11d%11d%11d" % ( MF, MT, len( data ) - 1, 0 ) )
    directory.insert( 0, "%33d%11d%11d%11d" % ( 1, 451, len( directory ) + len( endfMFList[1][451] ) + 1, 0 ) )  # Add current line of directory
    directory.append( 99999 )
    endfMFList[1][451][3] = endfMFList[1][451][3][:55] + '%11d' % len(directory[:-1])
    endfMFList[1][451] += directory

    endfList = [ "%66s%s" % ( " ", endfFENDLine( 1 )[66:75] ) ]
    for MF in MFs :
        MFData = endfMFList[MF]
        MTs = sorted( MFData.keys( ) )
        for MT in MTs :
            data = MFData[MT]
            for i1, datum in enumerate( data ) :
                if( datum == 99999 ) :
                    endfList.append( endfSENDLine( MAT, MF, lineNumbers ) )
                else :
                    if lineNumbers:
                        endfList.append( '%-66s%4d%2d%3d%5d' % ( datum, MAT, MF, MT, i1 + 1 ) )
                    else:
                        endfList.append( '%-66s%4d%2d%3d' % ( datum, MAT, MF, MT ) )
        if( len( MTs ) > 0 ) : endfList.append( endfFENDLine( MAT, lineNumbers ) )
    endfList.append( endfMENDLine( lineNumbers ) )
    endfList.append( endfTENDLine( lineNumbers ) )
    endfList.append( '' )
    return( '\n'.join( endfList ) )
