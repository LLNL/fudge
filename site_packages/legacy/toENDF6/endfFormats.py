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

"""Routines for writing an ENDF file"""

from pqu import PQU as PQUModule

from xData import XYs as XYsModule
from xData import regions as regionsModule

from . import gndsToENDF6 as gndsToENDF6Module

useRedsFloatFormat = False

def floatToFunky( value ) :

    if( useRedsFloatFormat ) :
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
    valueStr = str( value )
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

    if( isinstance( self, XYsModule.XYs1d ) ) :
        ENDFDataList = [ endfContLine( C1, C2, 0, 0, 1, len( self ) ) ] + \
            endfInterpolationList( [ len( self ), gndsToENDF6Module.gndsToENDFInterpolationFlag( self.interpolation ) ] )
        ENDFDataList += endfNdDataList( self, xUnit = xUnitTo, yUnit = yUnitTo )
    elif( isinstance( self, regionsModule.regions1d ) ) :
        interpolations, data = [], []
        for region in self :
            subData = region.copyDataToXYs( )
            if( len( data ) > 0 ) :
                if( subData[0] == data[-1] ) : subData.pop( 0 )
            data += subData
            interpolations += [ len( data ), gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation ) ]
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
