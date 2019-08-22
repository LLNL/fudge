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

"""Routines for writing an ENDF file"""

from fudge.core.math.xData import axes, XYs
from fudge.core.math import fudgemath
from pqu import physicalQuantityWithUncertainty as pQU

GNDToENDFInterpolations = { axes.linearToken : 0, axes.logToken : 1 }

def floatToFunky( f ) :

    s = '%13.6e' % f
    if s[ 11 ] == '0' :
        s = s[:9] + s[10] + s[-1]
    else :
        s = '%12.5e' % f
        s = s.replace( 'e', '' )
    return s

def endfContLine( C1, C2, L1, L2, N1, N2 ) :
    # This is the basic 'control' line in ENDF
    s = '%11s%11s%11d%11d%11d%11d' % ( floatToFunky( C1 ), floatToFunky( C2 ), L1, L2, N1, N2 )
    return s

def endfContLine2( C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS ) :
    # This is the basic 'control' line in ENDF
    s = '%11s%11s%11d%11d%11d%11d%4d%2d%3d%5d' % ( floatToFunky( C1 ), floatToFunky( C2 ),
                                                   L1, L2, N1, N2, MAT, MF, MT, NS )
    return s

def endfHeadLine( ZA, AWR, L1, L2, N1, N2 ) :
    # Indicates the start of an ENDF data section
    return endfContLine( ZA, AWR, L1, L2, N1, N2 )

def endfSENDLineNumber( ) :
    # Indicates the end of an ENDF data section for one (MF, MT) pair
    return( 99999 )

def endfSENDLine( MAT, MF ) :
    # Indicates the end of an ENDF data section for one (MF, MT) pair
    return( endfContLine2( 0, 0, 0, 0, 0, 0, MAT, MF, 0, endfSENDLineNumber( ) ) )

def endfFENDLine( MAT ) :
    # Indicates the end of an ENDF data block for one MF
    return( endfContLine2( 0, 0, 0, 0, 0, 0, MAT, 0, 0, 0 ) )

def endfMENDLine( ) :
    # Indicates the end of ENDF data for one material
    return endfContLine2( 0, 0, 0, 0, 0, 0,  0, 0, 0, 0 )

def endfTENDLine( ) :
    # Indicates the end of ENDF data
    return endfContLine2( 0, 0, 0, 0, 0, 0, -1, 0, 0, 0 )

def endfComment( text, MAT, MF, MT, NS ) :
    # used for writing the introductory comments
    s = '%66s%4d%2d%3d%5d' % ( text, MAT, MF, MT, NS )
    return s

def endfDataLine( data ) :
    dData = dataListToSupportDimensionlessPQ(data)
    # makes one ENDF data line
    s = ''
    for d in dData : s += floatToFunky( d )
    return ( "%-66s" % ( s ) )

def endfDataList( data ) :
    dData = dataListToSupportDimensionlessPQ(data)
    # writes the data in ENDF format
    dataOut = []
    for i in xrange( 0, len( dData ), 6 ): dataOut.append( endfDataLine( dData[ i:i+6 ] ) )
    return( dataOut )

def endfNdDataList( nDdata, xUnit = 'eV', yUnit = '' ) :
    dnDdata = dataListToSupportDimensionlessPQ(nDdata)
    data = []
    for x, y in dnDdata : 
        if( fudgemath.isNumber( x ) ) :
            data.append( x )
        else :
            data.append( XYs.evaluateValueAsUnit( xUnit, x ) )
        if( fudgemath.isNumber( y ) ) :
            data.append( y )
        else :
            data.append( XYs.evaluateValueAsUnit( yUnit, y ) )
    return( endfDataList( data ) )

#----------------
# modified by nrp
# support dimensionless PhysicalQuantityWithUncertainty
def dataListToSupportDimensionlessPQ(data):
	dataOut = []

	for d in data:
		if isinstance(d, pQU.PhysicalQuantityWithUncertainty) and d.unit.isDimensionless():
			dataOut.append(d.value)
		elif type(d) == type([]):
			ddataOut = []

			for dd in d:
				if isinstance(dd, pQU.PhysicalQuantityWithUncertainty) and dd.unit.isDimensionless():
					ddataOut.append(dd.value)
				else:
					ddataOut.append(dd)

			dataOut.append(ddataOut)
		else:
			dataOut.append(d)

	return dataOut
#----------------

def endfInterpolationLine( interpolation ) :
    # makes one ENDF interpolation line
    s = ''
    for d in interpolation : s += "%11d" % d
    for i in xrange( len( interpolation ), 6 ) : s += "%11d" % 0
    return( s )

def endfInterpolationList( interpolation ) :
    # writes the interpolation in ENDF format
    interpolationOut = []
    for i in xrange( 0, len( interpolation ), 6 ):
        interpolationOut.append( endfInterpolationLine( interpolation[ i:i+6 ] ) )
    return( interpolationOut )

def toEndfStringList( documentation ) :
    
    EndfString = []
    docList = documentation.getLines( )
    for docEntry in docList : EndfString.append( '%-66s' % ( docEntry ) )
    return( EndfString )

def XYToENDFInterpolation( x, y ) :

    if( y == axes.flatToken ) : return( 1 )
    elif( y == axes.chargedParticleToken ) : return( 6 )
    interpolation = 2 + GNDToENDFInterpolations[x] + 2 * GNDToENDFInterpolations[y]
    return( interpolation )

def XYStringToENDFInterpolation( s ) :

    xPrefix, yPrefix = None, None
    x, y = s.split( ',' )
    if( ':' in x ) : xPrefix, x = x.split( ':' )
    if( ':' in y ) : yPrefix, y = y.split( ':' )
    if( not( xPrefix is None ) ) : raise Exception( 'xPrefix = "%s"' % xPrefix )
    if( not( yPrefix is None ) ) : raise Exception( 'yPrefix = "%s"' % yPrefix )
    return( XYToENDFInterpolation( x, y ) )

def twoAxesToENDFInterpolation( axes, indexX ) :

    independent, dependent, dummy = axes[indexX].interpolation.getInterpolationTokens( )
    return( XYToENDFInterpolation( independent, dependent ) )

def writeEndfINTG( row, col, vals, ndigit ):
    """ special INTG format, only used for sparse correlation matrices """
    def toStr(val):
        if val == 0: return ' ' * (ndigit+1)
        if abs(val) > 10**ndigit:
            raise ValueError ("%i too large for INTG format with ndigit=%i!" % (val,ndigit))
        return ("%%%ii" % (ndigit+1)) % val

    linelength = 56 # space available for writing integers
    nints = linelength // (ndigit+1)
    if ndigit==3: nints = 13    # special case

    while len(vals) < nints: vals.append(0)

    rets = "%5i%5i" % (row,col)
    # number of padding spaces depends on ndigit:
    padleft, padright = {2: (' ',' '), 3: (' ','   '), 4: (' ',''), 5: (' ',' '), 6: ('','')}[ndigit]
    rets += padleft
    for a in vals:
        rets += toStr(a)
    rets += padright
    return rets

def endfMFListToFinalFile( endfMFList, MAT, lineNumbers=True ):
	""" from dictionary of MF/MTs, build final ENDF file as a string """
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
	endfMFList[1][451][3] = endfMFList[1][451][3][:55] + '%11i' % len(directory[:-1])
	endfMFList[1][451] += directory

	endfList = [ "%66s%s" % ( " ", endfFENDLine( 1 )[66:75] ) ]
	for MF in MFs :
		MFData = endfMFList[MF]
		MTs = sorted( MFData.keys( ) )
		for MT in MTs :
			data = MFData[MT]
			for i, datum in enumerate( data ) :
				if( datum == 99999 ) :
					endfList.append( endfSENDLine( MAT, MF ) )
				else :
					if lineNumbers:
						endfList.append( '%-66s%4d%2d%3d%5d' % ( datum, MAT, MF, MT, i + 1 ) )
					else:
						endfList.append( '%-66s%4d%2d%3d' % ( datum, MAT, MF, MT ) )
		if( len( MTs ) > 0 ) : endfList.append( endfFENDLine( MAT ) )
	endfList.append( endfMENDLine( ) )
	endfList.append( endfTENDLine( ) )
	endfList.append( '' )
	return( '\n'.join( endfList ) )

