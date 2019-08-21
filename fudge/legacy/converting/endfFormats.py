# <<BEGIN-copyright>>
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

