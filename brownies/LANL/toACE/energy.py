# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.distributions.energy module.
"""

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule

from fudge.productData.distributions import energy as energyModule

def ACEInterpolation( function, label ) :

    interpolation = function.interpolation
    INTE = -1
    if interpolation == xDataEnumsModule.Interpolation.flat:
        INTE = 1
    elif interpolation == xDataEnumsModule.Interpolation.linlin:
        INTE = 2
    if( INTE == -1 ) :
        INTE = 2
        print('    WARNING: for %s changing interpolation from "%s" to "%s"' %
              (label, interpolation, xDataEnumsModule.Interpolation.linlin))

    return( INTE )
#
#   XYs2d energy (i.e., f(E'|E)) logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 4, offset + len( weight ) + 4 ] + weight

    INTE = ACEInterpolation( self, label )
    INTT = ACEInterpolation( self[0], label )

    NE = len( self )
    e_ins = []
    Ls = []
    epData = []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header plus NR, NE, Es, Ls, (1-based).
    for xs_pdf_cdf in self :
        e_ins.append( xs_pdf_cdf.outerDomainValue )
        Ls.append( offset + len( epData ) )
        eps = xs_pdf_cdf.xs.values.values
        pdf = xs_pdf_cdf.pdf.values.values
        cdf = xs_pdf_cdf.cdf.values.values
        epData += [ INTT, len( eps ) ] + eps + pdf + cdf

    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + epData )

energyModule.XYs2d.toACE = toACE
#
#   XYs2d energy (i.e., f(E'|E)) logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    XYs2d = energyModule.XYs2d( )
    for i1, region in enumerate( self ) :
        for i2, xs_pdf_cdf1d in enumerate( region ) :
            xs_pdf_cdf1d = xs_pdf_cdf1d.copy( )
            if( ( i1 > 0 ) and( i2 == 0 ) ) : xs_pdf_cdf1d.outerDomainValue *= ( 1.0 + 1e-8 )
            XYs2d.append( xs_pdf_cdf1d )
    return( XYs2d.toACE( label, offset, weight, **kwargs ) )

energyModule.Regions2d.toACE = toACE
#
#   NBodyPhaseSpace logic
#
def toACE( self, label, offset, weight, massUnit = None, neutronMass = None, **kwargs ) :

    return( [ 0, 66, offset + len( weight ) + 4 ] + weight + [ self.numberOfProducts, self.mass.getValueAs( massUnit ) / neutronMass  ] )

energyModule.NBodyPhaseSpace.toACE = toACE

#
#   Evaporation and SimpleMaxwellian logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, self.LF, offset + len( weight ) + 4 ] + weight

    theta = self.parameter1.data.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-6 )

    NE, e_ins, Ts = len( theta ), [], []
    for x1, y1 in theta :
        e_ins.append( x1 )
        Ts.append( y1  )

    return( header + [ 0, NE ] + e_ins + Ts + [ self.U.getValueAs( 'MeV' ) ] )

energyModule.Evaporation.toACE = toACE
energyModule.SimpleMaxwellianFission.toACE = toACE

#
#   WeightedFunctionals logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    DLWs = []
    for functional in self :
        weightEs = []
        weightPs = []
        for E1, P1 in functional.weight :
            weightEs.append( E1 )
            weightPs.append( P1 )
        weight = [ 0, len( weightEs ) ] + weightEs + weightPs
        DLW = functional.functional.toACE( label, offset, weight, **kwargs )
        offset += len( DLW )
        if( functional is not self[-1] ) : DLW[0] = offset + 1
        DLWs += DLW

    return( DLWs )        

energyModule.WeightedFunctionals.toACE = toACE

#
#   Watt spectrum logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, self.LF, offset + len( weight ) + 4 ] + weight

    A = self.parameter1.data
    A = A.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-6 )
    xs, ys = A.copyDataToXsAndYs( )
    data = [ 0, len( A ) ] + xs + ys

    B = self.parameter2.data
    B = B.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-6 )
    xs, ys = B.copyDataToXsAndYs( )
    data += [ 0, len( B ) ] + xs + ys
    data.append( self.U.getValueAs( 'MeV' ) )

    return( header + data )

energyModule.Watt.toACE = toACE

#
#   MadlandNix spectrum logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    return( self.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-6 ).toACE( label, offset, weight, **kwargs ) )

energyModule.MadlandNix.toACE = toACE

#
#   GeneralEvaporation logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, self.LF, offset + len( weight ) + 4 ] + weight

    theta = self.parameter1.data
    if( not( isinstance( theta, XYs1dModule.XYs1d ) ) ) : raise TypeError( 'unsupported theta form: moniker = "%s".' % theta.moniker )

    function = self.parameter2.data
    xys1d = energyModule.XYs1d( [ function.xs.values.values, function.pdf.values.values ], dataForm = 'xsandys',
            interpolation = function.interpolation )

    xys2d = energyModule.XYs2d( )
    for e_in, T in theta :
        EpP = [ [ e_out * T, P / T ] for e_out, P in xys1d ]
        xys2d.append( energyModule.XYs1d( EpP, outerDomainValue = e_in, interpolation = xys1d.interpolation ) )
    xs_pdf_cdf1d = xys2d.to_xs_pdf_cdf1d( None, None, None )
    return( xs_pdf_cdf1d.toACE( label, offset, weight, **kwargs ) )

energyModule.GeneralEvaporation.toACE = toACE
