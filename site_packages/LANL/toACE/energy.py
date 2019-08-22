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

"""
This module adds the method toACE to the classes in the fudge.gnd.productData.distributions.energy module.
"""

from xData import axes
from fudge.gnd.productData.distributions import energy

#
#   pointwise energy (i.e., f(E'|E)) logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 4, offset + len( weight ) + 4 ] + weight
    e_inFactor, e_outFactor = self.domainUnitConversionFactor( 'MeV' ), self[0].domainUnitConversionFactor( 'MeV' )
    NE = len( self )

    INTE = -1
    independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTE = 1
    elif( independent == dependent == axes.linearToken ) :
        INTE = 2
    if( INTE == -1 ) :
        INTE = 2
        print '    WARNING: for %s changing interpolation from "%s,%s" to "%s,%s"' % ( label, independent, dependent, axes.linearToken, axes.linearToken )

    INTT = -1
    independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTT = 1
    elif( independent == dependent == axes.linearToken ) :
        INTT = 2

    distribution = self
    if( INTT == -1 ) : distribution = self.toPointwise_withLinearXYs( )

    e_ins, Ls, epData = [], [], []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header plus NR, NE, Es, Ls, (1-based).
    for xys_ in distribution :
        e_ins.append( xys_.value * e_inFactor )
        Ls.append( offset + len( epData ) )
        xys = xys_.normalize( )
        cdf = xys.runningIntegral( )
        eps, pdf = [], []
        for x1, y1 in xys :
            eps.append( e_outFactor * x1 )
            pdf.append( y1 / e_outFactor )
        cdf[-1] = 1.
        epData += [ INTT, len( eps ) ] + eps + pdf + cdf
    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + epData )

energy.pointwise.toACE = toACE

#
#   NBodyPhaseSpace logic
#
def toACE( self, label, offset, weight, massUnit = None, neutronMass = None, **kwargs ) :

    return( [ 0, 66, offset + len( weight ) + 4 ] + weight + [ self.numberOfProducts, self.numberOfProductsMasses.getValueAs( massUnit ) / neutronMass  ] )

energy.NBodyPhaseSpace.toACE = toACE

#
#   evaporationSpectrum and simpleMaxwellianFissionSpectrum logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, self.LF, offset + len( weight ) + 4 ] + weight
    theta = self.parameter1.toPointwise_withLinearXYs( )
    e_inFactor, e_outFactor = theta.axes[0].unitConversionFactor( 'MeV' ), theta.axes[1].unitConversionFactor( 'MeV' )

    NE, e_ins, Ts = len( theta ), [], []
    for x1, y1 in theta :
        e_ins.append( x1 * e_inFactor )
        Ts.append( y1 * e_outFactor )
    return( header + [ 0, NE ] + e_ins + Ts + [ self.U.getValueAs( 'MeV' ) ] )

energy.evaporationSpectrum.toACE = toACE
energy.simpleMaxwellianFissionSpectrum.toACE = toACE

#
#   weightedFunctionals logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    DLWs = []
    for functional in self :
        weightEs = []
        weightPs = []
        e_inFactor = functional.weight.axes[0].unitConversionFactor( 'MeV' )
        for E1, P1 in functional.weight :
            weightEs.append( e_inFactor * E1 )
            weightPs.append( P1 )
        weight = [ len( weightEs ) ] + weightEs + weightPs
        DLW = functional.functional.toACE( label, offset, weight, **kwargs )
        offset += len( DLW )
        if( functional is not self[-1] ) : DLW[0] = offset + 1
        DLWs += DLW

    return( DLWs )        

energy.weightedFunctionals.toACE = toACE

#
#   weightedFunctionals logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    linear = self.toPointwise_withLinearXYs( 1e-3, 1e-6, 1e-6 )
    return( linear.toACE( label, offset, weight, **kwargs ) )

energy.semiPiecewise.toACE = toACE

#
#   Watt spectrum logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, self.LF, offset + len( weight ) + 4 ] + weight

    e_inFactor, aFactor = self.parameter1.axes[0].unitConversionFactor( 'MeV' ), self.parameter1.axes[1].unitConversionFactor( 'MeV' )
    A = self.parameter1.toPointwise_withLinearXYs( )
    data = [ 0, len( A ) ]
    for E1, a1 in A : data += [ e_inFactor * E1, aFactor * a1 ]

    e_inFactor, bFactor = self.parameter2.axes[0].unitConversionFactor( 'MeV' ), self.parameter2.axes[1].unitConversionFactor( '1/MeV' )
    B = self.parameter2.toPointwise_withLinearXYs( )
    data += [ 0, len( B ) ]
    for E1, b1 in B : data += [ e_inFactor * E1, bFactor * b1 ]
    data.append( self.U.getValueAs( 'MeV' ) )

    return( header + data )

energy.WattSpectrum.toACE = toACE

#
#   Watt spectrum logic
#
def toACE( self, label, offset, weight, **kwargs ) :

    return( self.toPointwise_withLinearXYs( ).toACE( label, offset, weight, **kwargs ) )

energy.MadlandNix.toACE = toACE
