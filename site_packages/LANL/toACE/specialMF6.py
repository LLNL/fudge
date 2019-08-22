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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
This module handles a few special case for MF 6 data that are not easily converted from the ENDF format to the ACE format.
"""

from fudge.gnd.productData.distributions import angular, energy, energyAngular
from fudge.core.math.xData import axes

class multipleNeutronDistributions :
    """
    This class simulates gnd class productData.distributions.energy.weightedFunctionals for MF 6 distributions.
    """

    def __init__( self ) :

        self.neutrons = []

    def __len__( self ) :

        return( len( self.neutrons ) )

    def __getitem__( self, index ) :

        return( self.neutrons[index] )

    def append( self, neutronData ) :

        self.neutrons.append( neutronData )

    def toACE( self, label, offset, weight, **kwargs  ) :

        DLWs = []
        for i1, functional in enumerate( self ) :
            weight = [ 2 ] + [ E1 for E1 in functional.domain( 'MeV' ) ] + [ 1., 1. ]
            DLW = functional.toACE( label, offset, weight, **kwargs )
            offset += len( DLW )
            if( functional is not self[-1] ) : DLW[0] = offset + 1
            DLWs += DLW
        return( DLWs )

class energyToEnergyAngular :

    def __init__( self, energyData ) :

        self.energyData = energyData

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.energyData.domain( unitTo = unitTo, asPQU = asPQU ) )

    def toACE( self, label, offset, weight, **kwargs  ) :

        energyData = self.energyData
        header = [ 0, 61, len( weight ) + 4 ] + weight
        e_inFactor, e_outFactor = energyData.axes[0].unitConversionFactor( 'MeV' ), energyData.axes[1].unitConversionFactor( 'MeV' )

        INTE = -1
        independent, dependent, qualifier = energyData.axes[0].interpolation.getInterpolationTokens( )
        if( dependent == axes.flatToken ) :
            INTE = 1
        elif( independent == dependent == axes.linearToken ) :
            INTE = 2
        if( INTE == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for incident energy' % ( independent, dependent ) )

        INTT = -1
        independent, dependent, qualifier = energyData.axes[1].interpolation.getInterpolationTokens( )
        if( dependent == axes.flatToken ) :
            INTT = 1
        elif( independent == dependent == axes.linearToken ) :
            INTT = 2
        if( INTT == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for outgoing energy' % ( independent, dependent ) )

        NE, e_ins, Ls, epData = len( energyData ), [], [], []
        offset += len( header ) + 3 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).

        for xys in energyData :
            e_ins.append( xys.value * e_inFactor )
            Ls.append( offset + len( epData ) )
            EOuts, pdfOfEOuts, cdfOfEOuts, LCs, muPData = [], [], [], [], []
            NP = len( xys )
            offset_LC = Ls[-1] + 1 + 4 * NP
            x1 = None
            for Ep2, P2, in xys :
                x2 = Ep2 * e_outFactor
                EOuts.append( x2 )

                y2 = P2 / e_outFactor
                pdfOfEOuts.append( y2 )

                if( x1 is None ) :
                    runningIntegral = 0
                else :
                    if( INTT == 1 ) :
                        runningIntegral += y1 * ( x2 - x1 )
                    else :
                        runningIntegral += 0.5 * ( y1 + y2 ) * ( x2 - x1 )
                cdfOfEOuts.append( runningIntegral )
                x1, y1 = x2, y2

                LCs.append( offset_LC )
                muData = [ 1, 2, -1.0, 1.0, 0.5, 0.5, 0.0, 1.0 ]
                offset_LC += len( muData )
                muPData += muData
            for i1 in range( NP ) :
                pdfOfEOuts[i1] /= cdfOfEOuts[-1]
                cdfOfEOuts[i1] /= cdfOfEOuts[-1]
            epData += [ INTT, NP ] + EOuts + pdfOfEOuts + cdfOfEOuts + LCs + muPData

        return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + epData )

def neutrons( neutronDatas ) :

    LAW61 = multipleNeutronDistributions( )
    frame = neutronDatas[0]['frame']
    for neutronData in neutronDatas :
        if( neutronData['multiplicity'] != 1 ) : raise Exception( 'multiplicity = %s != 1 is not currently supported' % neutronData['multiplicity'] )
        angularData = neutronData['angularData']
        if( neutronData['frame'] != frame ) : raise Exception( 'All neutrons must have the same frame' )
        energyData = neutronData['energyData']
        if( isinstance( energyData, energyAngular.pointwise ) ) :
            LAW61.append( energyData )
        elif( isinstance( energyData, energy.pointwise ) ) :
            if( not( isinstance( angularData, angular.isotropic ) ) ) : raise Exception( 'Only isotropic angular data is currently supported' )
            angularData = None
            LAW61.append( energyToEnergyAngular( energyData ) )
        else :
            raise Exception( "Unsupported energyData %s" % repr( energyData ) )
        if( angularData is not None ) : raise Exception( 'angularData must be None, not %s' % neutronData['angularData'] )

    return( len( LAW61 ), None, LAW61 )
