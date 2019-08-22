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
This module adds the method toACE to the classes in the fudge.gnd.productData.distributions.energyAngular module.
"""

from fudge.core.math.xData import axes
from fudge.gnd.productData.distributions import energyAngular

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 44, offset + len( weight ) + 4 ] + weight
    e_inFactor, e_outFactor = self.domainUnitConversionFactor( 'MeV' ), self.axes[1].unitConversionFactor( 'MeV' )

    independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
    if( independent != dependent != axes.linearToken ) : raise Exception( 'interpolation = %s and %s not supported' % ( independent, dependent ) )

    NE, e_ins, Ls, epData = len( self ), [], [], []
    offset += len( header ) + 1 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for energy_coefficients in self :
        pdf, Rs, As = self.getFRAatEnergy_asLinearPointwise( energy_coefficients.value )
        pdf_, Rs, As = pdf.commonXGrid( [ Rs, As ] )
        e_ins.append( energy_coefficients.value * e_inFactor )
        Ls.append( offset + len( epData ) )
        pdf_ = pdf_.normalize( )
        cdf = pdf_.runningIntegral( )
        eps, pdf = [], []
        for x1, y1 in pdf_ :
            eps.append( e_outFactor * x1 )
            pdf.append( y1 / e_outFactor )
        Rs = [ r for x, r in Rs ]
        As = [ a for x, a in As ]
        cdf[-1] = 1.
        epData += [ 2, len( eps ) ] + eps + pdf + cdf + Rs + As
    return( header + [ 0, NE ] + e_ins + Ls + epData )

energyAngular.KalbachMann.toACE = toACE

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 61, offset + len( weight ) + 4 ] + weight
    e_inFactor, e_outFactor = self.domainUnitConversionFactor( 'MeV' ), self.axes[1].unitConversionFactor( 'MeV' )

    INTE = -1
    independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTE = 1
    elif( independent == dependent == axes.linearToken ) :
        INTE = 2
    if( INTE == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for incident energy' % ( independent, dependent ) )

    INTT = -1
    independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTT = 1
    elif( independent == dependent == axes.linearToken ) :
        INTT = 2
    if( INTT == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for outgoing energy' % ( independent, dependent ) )

    NE, e_ins, Ls, epData = len( self ), [], [], []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for w_xys in self :
        e_ins.append( w_xys.value * e_inFactor )
        Ls.append( offset + len( epData ) )
        EOuts, pdfOfEOuts, cdfOfEOuts, LCs, muPData = [], [], [], [], []
        NP = len( w_xys )
        offset_LC = Ls[-1] + 1 + 4 * NP
        for i1, xys_ in enumerate( w_xys ) :
            x2 = xys_.value * e_outFactor
            EOuts.append( x2 )

            norm = xys_.integrate( )
            y2 = norm / e_outFactor
            pdfOfEOuts.append( y2 )

            if( i1 == 0 ) :
                runningIntegral = 0
            else :
                if( INTT == 1 ) :
                    runningIntegral += y1 * ( x2 - x1 )
                else :
                    runningIntegral += 0.5 * ( y1 + y2 ) * ( x2 - x1 )
            cdfOfEOuts.append( runningIntegral )
            x1, y1 = x2, y2

            LCs.append( offset_LC )
            if( norm == 0 ) :
                muData = [ 1, 2, -1.0, 1.0, 0.5, 0.5, 0.0, 1.0 ]
            else :
                xys = xys_.normalize( )
                cdfOfMus = xys.runningIntegral( )
                mus, pdfOfMus = [], []
                for mu1, pdf1 in xys :
                    mus.append( mu1 )
                    pdfOfMus.append( pdf1 )
                cdfOfMus[-1] = 1.
                muData = [ 1, len( mus ) ] + mus + pdfOfMus + cdfOfMus
            offset_LC += len( muData )
            muPData += muData
        for i1 in range( NP ) :
            pdfOfEOuts[i1] /= cdfOfEOuts[-1]
            cdfOfEOuts[i1] /= cdfOfEOuts[-1]
        epData += [ INTT, NP ] + EOuts + pdfOfEOuts + cdfOfEOuts + LCs + muPData

    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + epData )

energyAngular.pointwise.toACE = toACE
