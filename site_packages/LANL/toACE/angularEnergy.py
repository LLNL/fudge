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
This module adds the method toACE to the classes in the fudge.gnd.productData.distributions.angularEnergy module.
"""

from xData import axes, XYs, W_XYs
from fudge.gnd.productData.distributions import angularEnergy


class angularFor_angularEnergy( W_XYs.W_XYs ) :

    def __init__( self, angularEnergy_ ) :

        self.productFrame = angularEnergy_.getProductFrame( )
        axes_ = axes.defaultAxes( 3 )
        axes_[0] = angularEnergy_.axes[0]
        axes_[1] = angularEnergy_.axes[1]
        axes_[2] = axes.axis( "P(mu|energy_in)", 2, "" )
        W_XYs.W_XYs.__init__( self, axes_ )

        for w_xys in angularEnergy_ :
            P_mu = [ [ xys.value, xys.integrate( ) ] for xys in w_xys ]
            self.append( XYs.XYs( XYs.XYs.defaultAxes( ), P_mu, 1e-3, value = w_xys.value ) )

    def getProductFrame( self ) :

        return( self.productFrame )

    def isIsotropic( self ) :

        return( False )

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 67, offset + len( weight ) + 4 ] + weight
    e_inFactor, e_outFactor = self.domainUnitConversionFactor( 'MeV' ), self.axes[2].unitConversionFactor( 'MeV' )

    INTE = -1
    independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTE = 1
    elif( independent == dependent == axes.linearToken ) :
        INTE = 2
    if( INTE == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for incident energy' % ( independent, dependent ) )

    INTMU = -1
    independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTMU = 1
    elif( independent == dependent == axes.linearToken ) :
        INTMU = 2
    if( INTMU == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for outgoing energy' % ( independent, dependent ) )

    INTEP = -1
    independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        INTEP = 1
    elif( independent == dependent == axes.linearToken ) :
        INTEP = 2
    if( INTEP == -1 ) : raise Exception( 'Interpolation "%s, %s" not supported for outgoing energy' % ( independent, dependent ) )

    NE, e_ins, Ls, MuData = len( self ), [], [], []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for w_xys in self :
        e_ins.append( w_xys.value * e_inFactor )
        Ls.append( offset + len( MuData ) )

        NMU, XMU, LMU, EpPData = len( w_xys ), [], [], []
        offset_LC = Ls[-1] + 1 + 2 * NMU
        for i1, xys_ in enumerate( w_xys ) :
            XMU.append( xys_.value )
            LMU.append( offset_LC + len( EpPData ) )

            xys = xys_.normalize( )
            cdfOfEps = xys.runningIntegral( )
            cdfOfEps[-1] = 1.
            Eps, pdfOfEps = [], []
            for Ep1, pdf1 in xys :
                Eps.append( Ep1 * e_outFactor )
                pdfOfEps.append( pdf1 / e_outFactor )
            EpPData += [ INTEP, len( Eps ) ] + Eps + pdfOfEps + cdfOfEps
        MuData += [ INTMU, NMU ] + XMU + LMU + EpPData

    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + MuData )

angularEnergy.pointwise.toACE = toACE
