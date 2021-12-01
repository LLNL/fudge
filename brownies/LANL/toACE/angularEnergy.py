# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.distributions.angularEnergyMC module.
"""

from xData import standards as standardsModule

from fudge.productData.distributions import angularEnergyMC as angularEnergyMCModule

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 67, offset + len( weight ) + 4 ] + weight

    INTE = -1
    interpolation = self.interpolation
    if( interpolation == standardsModule.interpolation.flatToken ) :
        INTE = 1
    elif( interpolation == standardsModule.interpolation.linlinToken ) :
        INTE = 2
    if( INTE == -1 ) : raise Exception( 'Interpolation "%s" not supported for incident energy' % interpolation )

    INTMU = -1
    interpolation = self[0].interpolation
    if( interpolation == standardsModule.interpolation.flatToken ) :
        INTMU = 1
    elif( interpolation == standardsModule.interpolation.linlinToken ) :
        INTMU = 2
    if( INTMU == -1 ) : raise Exception( 'Interpolation "%s" not supported for outgoing energy' % interpolation )

    INTEP = -1
    interpolation = self[0][0].interpolation
    if( interpolation == standardsModule.interpolation.flatToken ) :
        INTEP = 1
    elif( interpolation == standardsModule.interpolation.linlinToken ) :
        INTEP = 2
    if( INTEP == -1 ) : raise Exception( 'Interpolation "%s" not supported for outgoing energy' % interpolation )

    NE = len( self )
    e_ins = []
    Ls = []
    MuData = []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for w_xys in self :
        e_ins.append( w_xys.outerDomainValue )
        Ls.append( offset + len( MuData ) )

        NMU = len( w_xys )
        XMU = []
        LMU = []
        EpPData = []
        offset_LC = Ls[-1] + 1 + 2 * NMU + 1
        for i1, xys in enumerate( w_xys ) :
            XMU.append( xys.outerDomainValue )
            LMU.append( offset_LC + len( EpPData ) )
            Eps = xys.xs.values.values
            EpPData += [ INTEP, len( Eps ) ] + Eps + xys.pdf.values.values + xys.cdf.values.values
        MuData += [ INTMU, NMU ] + XMU + LMU + EpPData

    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + MuData )

angularEnergyMCModule.XYs3d.toACE = toACE
