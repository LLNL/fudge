# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.distributions.energyAngularMC module.
"""

from xData import standards as standardsModule

from fudge.productData.distributions import energyAngularMC as energyAngularMCModule

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 61, offset + len( weight ) + 4 ] + weight

    energyData = self.energy.data
    energyAngularData = self.energyAngular.data

    INTE = -1
    interpolation = energyData.interpolation
    if( interpolation == standardsModule.interpolation.flatToken ) :
        INTE = 1
    elif( interpolation == standardsModule.interpolation.linlinToken ) :
        INTE = 2
    if( INTE == -1 ) : raise Exception( 'Interpolation "%s" not supported for incident energy' % interpolation )
    if( energyData.interpolationQualifier == standardsModule.interpolation.unitBaseToken ) : INTE += 20

    INTT = -1
    interpolation = energyData[0].interpolation
    if( interpolation == standardsModule.interpolation.flatToken ) :
        INTT = 1
    elif( interpolation == standardsModule.interpolation.linlinToken ) :
        INTT = 2
    if( INTT == -1 ) : raise Exception( 'Interpolation "%s" not supported for outgoing energy' % interpolation )

    NE = len( energyData )
    e_ins = []
    Ls = []
    EpData = []
    offset += len( header ) + 3 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for i1, _EpData in enumerate( energyData ) :
        e_ins.append( _EpData.outerDomainValue  )
        Ls.append( offset + len( EpData ) )

        LCs = []
        muPData = []
        NP = len( _EpData )
        offset_LC = Ls[-1] + 1 + 4 * NP + 1
        for i2, _muData in enumerate( energyAngularData[i1] ) :
            LCs.append( offset_LC )
            mus = _muData.xs.values.values
            interpolation = { standardsModule.interpolation.flatToken : 1, standardsModule.interpolation.linlinToken : 2 }[_muData.interpolation]
            muData = [ interpolation, len( mus ) ] + mus + _muData.pdf.values.values + _muData.cdf.values.values
            offset_LC += len( muData )
            muPData += muData
        EpData += [ INTT, NP ] + _EpData.xs.values.values + _EpData.pdf.values.values + _EpData.cdf.values.values + LCs + muPData

    return( header + [ 1, NE, INTE, NE ] + e_ins + Ls + EpData )

energyAngularMCModule.form.toACE = toACE
