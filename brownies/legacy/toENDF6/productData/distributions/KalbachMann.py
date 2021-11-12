# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import XYs as XYsModule

from fudge.productData.distributions import KalbachMann as KalbachMannModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

#
# KalbachMann
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    fSubform = self.fSubform.data
    rSubform = self.rSubform.data
    aSubform = self.aSubform
    if( not( aSubform.isEmptyASubform( ) ) ) : raise 'hell - FIXME'

    outgoingInterpolation = set( [val.interpolation for val in fSubform] + [val.interpolation for val in rSubform] )
    if len(outgoingInterpolation) != 1:
        raise NotImplementedError("Only one outgoing interpolation supported when writing Kalbach-Mann to ENDF-6")

    LEP = { 'flat': 1, 'lin-lin': 2 }[outgoingInterpolation.pop()]
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 2, LEP, 1, len( fSubform ) ) ]
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( fSubform ), 2 ] )
    EInUnit = fSubform.axes[2].unit
    EpUnit = fSubform.axes[1].unit
    fUnit = fSubform.axes[0].unit
    EInFactor = PQUModule.PQU( 1, EInUnit ).getValueAs( 'eV' )
    EpFactor = PQUModule.PQU( 1, EpUnit ).getValueAs( 'eV' )
    fFactor = PQUModule.PQU( 1, fUnit ).getValueAs( '1/eV' )
    for i1, fAtEnergy in enumerate( fSubform ) :
        rAtEnergy = rSubform[i1]
        if( not( isinstance( fAtEnergy, XYsModule.XYs1d ) ) ) : raise 'hell - FIXME'
        if( not( isinstance( rAtEnergy, XYsModule.XYs1d ) ) ) : raise 'hell - FIXME'
        value = fAtEnergy.outerDomainValue
        if( value != rAtEnergy.outerDomainValue ) : raise 'hell - FIXME'
        value *= EInFactor
        coefficients = []
        for i1, ( Ep1, f1 ) in enumerate( fAtEnergy ) :
            Ep2, r1 = rAtEnergy[i1]
            if( Ep1 != Ep2 ) : raise 'hell - FIXME'
            coefficients += [ EpFactor * Ep1, fFactor * f1, r1 ]
        length = len( coefficients )
        ENDFDataList += [ endfFormatsModule.endfContLine( 0, value, 0, 1, length, length / 3 ) ]
        ENDFDataList += endfFormatsModule.endfDataList( coefficients )
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 1, self.productFrame, ENDFDataList )

KalbachMannModule.form.toENDF6 = toENDF6
