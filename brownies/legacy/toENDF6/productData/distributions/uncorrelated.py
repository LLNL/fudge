# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from xData import axes as axesModule
from xData import regions as regionsModule
from xData import standards as standardsModule

import fudge.productData.distributions.uncorrelated as uncorrelatedModule
import fudge.productData.distributions.energy as energyModule
import fudge.productData.distributions.angular as angularModule
import fudge.productData.distributions.energyAngular as energyAngularModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :
    """
    In ENDF MF=6, some distributions should really be treated as uncorrelated: NBodyPhaseSpace, and also 
    Legendre expansions when only L=0 is listed.
    For GNDS we split these into uncorrelated angular (isotropic) and energy distributions.
    Must put back in original format when writing back to ENDF.
    """

    frame = self.productFrame
    energySubform = self.energySubform.data
    angularSubform = self.angularSubform.data
    if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
        energySubform.toENDF6( MT, endfMFList, flags, targetInfo )
    elif( targetInfo['doMF4AsMF6'] ) :
        if( isinstance( energySubform, ( energyModule.discreteGamma, energyModule.primaryGamma ) ) ) :
# BRB - this needs to be checked.
            if( targetInfo['product'].id != IDsPoPsModule.photon ) : raise ValueError( 'This logic is only for discrete gammas' )
            energyForm = energyModule.form( self.label, frame, energySubform )
            angularForm = angularModule.form( self.label, frame, angularSubform )
            energyForm.toENDF6( MT, endfMFList, flags, targetInfo )
            angularForm.toENDF6( MT, endfMFList, flags, targetInfo )
        elif( isinstance( angularSubform, angularModule.isotropic2d ) ) :                # Change to energyAngular with Legendre
            if( not( isinstance( energySubform, energyModule.XYs2d ) ) ) : raise 'hell - fix me'
            axes = axesModule.axes( rank = 4 )
            axes[3] = axesModule.axis( 'energy_in', 3, 'eV' )
            axes[2] = axesModule.axis( 'energy_out', 2, 'eV' )
            axes[1] = axesModule.axis( 'l', 1, '' )
            axes[0] = axesModule.axis( 'C_l(energy_out|energy_in)', 0, '1/eV' )
            energyAngularSubform = energyAngularModule.XYs3d( axes = axes, interpolation = energySubform.interpolation,
                    interpolationQualifier = energySubform.interpolationQualifier )
            EInFactor = PQUModule.PQU( 1, energySubform.axes[2].unit ).getValueAs( 'eV' )

            for EIn in energySubform :
                if isinstance( EIn, regionsModule.regions1d ):
                    # writing to MF6 LAW=1, LANG=1 doesn't support multiple regions, must recombine.
                    if( len( set( [ ein.interpolation for ein in EIn ] ) ) != 1 ) :
                        raise NotImplemented( "ENDF MF6 LAW=1 LANG=1 doesn't support multiple E' interpolations!" )
                    xyvals = EIn[0].copyDataToXYs()
                    for region in EIn[1:]:
                        xynew = region.copyDataToXYs()
                        xynew[0][0] *= 1.0000000001
                        xyvals.extend( xynew )
                    EIn_copy = EIn[0].copy( )
                    EIn_copy.outerDomainValue = EIn.outerDomainValue
                    EIn_copy.axes = EIn.axes.copy( )
                    EIn_copy.setData( xyvals )
                    EIn = EIn_copy
                multiD_2d = energyAngularModule.XYs2d( outerDomainValue = EIn.outerDomainValue * EInFactor, interpolation = EIn.interpolation )
                EpCls = EIn.copyDataToXYs( )
                for e_out, Cls in EpCls :
                    multiD_2d.append( energyAngularModule.Legendre( [ Cls ], outerDomainValue = e_out ) )
                energyAngularSubform.append( multiD_2d )
            form = energyAngularModule.form( '', self.productFrame, energyAngularSubform )
            if( targetInfo.get( "gammaToENDF6" ) ) : return( form )
            form.toENDF6( MT, endfMFList, flags, targetInfo )
        elif( isinstance( energySubform, energyModule.XYs2d ) and isinstance( angularSubform, angularModule.XYs2d ) ) :

            LANG, LEP = 12, 2
            if( energySubform.interpolation == standardsModule.interpolation.flatToken ) :
                LEP = 1  # interpolation for E_out
                LANG = 11
            elif( energySubform.interpolation in (standardsModule.interpolation.loglinToken,
                standardsModule.interpolation.loglogToken ) ) :
                LANG = 14
            MF6 = [ endfFormatsModule.endfContLine( 0, 0, LANG, LEP, 1, len( energySubform ) ) ]
            EInInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( energySubform.interpolation )
            MF6 +=  endfFormatsModule.endfInterpolationList( [ len( energySubform ), EInInterpolation ] )
            for indexE, EEpP in enumerate( energySubform ) :
                EMuP = angularSubform[indexE]
                if( EEpP.outerDomainValue != EMuP.outerDomainValue ) : raise Exception( "EEpP.outerDomainValue = %s != EMuP.outerDomainValue = %s" % ( EEpP.outerDomainValue, EMuP.outerDomainValue ) )
                NA, NEP = 2 * len( EMuP ), len( EEpP )
                MF6.append( endfFormatsModule.endfContLine( 0, EEpP.outerDomainValue, 0, NA, NEP * ( NA + 2 ), NEP ) )
                data = []
                for EpP in EEpP :
                    data = [ EpP[0], EpP[1] ]
                    for muP in EMuP : data += muP
                    MF6 += endfFormatsModule.endfDataList( data )
            LAW = 1
            gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
        else :
            raise Exception( 'uncorrelated.toENDF6 not supported for energy subform = %s and angular subform = %s' %
                ( energySubform.moniker, angularSubform.moniker ) )
    else :                          # original data is in uncorrelated form
        if( MT not in [ 527, 528 ] ) :
            angularForm = angularModule.form( "", frame, self.angularSubform.data.copy() )
            angularForm.toENDF6( MT, endfMFList, flags, targetInfo )
        energyForm = energyModule.form( "", frame, self.energySubform.data.copy() )
        energyForm.toENDF6( MT, endfMFList, flags, targetInfo )
        if( MT == 527 ) : endfMFList[26][MT][0]     # FIXME statement has no effect

uncorrelatedModule.form.toENDF6 = toENDF6
