# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from fudge.productData.distributions import LLNL_angularEnergy as LLNL_angularEnergyModule

from ... import gndsToENDF6 as gndsToENDF6Module


#
# LLNLAngularEnergyForm
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :    # FIXME appears to be broken code

    angularForm = self.angularForm
    angularEnergyForm = self.angularEnergyForm

    energy_inInterpolation, energy_inFunctionInterpolation, energy_inInterpolationQualifier = angularEnergyForm.axes[0].interpolation.getInterpolationTokens( )
    muInterpolation, muFunctionInterpolation, muQualifier = angularEnergyForm.axes[1].interpolation.getInterpolationTokens( )
    energy_outInterpolation, probabilityInterpolation, energy_outQualifier = angularEnergyForm.axes[2].interpolation.getInterpolationTokens( )
    frame = angularEnergyForm.getProductFrame( )
    axes = pointwise.defaultAxes( energyInterpolation = energy_inInterpolation, energyFunctionInterpolation = energy_inFunctionInterpolation, 
            energyInterpolationQualifier = energy_inInterpolationQualifier, muInterpolation = muInterpolation, 
            energy_outInterpolation = energy_outInterpolation, probabilityInterpolation = probabilityInterpolation )
    E_inRatio = PQUModule.PQU( 1, angularEnergyForm.axes[0].getUnit( ) ).getValueAs( 'eV' )
    E_outRatio = PQUModule.PQU( 1, angularEnergyForm.axes[2].getUnit( ) ).getValueAs( 'eV' )
    LAW7 = pointwise( axes, self.getProductFrame( ) )

    if( len( angularForm ) != len( angularEnergyForm ) ) :
        raise Exception( "len( angularForm ) = %s != len( angularEnergyForm ) = %s" % ( len( angularForm ), len( angularEnergyForm ) ) )
    for indexE, EMuP in enumerate( angularForm ) :
        EMuEpP = angularEnergyForm[indexE]
        if( EMuP.outerDomainValue != EMuEpP.outerDomainValue ) : raise Exception( "At indexE = %d, EMuP.outerDomainValue %s != EMuEpP.outerDomainValue = %s" % ( indexE, EMuP.outerDomainValue, EMuEpP.outerDomainValue ) )
        if( len( EMuP ) != len( EMuEpP ) ) :
            raise Exception( "At indexE = %d (E_in = %s), len( EMuP ) %s != len( EMuEpP ) = %s" % ( indexE, EMuP.outerDomainValue, len( EMuP ), len( EMuEpP ) ) )
        w_xys = W_XYs.W_XYs( axesW_XY, index = indexE, outerDomainValue = EMuP.outerDomainValue )
        for indexMu, muP in enumerate( EMuP ) :
            muEpP = EMuEpP[indexMu]
            if( muP[0] != muEpP.outerDomainValue ) : raise Exception( "At indexE = %d, mu = %s != muEpP.outerDomainValue = %s" % ( indexE, muP[0], muEpP.outerDomainValue ) )
            xys = [ [ E_outRatio * Ep, muP[1] * P / E_outRatio ] for Ep, P in muEpP ]
            xys = XYs.XYs( axesXY, xys, accuracy = muEpP.getAccuracy( ), outerDomainValue = muP[0], index = indexMu, parent = w_xys )
            w_xys.append( XYs.XYs( axesXY, muEpP * muP[1], accuracy = muEpP.getAccuracy( ), outerDomainValue = muP[0], index = indexMu, parent = w_xys ) )
        LAW7.append( w_xys )

    LAW, frame, MF6 = LAW7.toENDF6( flags, targetInfo )
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )

LLNL_angularEnergyModule.LLNLAngularEnergyForm.toENDF6 = toENDF6
