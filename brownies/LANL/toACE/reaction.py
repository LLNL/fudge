# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the reaction class.
"""

from PoPs import IDs as IDsPoPsModule

from fudge.reactions import reaction as reactionModule

from . import product as productModule

def toACE( self, styleName, cdf_style, ACE_data, delayedNeutronRateAndDatas, verbose ) :

    MT = self.ENDF_MT
    if( verbose > 1 ) : print( '   %s: MT = %d' % ( str( self ), MT ) )

    MTData = {}
    MTData['isFission'] = self.isFission( )
    MTData['ESZ'] = self.crossSection[styleName]
    MTData['Q_final'] = self.getQ( 'MeV' )
    MTData['Q_initial'] = self.getQ( 'MeV', final = False )

    MTData[IDsPoPsModule.neutron] = []
    MTData[IDsPoPsModule.photon] = []

    if( MTData['isFission'] ) :
        for delayedNeutron in self.outputChannel.fissionFragmentData.delayedNeutrons :
                delayedNeutronRate = 1e-8 * delayedNeutron.rate[0].pqu( ).getValueAs( '1/s' )
                productModule.toACE( delayedNeutron.product, cdf_style, MTData, MT, 0 )
                delayedNeutronRateAndDatas.append( [ delayedNeutronRate, MTData[IDsPoPsModule.neutron].pop( ) ] )

        fissionEnergyReleases = self.outputChannel.fissionFragmentData.fissionEnergyReleases
        if( len( fissionEnergyReleases ) > 0 ) :
            MTData['Q_final'] = fissionEnergyReleases[0].nonNeutrinoEnergy.data.evaluate( 0.0 )
            MTData['Q_initial'] = MTData['Q_final']

    self.outputChannel.toACE( cdf_style, MTData, MT, verbose )
    ACE_data.append( ( MT, MTData ) )

reactionModule.Reaction.toACE = toACE
