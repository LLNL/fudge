# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module handles a few special case for MF 6 data that are not easily converted from the ENDF format to the ACE format.
"""

from fudge.productData.distributions import angular as angularModule


class multipleNeutronDistributions :
    """
    This class simulates gnds class productData.distributions.energy.weightedFunctionals for MF 6 distributions.
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
            domainMin, domainMax = functional.domainMin, functional.domainMax
            weight = [ 0, 2, domainMin, domainMax, 1.0, 1.0 ]
            DLW = functional.toACE( label, offset, weight, **kwargs )
            offset += len( DLW )
            if( functional is not self[-1] ) : DLW[0] = offset + 1
            DLWs += DLW
        return( DLWs )

def neutrons( neutronDatas ) :

    weightedFunctionals = multipleNeutronDistributions( )
    frame = neutronDatas[0]['frame']
    for neutronData in neutronDatas :
        if( neutronData['multiplicity'] != 1 ) : raise Exception( 'multiplicity = %s != 1 is not currently supported' % neutronData['multiplicity'] )
        angularData = neutronData['angularData']
        if( neutronData['frame'] != frame ) : raise Exception( 'All neutrons must have the same frame' )
        energyData = neutronData['energyData']
        if( hasattr( energyData, 'toACE' ) ) :
            weightedFunctionals.append( energyData )
            if( angularData is not None ) :
                if( not( isinstance( angularData, angularModule.isotropic2d ) ) ) :
                    raise Exception( 'Only isotropic angular data is currently supported' )
            angularData = None
        else :
            raise Exception( "Unsupported energyData %s" % repr( energyData ) )
        if( angularData is not None ) : raise Exception( 'angularData must be None, not %s' % neutronData['angularData'] )

    return( len( weightedFunctionals ), None, weightedFunctionals )
