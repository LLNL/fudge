# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from fudge import abstractClasses as abstractClassesModule

from .photonScattering import coherent as coherentModule
from .photonScattering import incoherent as incoherentModule
from .chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from .thermalNeutronScatteringLaw import coherentElastic as TNSL_coherentElasticModule
from .thermalNeutronScatteringLaw import incoherentElastic as TNSL_incoherentElasticModule
from .thermalNeutronScatteringLaw import incoherentInelastic as TNSL_incoherentInelasticModule

#
# Double differential cross section component.
#
class Component( abstractClassesModule.Component ) :

    moniker = 'doubleDifferentialCrossSection'

    def __init__( self ) :

        abstractClassesModule.Component.__init__( self, ( coherentModule.Form, incoherentModule.Form, CoulombPlusNuclearElasticModule.Form,
                TNSL_coherentElasticModule.Form, TNSL_incoherentElasticModule.Form, TNSL_incoherentInelasticModule.Form ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainMin( self ) :

        return( self.evaluated.domainMin )

    @property
    def domainMax( self ) :

        return( self.evaluated.domainMax )

    @property
    def domainUnit( self ) :

        return( self.evaluated.domainUnit )

    def fixDomains(self, labels, energyMin, energyMax):
        """This method does nothing."""

        return 0

    def listOfProducts(self):
        """Returns, as a set, the list of PoP's ids for all products (i.e., outgoing particles) for *self*."""

        products = set()
        if len(self) > 0: products.add(self[0].pid)

        return products

