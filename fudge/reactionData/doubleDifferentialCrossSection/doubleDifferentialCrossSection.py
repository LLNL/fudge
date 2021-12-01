# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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
class component( abstractClassesModule.component ) :

    moniker = 'doubleDifferentialCrossSection'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( coherentModule.form, incoherentModule.form, CoulombPlusNuclearElasticModule.form,
                TNSL_coherentElasticModule.form, TNSL_incoherentElasticModule.form, TNSL_incoherentInelasticModule.form ) )

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
