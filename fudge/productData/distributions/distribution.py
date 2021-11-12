# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Distribution class."""

from PoPs import IDs as IDsPoPsModule

from fudge import abstractClasses as abstractClassesModule
from xData import standards as standardsModule

from . import angular as angularModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import energyAngularMC as energyAngularMCModule
from . import angularEnergyMC as angularEnergyMCModule
from . import KalbachMann as KalbachMannModule
from . import angularEnergy as angularEnergyModule
from . import LLNL_angularEnergy as LLNL_angularEnergyModule
from . import uncorrelated as uncorrelatedModule
from . import Legendre as LegendreModule
from . import photonScattering as photonScatteringModule
from . import reference as referenceModule
from . import multiGroup as multiGroupModule
from . import unspecified as unspecifiedModule
from . import branching3d as branching3dModule

# probably missing stuff from photonScattering.py.

__metaclass__ = type

class component( abstractClassesModule.component ) :

    moniker = 'distribution'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( angularModule.form, angularModule.twoBodyForm,
                KalbachMannModule.form,
                energyAngularModule.form, energyAngularMCModule.form,
                angularEnergyModule.form, angularEnergyMCModule.form,
                LLNL_angularEnergyModule.LLNLAngularEnergyForm,
                uncorrelatedModule.form, LegendreModule.form, referenceModule.form,
                referenceModule.CoulombPlusNuclearElastic, referenceModule.thermalNeutronScatteringLaw,
                photonScatteringModule.coherentPhotonScattering.form, photonScatteringModule.incoherentPhotonScattering.form, 
                multiGroupModule.form, unspecifiedModule.form, branching3dModule.form ) )

    def energySpectrumAtEnergy( self, energyIn, frame, **kwargs ) :
        """Returns the energy spectrum in the lab frame for the specified incident energy."""

        styleLabel = kwargs.get( 'styleLabel', self.evaluated.label )
        form = self[styleLabel]
        if( hasattr( form, 'energySpectrumAtEnergy' ) ) :
            if( frame == standardsModule.frames.centerOfMassToken ) :
                if( form.productFrame == standardsModule.frames.labToken ) : form = None
        else :
            form = None

        if( form is not None ) :
            return( form.energySpectrumAtEnergy( energyIn, frame ) )
        else :
            form = self[styleLabel]
            if( hasattr( form, 'energySpectrumAtEnergy' ) ) :
                print( '        WARNING: lab to center-of-mass translation not supported.' )
            else :
                print( '        WARNING: distribution "%s" does not have energySpectrumAtEnergy method.' % form.moniker )
            print( '            %s' % self.toXLink( ) )
            return( energyModule.XYs1d( axes = energyModule.defaultAxes( form.domainUnit ) ) )

    def getSpectrumAtEnergy( self, energy ) :
        """This method is deprecated, use energySpectrumAtEnergy instead. Returns the energy spectrum for self at projectile energy."""

        return( self.energySpectrumAtEnergy( energy, standardsModule.frames.labToken ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        form = style.findFormMatchingDerivedStyle( self )
        if( form is None ) : raise Exception( 'No matching style' )
        return( form.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def check( self, info ):
        """check all distribution forms"""
        from fudge import warning
        warnings = []

        for form in self:

            if info['isTwoBody']:
                if( form.productFrame != standardsModule.frames.centerOfMassToken ) :
                    warnings.append( warning.wrong2BodyFrame( form ) )

                if form.moniker not in (angularModule.twoBodyForm.moniker,
                                        referenceModule.form.moniker,
                                        referenceModule.CoulombPlusNuclearElastic.moniker,
                                        unspecifiedModule.form.moniker):
                    warnings.append( warning.wrongDistributionComponent( form.moniker, '2-body' ) )
            else:
                if form.moniker in (angularModule.twoBodyForm.moniker,
                                    angularModule.form.moniker,
                                    energyModule.form.moniker):
                    warnings.append( warning.wrongDistributionComponent( form.moniker, 'N-body' ) )

            def checkSubform( subform, contextMessage ):
                distributionErrors = []
                if hasattr(subform, 'domainMin') and (subform.domainMin, subform.domainMax) != info['crossSectionDomain']:
                    domain = (subform.domainMin, subform.domainMax)
                    # For gamma products, domainMin should be >= cross section start, upper bounds should match.
                    if( self.ancestor.id == IDsPoPsModule.photon ) :
                        startRatio = subform.domainMin / info['crossSectionDomain'][0]
                        endRatio = subform.domainMax / info['crossSectionDomain'][1]
                        if (startRatio < 1-standardsModule.floats.epsilon or endRatio < 1-standardsModule.floats.epsilon
                                or endRatio > 1+standardsModule.floats.epsilon):
                            distributionErrors.append( warning.domain_mismatch(
                                    *(domain + info['crossSectionDomain']), obj=subform ) )
                    # For all other products, check lower and upper edges: only warn if they disagree by > eps
                    else:
                        for e1,e2 in zip(domain, info['crossSectionDomain']):
                            ratio = e1 / e2
                            if (ratio < 1-standardsModule.floats.epsilon or ratio > 1+standardsModule.floats.epsilon):
                                distributionErrors.append( warning.domain_mismatch(
                                        *(domain + info['crossSectionDomain']), obj=subform ) )
                                break
                    if not hasattr(subform,'check'):
                        distributionErrors.append( warning.NotImplemented(subform.moniker, subform ) )
                        if info['failOnException']:
                            raise NotImplementedError("Checking distribution form '%s'" % subform.moniker)
                else:
                    distributionErrors += subform.check( info )
                if distributionErrors:
                    warnings.append( warning.context( contextMessage + " - %s:" % subform.moniker, distributionErrors) )

            if isinstance(form, uncorrelatedModule.form):
                for subformName in ('angularSubform','energySubform'):
                    subform = getattr(form, subformName ).data
                    checkSubform( subform, 'uncorrelated - ' + subformName.replace('Subform','') )

            elif isinstance(form, KalbachMannModule.form):
                checkSubform( form, form.moniker )

            else:
                for subform in form.subforms:
                    checkSubform( subform, form.moniker )

        return warnings

    def diff( self, other, diffResults ) :

        if( self.hasData( ) != other.hasData( ) ) :
            if( self.hasData( ) ) :
                diffResults.append( 'Distribution unspecified - 2', '', self.toXLink( ), other.toXLink( ) )
            else :
                diffResults.append( 'Distribution unspecified - 1', '', self.toXLink( ), other.toXLink( ) )

    def patch( self, other ) :

        pass

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. Need ability to find specific distribution component
        """
        if attribute is not None:
            for entity in self:
                if entity.moniker == entityName and getattr(entity,attribute) == value:
                    return entity
        else:
            for entity in self:
                if entity.moniker == entityName:
                    return entity
        return abstractClassesModule.component.findEntity( self, entityName, attribute, value )

    def hasData( self ) :
        """
        Returns False if self's only has unspecified form; otherwise, returns True.
        """

        for form in self :
            if( not( isinstance( form, unspecifiedModule.form ) ) ) : return( True )
        return( False )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = standardsModule.frames.productToken, LegendreOrder = 0 ) :

        if( len( self ) > 0 ) :
            form = self[0]
#            if( form.productFrame == standardsModule.frames.centerOfMassToken ) : return( 0.0 )
            if( hasattr( form, 'integrate' ) ) :
                return( form.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )
            else :
                print( 'missing integrate', type( form ) )

        return( 0.0 )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.evaluated.toPointwise_withLinearXYs( **kwargs ) )
