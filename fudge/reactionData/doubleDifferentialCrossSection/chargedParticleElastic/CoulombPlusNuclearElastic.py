# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
This module contains a form for storing the differential cross section for charged-particle elastic scattering.
Internally, data can be represented three ways:
 - pure Rutherford scattering
 - Rutherford scattering along with Legendre expansions for nuclear scattering
        and for real and imaginary nuclear/Coulomb interference
 - Rutherford scattering along with effective cross sections and distributions,
        obtained by summing the nuclear and interference terms
"""

import math

from pqu import PQU as PQUModule
from xData import enums as xDataEnumsModule

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideFamilyModule

from fudge.core.math import fudgemath as fudgemathModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData.distributions import angular as angularModule

from .. import base as baseModule

from . import misc as miscModule
from . import RutherfordScattering as RutherfordScatteringModule
from . import nuclearAmplitudeExpansion as nuclearAmplitudeExpansionModule
from . import nuclearPlusInterference as nuclearPlusInterferenceModule


class Form( baseModule.Form ):

    moniker = "CoulombPlusNuclearElastic"
    keyName = 'label'

    subformAttributes = ( 'RutherfordScattering', 'nuclearPlusInterference', 'nuclearAmplitudeExpansion' )

    def __init__(self, pid, label, productFrame=xDataEnumsModule.Frame.centerOfMass, RutherfordScattering=None,
                nuclearPlusInterference=None, nuclearAmplitudeExpansion=None, identicalParticles=False ):

        if( RutherfordScattering is None) : RutherfordScattering = RutherfordScatteringModule.RutherfordScattering( )
        if( not isinstance( RutherfordScattering, ( RutherfordScatteringModule.RutherfordScattering, type( None ) ) ) ) :
            raise TypeError( "Invalid nuclearPlusInterference type: '%s' not allowed in %s" % ( type( RutherfordScattering ), self.moniker ) )

        if( not isinstance( nuclearPlusInterference, ( nuclearPlusInterferenceModule.NuclearPlusInterference, type( None ) ) ) ) :
            raise TypeError( "Invalid NuclearPlusInterference type: '%s' not allowed in %s" % ( type( nuclearPlusInterference ), self.moniker ) )

        if( not isinstance( nuclearAmplitudeExpansion, ( nuclearAmplitudeExpansionModule.NuclearAmplitudeExpansion, type( None ) ) ) ) :
            raise TypeError( "Invalid nuclearAmplitudeExpansion type: '%s' not allowed in %s" % ( type( nuclearAmplitudeExpansion ), self.moniker ) )

        baseModule.Form.__init__( self, pid, label, productFrame, ( RutherfordScattering, nuclearPlusInterference, nuclearAmplitudeExpansion ), 
                identicalParticles = identicalParticles )

        self.__data = self.nuclearPlusInterference
        if( self.__data is None ) : self.__data = self.nuclearAmplitudeExpansion

        self.__etaCoefficient = None
        self.__spin = None                      # Only defined if identicalParticles is True.

    @property
    def spin( self ) :

        self.initialize( )
        return( self.__spin )

    @property
    def etaCoefficient( self ) :

        self.initialize( )
        return( self.__etaCoefficient  )

    @property
    def kCoefficient( self ) :

        self.initialize( )
        return( self.__kCoefficient )

    @property
    def data( self ) :

        return self.__data

    @property
    def domainMin( self) :

        if self.data is not None: return self.data.domainMin
        return self.RutherfordScattering.domainMin

    @property
    def domainMax( self ) :

        if self.data is not None: return self.data.domainMax
        return self.RutherfordScattering.domainMax

    @property
    def domainUnit( self ) :

        if self.data is not None: return self.data.domainUnit
        return self.RutherfordScattering.domainUnit

    def check( self, info ):

        from fudge import warning
        warnings = []

        RS = info['reactionSuite']
        target = RS.PoPs[ RS.target ]
        projectile = RS.PoPs[ RS.projectile ]
        identicalParticles = target is projectile
        if identicalParticles and not self.identicalParticles:
            warnings.append( warning.MissingCoulombIdenticalParticlesFlag() )
        elif not identicalParticles and self.identicalParticles:
            warnings.append( warning.IncorrectCoulombIdenticalParticlesFlag(
                RS.projectile, RS.target ) )

        if self.data is not None:
            dataWarnings = self.data.check( info )
            if dataWarnings:
                warnings.append( warning.Context('%s:' % self.data.moniker, dataWarnings) )

        return warnings

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """This method does nothing."""

        return 0

    def initialize( self ):
        """
        Pre-compute some factors used to calculate the Rutherford cross section.
        """

        if( self.__etaCoefficient is not None ) : return        # Already initialized.

        reactionSuite = self.rootAncestor

        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        if( isinstance( projectile, nuclideFamilyModule.Particle ) ) : projectile = projectile.nucleus

        targetID = reactionSuite.target
        if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]

        Z1 = chemicalElementMiscPoPsModule.ZAInfo( projectile )[0]
        Z2 = chemicalElementMiscPoPsModule.ZAInfo( target )[0]
        if( self.data is None ) :
            domainUnit = reactionSuite.reactions[0].domainUnit
        else :
            domainUnit = self.data.domainUnit
        mass1 = projectile.getMass( '%s / c**2' % domainUnit )
        mass2 = target.getMass( '%s / c**2' % domainUnit )
        if( self.identicalParticles ) : self.__spin = projectile.spin[0].value

        hbar_c = PQUModule.PQU( 1, 'hbar * c' ).getValueAs( 'm * %s' % domainUnit )
        alpha = PQUModule.PQU( '1', 'e * e / ( hbar*c * 4 * pi * eps0 )' ).getValueAs( '' )

        self.__etaCoefficient = Z1 * Z2 * alpha * math.sqrt( mass1 / 2 )
        A = mass2 / mass1
        self.__kCoefficient = (A / (A + 1)) * math.sqrt( 2 * mass1 ) / hbar_c * 1e-14        # 1e-14 = sqrt( barn )

    def dSigma_dMu(self, energy, muCutoff, accuracy=1e-3, epsilon=1e-6, excludeRutherfordScattering=False, probability=False):
        """
        Returns d(Sigma)/d(mu) at the specified incident energy if probability is **False** and P(mu) otherwise.

        :param energy:      Energy of the projectile.
        :param accuracy:    The accuracy of the returned *dSigma_dMu*.
        :param muMax:       Slices the upper domain mu to this value.
        :param probability: If **True** P(mu) is returned instead of d(Sigma)/d(mu).

        :return:            d(Sigma)/d(mu) at *energy*.
        """

        def dullPoint( mu, epsilon ) :

            if( mu < 0.0 ) : epsilon *= -1
            return( mu * ( 1 + epsilon ) )

        epsilon = max( 1e-15, min( 0.1, abs( epsilon ) ) )

        if( abs( muCutoff ) >= 1.0 ) : raise ValueError( 'muCutoff = %.17e must be in the range ( -1, 1 ).' % muCutoff )
        muMin = -1.0
        if( self.identicalParticles ) :
            muCutoff = abs( muCutoff )
            muMin = 0.0

        if( ( self.data is None ) or ( energy < self.data.domainMin ) ) :
            _dSigma_dMu = angularModule.XYs1d( [ [ -1.0, 0.0 ], [ 1.0, 0.0 ] ], axes = miscModule.defaultAxes( self.domainUnit ) )
        else :
            _dSigma_dMu = self.data.dSigma_dMu(energy, accuracy=accuracy, muMax=muCutoff, probability=probability)
        if not excludeRutherfordScattering:
            _dSigma_dMu += self.RutherfordScattering.dSigma_dMu(energy, accuracy=accuracy, muMax=muCutoff)
        _dSigma_dMu = _dSigma_dMu.thin( accuracy = 1e-3 )

        return( _dSigma_dMu )

    def evaluate( self, E, mu, phi = 0.0, excludeRutherfordScattering = False ) :
        """
        Compute the cross section at (E, mu), including Coulomb, nuclear and interference terms.

        :param E: incident energy in the lab frame.
        :param mu: scattering angle cosine in the center-of-mass.
        :param phi: scattering azimuthal angle in the center-of-mass.
        :param excludeRutherfordScattering: If True, only the nuclear and interference terms are included in the returned value.
        :return:
        """

        if( excludeRutherfordScattering ) :
            RS = 0.0
        else :
            RS = self.RutherfordScattering.evaluate( E, mu, phi )
        NS = 0
        if( self.data is not None ) : NS = self.data.evaluate( E, mu, phi )
        return RS + NS

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processCoulombPlusNuclearMuCutoff( self, style, energyMin = None, accuracy = 1e-3, epsilon = 1e-6, excludeRutherfordScattering = False ) :

        class Tester :

            def __init__( self, dSigma_dMu, muCutoff, relativeTolerance, absoluteTolerance ) :

                self.dSigma_dMu = dSigma_dMu
                self.muCutoff = muCutoff
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                dSigma_dMu = self.dSigma_dMu( energy, muCutoff, accuracy = self.relativeTolerance, excludeRutherfordScattering = excludeRutherfordScattering )
                return dSigma_dMu.integrate()

        nuclearPlusInterferenceCrossSection = None
        if( self.nuclearPlusInterference is None ) :
            if( self.nuclearAmplitudeExpansion is None ) :
                if( excludeRutherfordScattering ) : return( None, None )
                energies = [ self.RutherfordScattering.domainMin, self.RutherfordScattering.domainMax ]
                energies.insert( 1, math.sqrt( energies[0] * energies[1] ) )
            else :
                nuclearTerm = self.nuclearAmplitudeExpansion.nuclearTerm.data
                if( isinstance( nuclearTerm, angularModule.XYs2d ) ) :
                    energies = nuclearTerm.domainGrid
                elif( isinstance( nuclearTerm, angularModule.Regions2d ) ) :
                    energies = []
                    for region in nuclearTerm : energies += region.domainGrid
                    energies = sorted( set( energies ) )
                else :
                    raise Exception( 'distribution type "%s" not supported' % type( nuclearTerm ) )
        else :
            energies = self.nuclearPlusInterference.distribution.data.domainGrid
            nuclearPlusInterferenceCrossSection = self.nuclearPlusInterference.crossSection.data.toPointwise_withLinearXYs(lowerEps=1e-6)

        if not excludeRutherfordScattering:
            RutherfordEnergies = energies.copy( )
            RutherfordEnergyMin = self.RutherfordScattering.domainMin
            if( RutherfordEnergyMin is None ) : RutherfordEnergyMin = PQUModule.PQU( 1e-4, 'MeV' ).getValueAs( self.domainUnit )
            if( energyMin is not None ) :
                if( energyMin < RutherfordEnergyMin ) : RutherfordEnergyMin = energyMin

            nonRutherfordEnergyMin = energies[0]
            index = 0
            while( RutherfordEnergyMin < 0.8 * nonRutherfordEnergyMin ) :
                RutherfordEnergies.insert( index, RutherfordEnergyMin )
                index += 1
                RutherfordEnergyMin *= 1.4142135623731

            energies = RutherfordEnergies

        muCutoff = style.muCutoff
        crossSection = []
        for energy in energies :
            dSigma_dMu = self.dSigma_dMu( energy, muCutoff, accuracy = accuracy, excludeRutherfordScattering = excludeRutherfordScattering )
            crossSection.append([ energy, dSigma_dMu.integrate() ])
        _tester = Tester( self.dSigma_dMu, muCutoff, accuracy, accuracy * crossSection[-1][1] )
        crossSection = fudgemathModule.thickenXYList( crossSection, _tester, biSectionMax = 16 )

        crossSectionAxes = crossSectionModule.defaultAxes( self.domainUnit )
        crossSection = crossSectionModule.XYs1d( data = crossSection, axes = crossSectionAxes, label = style.label )

        xys2d = angularModule.XYs2d( axes = angularModule.defaultAxes( self.domainUnit ) )

        crossSectionData = []
        probability = nuclearPlusInterferenceCrossSection is not None
        for energy in energies :
            data = self.dSigma_dMu(energy, muCutoff, accuracy=accuracy, excludeRutherfordScattering=excludeRutherfordScattering, probability=probability)
            if( excludeRutherfordScattering ) : data = data.clip( rangeMin = 0.0 )

            xSec = data.integrate()
            crossSectionData.append( [ energy, xSec ] )

            if xSec == 0.0:
                data = [[-1.0, 0.5], [1.0, 0.5]]
                if self.identicalParticles:
                    data = [[-1.0, 0.0], [-1e-12, 0.0], [1e-12, 1.0], [1.0, 1.0]]
            else :
                data /= xSec

            xys1d = angularModule.XYs1d(data=data, axes=xys2d.axes, outerDomainValue=energy)
            xys2d.append(xys1d)

        if excludeRutherfordScattering:
            if nuclearPlusInterferenceCrossSection is not None: crossSectionData = nuclearPlusInterferenceCrossSection
            crossSection = crossSectionModule.XYs1d(data=crossSectionData, axes=crossSectionAxes, label=style.label)

        return( crossSection, xys2d )

    def processMultiGroup( self, style, tempInfo, indent ) :

        print('    processMultiGroup not implemented for distribution form %s.' % self.moniker)
        return( None )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        RutherfordScattering = None
        nuclearPlusInterference = None
        nuclearAmplitudeExpansion = None
        for child in element:
            if child.tag == RutherfordScatteringModule.RutherfordScattering.moniker:
                RutherfordScattering = RutherfordScatteringModule.RutherfordScattering.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == nuclearPlusInterferenceModule.NuclearPlusInterference.moniker:
                nuclearPlusInterference = nuclearPlusInterferenceModule.NuclearPlusInterference.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == nuclearAmplitudeExpansionModule.NuclearAmplitudeExpansion.moniker:
                nuclearAmplitudeExpansion = nuclearAmplitudeExpansionModule.NuclearAmplitudeExpansion.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise TypeError( "Encountered unexpected element '%s' in %s" % ( child.tag, element.tag ) )
        subForms = ( RutherfordScattering, nuclearPlusInterference, nuclearAmplitudeExpansion )
        identicalParticles = element.get( 'identicalParticles', '' ) == 'true'

        Coul = cls( element.get( 'pid' ), element.get( 'label' ), element.get( 'productFrame' ), identicalParticles = identicalParticles,
                RutherfordScattering = RutherfordScattering, nuclearPlusInterference = nuclearPlusInterference, nuclearAmplitudeExpansion = nuclearAmplitudeExpansion )

        xPath.pop( )

        return Coul

class CoulombDepositionNotSupported( Exception ):
    """
    Custom Exception, returned when calculateAverageProductData() is called for
    nuclearAmplitudeExpansion or nuclearPlusInterference
    """
    pass
