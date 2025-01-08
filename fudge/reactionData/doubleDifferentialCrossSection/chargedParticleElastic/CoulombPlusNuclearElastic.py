# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""     
This module contains a class storing the double differential cross section for charged-particle elastic scattering.
The double differential cross section contain a Rutherford scattering term and maybe nuclear plus interference terms.

This module contains the following classes:
        
    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Form                              | This class represents the double differential cross section for       |
    |                                   | charged-particle elastic scattering.                                  |
    +-----------------------------------+-----------------------------------------------------------------------+
    | CoulombDepositionNotSupported     | This is an exception class that is raised when a method of the        |
    |                                   | :py:class:`Form` cannot perform its calculation.                      |
    +-----------------------------------+-----------------------------------------------------------------------+
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
    """
    This class represents the double differential cross section for charged-particle elastic scattering in one of the following
    three cases:

        -) Only Rutherford (i.e., Coulomb) scattering is represented.
        -) Rutherford and nuclear plus interference are represented.
        -) Rutherford and nuclear amplitude expansion are represented.

    The following table list the primary members of this class:

    +-------------------------------+-------------------------------------------------------------------+
    | Member                        | Description                                                       |
    +===============================+===================================================================+
    | pid                           | The GNDS PoPs id of the product.                                  |
    +-------------------------------+-------------------------------------------------------------------+
    | label                         | The label for this form.                                          |
    +-------------------------------+-------------------------------------------------------------------+
    | productFrame                  | The frame the product data are specified in.                      |
    +-------------------------------+-------------------------------------------------------------------+
    | identicalParticles            | If True, both products in a two-body reaction are identical,      |
    |                               | otherwise they are not.                                           |
    +-------------------------------+-------------------------------------------------------------------+
    | RutherfordScattering          | The Rutherford scattering instance.                               |
    +-------------------------------+-------------------------------------------------------------------+
    | nuclearPlusInterference       | The nuclear + interference instance.                              |
    +-------------------------------+-------------------------------------------------------------------+
    | nuclearAmplitudeExpansion     | The nuclear amplitude expansion instance.                         |
    +-------------------------------+-------------------------------------------------------------------+
    """

    moniker = "CoulombPlusNuclearElastic"
    keyName = 'label'

    subformAttributes = ( 'RutherfordScattering', 'nuclearPlusInterference', 'nuclearAmplitudeExpansion' )

    def __init__(self, pid, label, productFrame=xDataEnumsModule.Frame.centerOfMass, RutherfordScattering=None,
                nuclearPlusInterference=None, nuclearAmplitudeExpansion=None, identicalParticles=False ):
        """
        :param pid:                         The PoPs id of the product.
        :param label:                       The label for the form.
        :param productFrame:                The frame the product data are specified in.
        :param RutherfordScattering:        The Rutherford scattering instance.
        :param nuclearPlusInterference:     The nuclear + interference instance.
        :param nuclearAmplitudeExpansion:   The nuclear amplitude expansion instance.
        :param identicalParticles:          If True, the projectile and target are the same type of particle, otherwise they are not.
        """

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
        """This function returns the spin of the product for identical particles."""

        self.initialize( )
        return( self.__spin )

    @property
    def etaCoefficient( self ) :
        """This function returns the parameter :math:`\\eta \\, \\sqrt{E} = Z_1 \\, Z_2 \\, sqrt{\\alpha^2 \\mu m_1 / 2}`."""

        self.initialize( )
        return( self.__etaCoefficient  )

    @property
    def kCoefficient( self ) :
        """This function returns the coefficient for the particle wave number (i.e., :math:`k(E) / \\sqrt{E}`)."""

        self.initialize( )
        return( self.__kCoefficient )

    @property
    def data( self ) :
        """This function returns the non-Rutherform term."""

        return self.__data

    @property
    def domainMin( self) :
        """Returns the minimum projectile energy for *self*."""

        if self.data is not None: return self.data.domainMin
        return self.RutherfordScattering.domainMin

    @property
    def domainMax( self ) :
        """Returns the maximum projectile energy for *self*."""

        if self.data is not None: return self.data.domainMax
        return self.RutherfordScattering.domainMax

    @property
    def domainUnit( self ) :
        """Returns the energy unit of the projectile."""

        if self.data is not None: return self.data.domainUnit
        return self.RutherfordScattering.domainUnit

    def check( self, info ):
        """This function does a check of *self*."""

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
        """This function does nothing."""

        return 0

    def initialize( self ):
        """
        This function pre-computes some factors used to calculate the Rutherford cross section.
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
        This function returns :math:`d\\sigma / d\\mu` at the specified incident energy if *probability* is False and :math:`P(\mu)` otherwise
        if True..  The :math:`\\mu` domain goes from muMin to *muCutoff*. For identical particles, muMin is set to -*muCutoff* otherwise it is -1.

        :param energy:                          Energy of the projectile.
        :param muCutoff:                        The maximum (and maybe minimum) value of :math:`\\mu` for the returned function.
        :param accuracy:                        The accuracy of the returned function.
        :param epsilon:                         This variable is not used.
        :param excludeRutherfordScattering:     If True, Rutherford scattering is not added to the returned function, otherwise it is.
        :param probability:                     If True :math:`P(\\mu)` is returned otherwise :math:`d\\sigma / d\\mu` is returned.

        :return:                                An instance of :py:class:`angularModule.XYs1d`.
        """

        def dullPoint( mu, epsilon ) :
            """This function is not used and should be deleted."""

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
        This method returns the cross section at (E, mu).

        :param E:                               Incident energy in the lab frame.
        :param mu:                              Scattering angle cosine in the center-of-mass.
        :param phi:                             Scattering azimuthal angle in the center-of-mass.
        :param excludeRutherfordScattering:     If True, only the nuclear and interference terms are included in the returned value.

        :return:                                A python float.
        """

        if( excludeRutherfordScattering ) :
            RS = 0.0
        else :
            RS = self.RutherfordScattering.evaluate( E, mu, phi )
        NS = 0
        if( self.data is not None ) : NS = self.data.evaluate( E, mu, phi )
        return RS + NS

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Ths function is not implemented and executes a raise.

        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this function does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processCoulombPlusNuclearMuCutoff( self, style, energyMin = None, accuracy = 1e-3, epsilon = 1e-6, excludeRutherfordScattering = False ) :
        """
        This function returns the cross section and angular distribution for :math:`\\mu` from muMin to muMax.
        For identical particles, muMin is set to -muMax otherwise it is -1. The value of muMax is the *muCufoff* 
        member of *style*.

        :param style:                           The style for the returned data.
        :param energyMin:                       The minimum projectile energy for calculating the data.
        :param accuracy:                        The accuracy of the returned function.
        :param epsilon:                         This variable is not used.
        :param excludeRutherfordScattering:     If True, only the nuclear and interference terms are included in the returned value.

        :returns:                               An instances :py:class:`angularModule.XYs1d`.
        """

        class Tester :

            def __init__( self, dSigma_dMu, muCutoff, mu_accuracy, relativeTolerance, absoluteTolerance ) :

                self.dSigma_dMu = dSigma_dMu
                self.muCutoff = muCutoff
                self.mu_accuracy = mu_accuracy
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                dSigma_dMu = self.dSigma_dMu( energy, muCutoff, accuracy = self.mu_accuracy, excludeRutherfordScattering = excludeRutherfordScattering )
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
        mu_accuracy = 1e-5
        crossSection = []
        for energy in energies :
            dSigma_dMu = self.dSigma_dMu( energy, muCutoff, accuracy = mu_accuracy, excludeRutherfordScattering = excludeRutherfordScattering )
            crossSection.append([ energy, dSigma_dMu.integrate() ])
        _tester = Tester( self.dSigma_dMu, muCutoff, mu_accuracy, accuracy, accuracy * crossSection[-1][1] )
        crossSection = fudgemathModule.thickenXYList( crossSection, _tester, biSectionMax = 16 )

        crossSectionAxes = crossSectionModule.defaultAxes( self.domainUnit )
        crossSection = crossSectionModule.XYs1d( data = crossSection, axes = crossSectionAxes, label = style.label )

        xys2d = angularModule.XYs2d( axes = angularModule.defaultAxes( self.domainUnit ) )

        crossSectionData = []
        probability = nuclearPlusInterferenceCrossSection is not None
        for energy in energies :
            data = self.dSigma_dMu(energy, muCutoff, accuracy=mu_accuracy, excludeRutherfordScattering=excludeRutherfordScattering, probability=probability)
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
        """
        This functin does nothing.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data. 
        :param indent:          The indentation for any verbosity.
            
        :return:                None.
        """

        print('    processMultiGroup not implemented for distribution form %s.' % self.moniker)
        return( None )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

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
    Custom Exception, returned when calculateAverageProductData().
    """

    pass
