# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from pqu import PQU as PQUModule
from xData import multiD_XYs as multiD_XYsModule

from LUPY import ancestry as ancestryModule

from fudge.core.math import fudgemath as fudgemathModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule

class RutherfordScattering( ancestryModule.AncestryIO ):
    """
    Stores the pure-Coulomb contribution to charged-particle elastic scattering
    """

    moniker = "RutherfordScattering"

    def __init__(self, domainMin=None, domainMax=None, domainUnit=None ):

        ancestryModule.AncestryIO.__init__(self)

        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

    @property
    def etaCoefficient( self ) :

        return( self.ancestor.etaCoefficient )

    @property
    def spin( self ) :

        return( self.ancestor.spin )

    @property
    def kCoefficient( self ) :

        return( self.ancestor.kCoefficient )

    @property
    def identicalParticles( self ) :

        return( self.ancestor.identicalParticles )

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[ self.domainUnit ]
            factor = PQUModule.PQU( 1., self.domainUnit ).getValueAs( newUnit )
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit

    def crossSectionVersusEnergy(self, muMax, accuracy=1e-3, energyMin=None, energyMax=None):
        """
        Returns the Rutherford cross section.
        """

        class Tester:

            def __init__(self, dSigma_dMu, muMax, relativeTolerance, absoluteTolerance):

                self.dSigma_dMu = dSigma_dMu
                self.muMax = muMax
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX(self, energy):

                dSigma_dMu = self.dSigma_dMu(energy, muMax=muMax, accuracy=self.relativeTolerance)
                return dSigma_dMu.integrate()

        energyUnit, energyMin, energyMax = self.energyDomain(energyMin, energyMax)
        energies = [energyMin, math.sqrt(energyMin * energyMax), energyMax]

        crossSection = []
        for energy in energies:
            dSigma_dMu = self.dSigma_dMu(energy, muMax=muMax, accuracy=accuracy)
            crossSection.append([energy, dSigma_dMu.integrate()])

        tester = Tester(self.dSigma_dMu, muMax, accuracy, accuracy * crossSection[-1][1])
        crossSection = fudgemathModule.thickenXYList(crossSection, tester, biSectionMax=16)

        return crossSectionModule.XYs1d(data=crossSection, axes=crossSectionModule.defaultAxes(energyUnit))


    def dSigma_dMuVersusEnergy(self, muMax, accuracy=1e-3, energyMin=None, energyMax=None):
        """
        Calls **dSigma_dMu** on self at various energies.
        """

        energyUnit, energyMin, energyMax = self.energyDomain(energyMin, energyMax)

        energyFactor = math.sqrt(10.)
        energies = [energyMin]
        energy = energyMin
        while True:
            energy *= energyFactor
            if energy > 0.9 * energyMax: break
            energies.append(energy)
        energies.append(energyMax)

        xys2d = multiD_XYsModule.XYs2d(axes=self.defaultAxes(energyUnit))
        for energy in energies:
            xys2d.append(self.dSigma_dMu(energy, accuracy=accuracy, muMax=muMax))

        return xys2d

    def dSigma_dMu(self, energy, accuracy=1e-3, muMax=0.999, probability=False):
        """
        Returns d(Sigma)/d(mu) at the specified incident energy.

        :param energy:      Energy of the projectile.
        :param accuracy:    The accuracy of the returned *dSigma_dMu*.
        :param muMax:       Slices the upper domain mu to this value.
        :param probability: If **True** P(mu) is returned instead of d(Sigma)/d(mu).

        :return:            d(Sigma)/d(mu) at *energy*.
        """

        class Tester :

            def __init__( self, evaluate, energy, relativeTolerance, absoluteTolerance ) :

                self.evaluate = evaluate
                self.energy = energy
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, mu ) :

                return( self.evaluate( self.energy, mu ) )

        muMin = -1.0
        if self.identicalParticles: muMin = 1e-8
        _dSigma_dMu = [ [ muMin, self.evaluate( energy, muMin ) ], [ muMax, self.evaluate( energy, muMax ) ] ]
        average = sum( y for x, y in _dSigma_dMu ) / len( _dSigma_dMu )
        _tester = Tester( self.evaluate, energy, accuracy, accuracy * average )
        _dSigma_dMu = fudgemathModule.thickenXYList(_dSigma_dMu, _tester, biSectionMax=16)
        if self.identicalParticles: _dSigma_dMu = [[-1.0, 0.0], [0.0, 0.0]] + _dSigma_dMu
        if _dSigma_dMu[-1][0] < 1:
            _dSigma_dMu += [[_dSigma_dMu[-1][0] * (1 + 1e-5 ), 0.0], [1.0, 0.0]]
        _dSigma_dMu = 2.0 * math.pi * angularModule.XYs1d(_dSigma_dMu, axes=self.defaultAxes(self.domainUnit))
        _dSigma_dMu.axes = self.defaultAxes( self.domainUnit )

        if probability: _dSigma_dMu = _dSigma_dMu.normalize()

        _dSigma_dMu.outerDomainValue = energy
        return( _dSigma_dMu )

    def energyDomain(self, energyMin=None, energyMax=None):
        """
        Returns the tuple (energyUnit, energyMin, energyMax).
        """

        energyUnit = self.domainUnit
        if energyUnit is None: energyUnit = self.ancestor.domainUnit

        if energyMin is None: energyMin = self.domainMin
        if energyMin is None: energyMin = PQUModule.PQU('1e-4 MeV').getValueAs(energyUnit)

        if energyMax is None: energyMax = self.domainMax
        if energyMax is None: energyMax = PQUModule.PQU('30 MeV').getValueAs(energyUnit)

        return energyUnit, energyMin, energyMax

    def evaluate( self, energy, mu, phi = 0.0 ) :
        """
        Returns the dSigma/dOmega(energy,mu,phi) for Rutherford scattering. Note, Rutherford scattering is independent of phi
        but phi is listed as an argument to be consistent with others dSigma/dOmega that may depend on phi.

        :param energy: Incident energy.
        :param mu: Scattering angle cosine.
        :param phi: The scattering azimuthal angle.

        :return: Differential cross section dSigma/dOmega at point energy,mu,phi (in b/sr).
        """

        eta = self.etaCoefficient / math.sqrt( energy )
        k = self.kCoefficient * math.sqrt( energy )
        if( self.identicalParticles ) :
            twiceSpin = 2.0 * float( self.spin )
            term1 = ( 2 * eta**2 ) / ( k**2 * ( 1.0 - mu**2 ) )
            term2 = ( 1 + mu**2 ) / ( 1 - mu**2 ) + (-1)**twiceSpin / ( twiceSpin + 1 ) * math.cos( eta * math.log( ( 1 + mu ) / ( 1 - mu ) ) )
            return term1 * term2
        else:
            return eta**2 / (k**2 * (1-mu)**2 )

    def toXML_strList( self, indent='', **kwargs ):

        attrStr = ''
        for attribute in ('domainMin','domainMax','domainUnit'):
            if getattr(self,attribute) is not None: attrStr += ' %s="%s"' % (attribute, getattr(self,attribute))
        return ['%s<%s%s/>' % (indent, self.moniker, attrStr)]

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        domainMin, domainMax = element.get('domainMin'), element.get('domainMax')
        if domainMin is not None:
            domainMin, domainMax = map(float, (domainMin, domainMax) )
        return cls( domainMin, domainMax, element.get('domainUnit') )

    @staticmethod
    def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
