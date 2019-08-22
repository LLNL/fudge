# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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
import numpy

from pqu import PQU as PQUModule

from xData import standards as standardsModule
from xData import axes as axesModule

from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideFamilyModule

from fudge.gnds import abstractClasses as abstractClassesModule
from fudge.core.math import fudgemath as fudgemathModule
from fudge.gnds.reactionData import crossSection as crossSectionModule
from fudge.gnds.productData.distributions import angular as angularModule

from .RutherfordScattering import RutherfordScattering
from .nuclearAmplitudeExpansion import nuclearAmplitudeExpansion
from .nuclearPlusInterference import nuclearPlusInterference

__metaclass__ = type

class form( abstractClassesModule.form ):

    moniker = "CoulombPlusNuclearElastic"

    def __init__( self, label, pid, data, productFrame = standardsModule.frames.centerOfMassToken,
                  RutherfordTerm = None, identicalParticles = False ):

        abstractClassesModule.form.__init__( self )
        self.label = label
        self.__pid = pid
        self.__productFrame = productFrame
        self.__identicalParticles = identicalParticles

        if RutherfordTerm is not None:
            self.RutherfordScattering = RutherfordTerm
        else:
            self.RutherfordScattering = RutherfordScattering()
        self.RutherfordScattering.setAncestor( self )
        self.data = data

        self.__etaCoefficient = None
        self.__spin = None                      # Only defined if identicalParticles is True.

    @property
    def pid(self):

        return self.__pid

    @property
    def productFrame(self):

        return self.__productFrame

    @property
    def identicalParticles( self ) :

        return self.__identicalParticles

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
    def data(self):

        return self.__data

    @data.setter
    def data(self, value):
        if not isinstance(value, (nuclearPlusInterference,
                                  nuclearAmplitudeExpansion, type(None))):
            raise TypeError("Data type '%s' not allowed in %s" % (type(value),self.moniker))
        if value is not None: value.setAncestor(self)
        self.__data = value

    # two additional ways to access data (meant be user-friendly)
    @property
    def nuclearPlusInterference(self):
        if isinstance( self.data, nuclearPlusInterference ): return self.data
        return None

    @property
    def nuclearAmplitudeExpansion(self):
        if isinstance( self.data, nuclearAmplitudeExpansion ): return self.data
        return None

    @property
    def domainMin(self):
        if self.data is not None: return self.data.domainMin
        return self.RutherfordScattering.domainMin

    @property
    def domainMax(self):
        if self.data is not None: return self.data.domainMax
        return self.RutherfordScattering.domainMax

    @property
    def domainUnit(self):
        if self.data is not None: return self.data.domainUnit
        return self.RutherfordScattering.domainUnit

    def check( self, info ):

        from fudge.gnds import warning
        warnings = []

        RS = info['reactionSuite']
        target = RS.PoPs[ RS.target ]
        projectile = RS.PoPs[ RS.projectile ]
        identicalParticles = target is projectile
        if identicalParticles and not self.identicalParticles:
            warnings.append( warning.missingCoulombIdenticalParticlesFlag() )
        elif not identicalParticles and self.identicalParticles:
            warnings.append( warning.incorrectCoulombIdenticalParticlesFlag(
                RS.projectile, RS.target ) )

        if self.data is not None:
            dataWarnings = self.data.check( info )
            if dataWarnings:
                warnings.append( warning.context('%s:' % self.data.moniker, dataWarnings) )

        return warnings

    def convertUnits( self, unitMap ):

        self.RutherfordScattering.convertUnits( unitMap )
        if self.data is not None:
            self.data.convertUnits( unitMap )

    def initialize( self ):
        """
        Pre-compute some factors used to calculate the Rutherford cross section.
        """

        if( self.__etaCoefficient is not None ) : return        # Already initialized.

        reactionSuite = self.getRootAncestor()

        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        if( isinstance( projectile, nuclideFamilyModule.particle ) ) : projectile = projectile.nucleus

        targetID = reactionSuite.target
        if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]

        Z1 = chemicalElementMiscPoPsModule.ZAInfo(projectile)[0]
        Z2 = chemicalElementMiscPoPsModule.ZAInfo(target)[0]
        mass1 = projectile.getMass( '%s / c**2' % self.data.domainUnit )
        mass2 = target.getMass( '%s / c**2' % self.data.domainUnit )
        if( self.identicalParticles ) : self.__spin = projectile.spin[0].value

        hbar_c = PQUModule.PQU( 1, 'hbar * c' ).getValueAs( 'm * %s' % self.data.domainUnit )
        alpha = PQUModule.PQU( '1', 'e * e / ( hbar*c * 4 * pi * eps0 )' ).getValueAs( '' )

        self.__etaCoefficient = Z1 * Z2 * alpha * numpy.sqrt( mass1 / 2 )
        A = mass2 / mass1
        self.__kCoefficient = (A / (A + 1)) * numpy.sqrt( 2 * mass1 ) / hbar_c * 1e-14        # 1e-14 = sqrt( barn )

    def dSigma_dmu( self, energy, muCutoff, accuracy = 1e-3, epsilon = 1e-6 ) :

        class tester :

            def __init__( self, evaluate, energy, relativeTolerance, absoluteTolerance ) :

                self.evaluate = evaluate
                self.energy = energy
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, mu ) :

                return( 2 * math.pi * self.evaluate( self.energy, mu ) )

        def dullPoint( mu, epsilon ) :

            if( mu < 0.0 ) : epsilon *= -1
            return( mu * ( 1 + epsilon ) )


        epsilon = max( 1e-15, min( 0.1, abs( epsilon ) ) )

        if( abs( muCutoff ) >= 1.0 ) : raise ValueError( 'muCutoff = %.17e must be in the range ( -1, 1 ).' % muCutoff )
        muMin = -1.0
        if( self.identicalParticles ) :
            if( muCutoff < 0.0 ) : muCutoff = -muCutoff
            muMin = -muCutoff

        _dSigma_dmu = [ [ muMin, 2 * math.pi * self.evaluate( energy, muMin ) ], [ muCutoff, 2 * math.pi * self.evaluate( energy, muCutoff ) ] ]
        average = sum( [ y for x, y in _dSigma_dmu ] ) / len( _dSigma_dmu )

        _tester = tester( self.evaluate, energy, accuracy, accuracy * average )
        _dSigma_dmu = angularModule.XYs1d( fudgemathModule.thickenXYList( _dSigma_dmu, _tester, biSectionMax = 10 ) )

        xys1d = _dSigma_dmu.copyDataToXYs( )

        muStart = dullPoint( xys1d[0][0], -epsilon )
        if( muStart > -1.0 ) : xys1d.insert( 0, [ muStart, 0.0 ] )
        if( xys1d[0][0] > -1.0 ) : xys1d.insert( 0, [ -1.0, 0.0 ] )

        muEnd = dullPoint( xys1d[-1][0], epsilon )
        if( muEnd < 1.0 ) : xys1d.append( [ muEnd, 0.0 ] )
        if( xys1d[-1][0] < 1.0 ) : xys1d.append( [ 1.0, 0.0 ] )

        domainUnit = self.domainUnit
        if( self.nuclearPlusInterference is None ) :
            rangeUnit = 'b'
        else :
            rangeUnit = self.nuclearPlusInterference.crossSection.data.rangeUnit
        axes = axesModule.axes( rank = 3 )
        axes[0] = axesModule.axis( 'dSigma_dmu', 0, rangeUnit )
        axes[1] = axesModule.axis( 'mu', 1, '' )
        axes[2] = axesModule.axis( 'energy_in', 2, domainUnit )

        _dSigma_dmu = angularModule.XYs1d( xys1d, value = energy, axes = axes )

        return( _dSigma_dmu )

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
        Compute the cross section at (E, mu), including Coulomb, nuclear and interference terms.

        :param E: incident energy (lab frame).
        :param mu: scattering angle cosine (in center-of-mass)
        :return:
        """

        RS = self.RutherfordScattering.evaluate( E, mu, phi )
        NS = 0
        if( self.data is not None ) : NS = self.data.evaluate( E, mu, phi )
        return RS + NS

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processCoulombPlusNuclearMuCutoff( self, style, energyMin = None, accuracy = 1e-3, epsilon = 1e-6 ) :

        class tester :

            def __init__( self, dSigma_dmu, muCutoff, relativeTolerance, absoluteTolerance ) :

                self.dSigma_dmu = dSigma_dmu
                self.muCutoff = muCutoff
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                dSigma_dmu = self.dSigma_dmu( energy, muCutoff, accuracy = self.relativeTolerance )
                return( float( dSigma_dmu.integrate( ) ) )

        muCutoff = style.muCutoff

        if( self.nuclearPlusInterference is None ) :
            nuclearTerm = self.nuclearAmplitudeExpansion.nuclearTerm.data
            if( isinstance( nuclearTerm, angularModule.XYs2d ) ) :
                energies = nuclearTerm.domainGrid
            elif( isinstance( nuclearTerm, angularModule.regions2d ) ) :
                energies = []
                for region in nuclearTerm : energies += region.domainGrid
                energies = sorted( set( energies ) )
            else :
                raise Exception( 'distribution type "%s" not supported' % type( nuclearTerm ) )
        else :
            energies = self.nuclearPlusInterference.crossSection.data.domainGrid

        if( energyMin is not None ) :
            if( energyMin < energies[0] ) : energies.insert( energyMin )

        crossSection = []
        for energy in energies :
            crossSection.append( [ energy, float( self.dSigma_dmu( energy, muCutoff, accuracy = accuracy ).integrate( ) ) ] )
        _tester = tester( self.dSigma_dmu, muCutoff, accuracy, accuracy * crossSection[-1][1] )
        crossSection = fudgemathModule.thickenXYList( crossSection, _tester, biSectionMax = 10 )

        xys2d = angularModule.XYs2d( axes = angularModule.defaultAxes( self.domainUnit ) )

        for energy, xSec in crossSection :
            xys1d = self.dSigma_dmu( energy, muCutoff, accuracy = accuracy )
            xys2d.append( angularModule.XYs1d( data = xys1d / xSec, axes = xys2d.axes, value = energy ) )

        axes = crossSectionModule.defaultAxes( self.domainUnit )
        axes[0].unit = xys1d.axes[0].unit
        crossSection = crossSectionModule.XYs1d( data = crossSection, axes = axes, label = style.label )

        return( crossSection, xys2d )

    def processMultiGroup( self, style, tempInfo, indent ) :

        print '    processMultiGroup not implemented for distribution form %s.' % self.moniker
        return( None )

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        attrs = ' label="%s" pid="%s" productFrame="%s"' % (self.label, self.pid, self.productFrame)
        if self.identicalParticles: attrs += ' identicalParticles="true"'
        xml = ['%s<%s%s>' % (indent,self.moniker,attrs)]
        xml += self.RutherfordScattering.toXMLList( indent2, **kwargs )
        if self.data is not None:
            xml += self.data.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        data = None
        for child in element:
            if child.tag == RutherfordScattering.moniker:
                RutherfordTerm = RutherfordScattering.parseXMLNode( child, xPath, linkData )
            elif child.tag == nuclearPlusInterference.moniker:
                data = nuclearPlusInterference.parseXMLNode( child, xPath, linkData )
            elif child.tag == nuclearAmplitudeExpansion.moniker:
                data = nuclearAmplitudeExpansion.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Encountered unexpected element '%s' in %s" % (child.tag, element.tag))
        Coul = form( element.get('label'), element.get('pid'), data, productFrame=element.get('productFrame'),
            RutherfordTerm=RutherfordTerm, identicalParticles=element.get('identicalParticles','')=='true')
        xPath.pop()
        return Coul

class CoulombDepositionNotSupported( Exception ):
    """
    Custom Exception, returned when calculateAverageProductData() is called for
    nuclearAmplitudeExpansion or nuclearPlusInterference
    """
    pass
