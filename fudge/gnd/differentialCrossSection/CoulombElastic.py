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

import abc
import numpy

from fudge.gnd import abstractClasses as abstractClassesModule
from xData import ancestry as ancestryModule
from xData import standards as standardsModule

from PoPs import misc as miscPoPsModule

from fudge.gnd.reactionData import crossSection as crossSectionModule
from fudge.gnd.productData.distributions import angular as angularModule

__metaclass__ = type

class form( abstractClassesModule.form ):

    moniker = "CoulombElastic"

    def __init__( self, label, data, productFrame = standardsModule.frames.centerOfMassToken,
                  RutherfordTerm = None, identicalParticles = False ):

        abstractClassesModule.form.__init__( self )
        self.label = label
        self.__productFrame = productFrame
        self.__identicalParticles = identicalParticles

        if RutherfordTerm is not None:
            self.RutherfordScattering = RutherfordTerm
        else:
            self.RutherfordScattering = RutherfordScattering()
        self.RutherfordScattering.setAncestor( self )
        self.data = data

    @property
    def productFrame(self): return self.__productFrame

    @property
    def identicalParticles(self): return self.__identicalParticles

    @property
    def data(self): return self.__data

    @data.setter
    def data(self, value):
        if not isinstance(value, (NuclearPlusCoulombInterference, CoulombExpansion, type(None))):
            raise TypeError("Data type '%s' not allowed in %s" % (type(value),self.moniker))
        if value is not None: value.setAncestor(self)
        self.__data = value

    # two additional ways to access data (meant be user-friendly)
    @property
    def NuclearPlusCoulombInterference(self):
        if isinstance( self.data, NuclearPlusCoulombInterference ): return self.data
        return None

    @property
    def CoulombExpansion(self):
        if isinstance( self.data, CoulombExpansion ): return self.data
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

        from fudge.gnd import warning
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

    def evaluate(self, E, mu):
        """
        Compute the cross section at E / mu, including Coulomb, nuclear and interference terms
        :param E: incident energy (lab frame).
        :param mu: scattering angle cosine (in center-of-mass)
        :return:
        """
        RS = self.RutherfordScattering.evaluate(E, mu)
        NS = 0
        if self.data is not None:
            NS = self.data.evaluate(E, mu)
        return RS + NS

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processMultiGroup( self, style, tempInfo, indent ) :

        print '    processMultiGroup not implemented for distribution form %s.' % self.moniker
        return( None )

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        attrs = ' label="%s" productFrame="%s"' % (self.label, self.productFrame)
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
            elif child.tag == NuclearPlusCoulombInterference.moniker:
                data = NuclearPlusCoulombInterference.parseXMLNode( child, xPath, linkData )
            elif child.tag == CoulombExpansion.moniker:
                data = CoulombExpansion.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Encountered unexpected element '%s' in %s" % (child.tag, element.tag))
        Coul = form( element.get('label'), data, productFrame=element.get('productFrame'),
            RutherfordTerm=RutherfordTerm, identicalParticles=element.get('identicalParticles','')=='true')
        xPath.pop()
        return Coul


class RutherfordScattering( ancestryModule.ancestry ):
    """
    Stores the pure-Coulomb contribution to charged-particle elastic scattering
    """

    moniker = "RutherfordScattering"

    def __init__(self, domainMin=None, domainMax=None, domainUnit=None ):

        ancestryModule.ancestry.__init__(self)

        self.__initialized = False      # initialize when evaluate() is first called
        self.identicalParticles = None
        self.spin = None
        self.etaCoefficient = None
        self.kCoefficient = None
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

    def initialize(self):
        """
        Pre-compute some factors used to calculate the Rutherford cross section:
        """

        reactionSuite = self.getRootAncestor( )
        self.identicalParticles = self.findAttributeInAncestry('identicalParticles')

        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        target = reactionSuite.PoPs[reactionSuite.target]
        Z1 = miscPoPsModule.ZAInfo( projectile )[0]
        Z2 = miscPoPsModule.ZAInfo( target )[0]
# BRB6 hardwired
        mass1 = projectile.mass[0].float( 'amu' )
        mass2 = target.mass[0].float( 'amu' )
        if self.identicalParticles:
            self.spin = projectile.spin[0].value

        # FIXME check against CODATA. Should be defined in PQU or fudge.core instead?
        # also, equations below are taken from ENDF manual (section 6.2.6)... need to be checked!
        amu = 931.494013e+6     # in eV/c**2
        hbar = 6.58211889e-16   # in eV*s
        c = 2.99792458e+8       # m/s
        alpha = 1/137.035999    # unitless

        if mass1 > mass2:
            mass1, mass2 = mass2, mass1
        self.etaCoefficient = Z1*Z2 * numpy.sqrt(alpha**2*amu/2*mass1)
        A = mass2 / mass1
        self.kCoefficient = (A / (A+1) ) * numpy.sqrt((2/hbar**2) * (amu/c**2) * mass1) * 1e-14

        self.__initialized = True

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            from pqu import PQU
            newUnit = unitMap[ self.domainUnit ]
            factor = PQU.PQU( 1., self.domainUnit ).getValueAs( newUnit )
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit

    def evaluate(self, E, mu):
        if not self.__initialized: self.initialize()

        eta = self.etaCoefficient / numpy.sqrt(E)
        k = self.kCoefficient * numpy.sqrt(E)
        if self.identicalParticles:
            term1 = (2*eta**2) / (k**2 * (1-mu**2))
            term2 = (1+mu**2) / (1-mu**2) + -1**(2*self.spin) / (2*self.spin + 1) * numpy.cos(
                eta * numpy.log( (1+mu) / (1-mu) ) )
            return term1 * term2
        else:
            return eta**2 / (k**2 * (1-mu)**2 )

    def toXMLList( self, indent='', **kwargs ):

        attrStr = ''
        for attribute in ('domainMin','domainMax','domainUnit'):
            if getattr(self,attribute) is not None: attrStr += ' %s="%s"' % (attribute, getattr(self,attribute))
        return ['%s<%s%s/>' % (indent, self.moniker, attrStr)]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        domainMin, domainMax = element.get('domainMin'), element.get('domainMax')
        if domainMin is not None:
            domainMin, domainMax = map(float, (domainMin, domainMax) )
        return RutherfordScattering( domainMin, domainMax, element.get('domainUnit') )


class term( ancestryModule.ancestry ):
    """
    In addition to the Rutherford scattering term, charged-particle elastic sections
    may have additional terms: 'effectiveCrossSection', 'effectiveDistribution',
    'nuclearTerm', etc. Each is just a context layer that contains the actual data.
    This serves as the base class for all terms.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def allowedDataForms(self): pass

    def __init__(self, data):
        ancestryModule.ancestry.__init__(self)
        self.data = data

    @property
    def data(self): return self.__data

    @data.setter
    def data(self, value):
        if not isinstance(value, self.allowedDataForms):
            raise TypeError("'%s' cannot contain data of type '%s'" % (self.moniker, type(value)))
        value.setAncestor( self )
        self.__data = value

    def convertUnits( self, unitMap ):

        self.data.convertUnits( unitMap )

    def toXMLList(self, indent='', **kwargs):

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXMLList(indent+'  ', **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        child = element[0]
        parseClass = None
        for allowedForm in cls.allowedDataForms:
            if child.tag == allowedForm.moniker:
                parseClass = allowedForm
                break
        if parseClass is None:
            raise TypeError("Unexpected element '%s' encountered inside %s" % (child.tag, cls.moniker))
        data = parseClass.parseXMLNode( child, xPath, linkData )
        term_ = cls( data )
        xPath.pop()
        return term_


class effectiveCrossSection( term ):

    moniker = 'effectiveCrossSection'

    allowedDataForms = (crossSectionModule.XYs1d, crossSectionModule.regions1d)

class effectiveDistribution( term ):

    moniker = 'effectiveDistribution'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)


class NuclearPlusCoulombInterference( ancestryModule.ancestry ):

    moniker = "NuclearPlusCoulombInterference"

    def __init__( self, muCutoff, effectiveCrossSection = None, effectiveDistribution = None ):

        ancestryModule.ancestry.__init__( self )
        self.__muCutoff = muCutoff
        self.effectiveCrossSection = effectiveCrossSection
        self.effectiveDistribution = effectiveDistribution

    @property
    def muCutoff(self): return self.__muCutoff

    @property
    def effectiveCrossSection( self ):
        return self.__effectiveCrossSection

    @effectiveCrossSection.setter
    def effectiveCrossSection(self, value):
        if not isinstance(value, effectiveCrossSection):
            raise TypeError( '%s cross section must be instance of %s' % (self.moniker, effectiveCrossSection.moniker))
        value.setAncestor( self )
        self.__effectiveCrossSection = value

    @property
    def effectiveDistribution( self ):
        return self.__effectiveDistribution

    @effectiveDistribution.setter
    def effectiveDistribution(self, value):
        if not isinstance(value, effectiveDistribution):
            raise TypeError( '%s distribution must be instance of %s' % (self.moniker, effectiveDistribution.moniker))
        value.setAncestor( self )
        self.__effectiveDistribution = value

    @property
    def domainMin(self): return self.effectiveCrossSection.data.domainMin

    @property
    def domainMax(self): return self.effectiveCrossSection.data.domainMax

    @property
    def domainUnit(self): return self.effectiveCrossSection.data.domainUnit

    def check( self, info ):
        """
        Things to check: mu_cutoff should be between 0.99 and 1.0,
        domain of effective xsc / distribution should match,
        if identical particles distribution should only be for 0->muCutoff, otherwise -1->muCutoff
        effectiveDistribution should be normalized
        If the cross section goes negative, also check that the pure Coulomb portion is large enough to 'win'
        :return:
        """
        return []
        raise NotImplementedError

    def convertUnits( self, unitMap ):

        self.effectiveCrossSection.convertUnits( unitMap )
        self.effectiveDistribution.convertUnits( unitMap )

    def evaluate(self, E, mu):
        """
        :param E: incident energy (or numpy array of incident energies)
        :param mu: scattering angle cosine
        :return: differential cross section at E,mu (in b)
        """
        return self.effectiveCrossSection.data.evaluate( E ) * self.effectiveDistribution.data.evaluate( E, mu )

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s muCutoff="%s">' % (indent, self.moniker, self.muCutoff) ]
        xml += self.effectiveCrossSection.toXMLList( indent2, **kwargs )
        xml += self.effectiveDistribution.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag == effectiveCrossSection.moniker:
                effXsc = effectiveCrossSection.parseXMLNode( child, xPath, linkData )
            elif child.tag == effectiveDistribution.moniker:
                effDist = effectiveDistribution.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        NPCI = NuclearPlusCoulombInterference(
                element.get('muCutoff'), effectiveCrossSection=effXsc, effectiveDistribution=effDist )

        xPath.pop()
        return NPCI


class nuclearTerm( term ):

    moniker = 'nuclearTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class realInterferenceTerm( term ):

    moniker = 'realInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class imaginaryInterferenceTerm( term ):

    moniker = 'imaginaryInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)


class CoulombExpansion( ancestryModule.ancestry ):

    moniker = 'CoulombExpansion'

    def __init__( self, nuclearTerm = None, realInterference = None, imaginaryInterference = None ):

        ancestryModule.ancestry.__init__( self )
        self.nuclearTerm = nuclearTerm
        self.realInterferenceTerm = realInterference
        self.imaginaryInterferenceTerm = imaginaryInterference

    @property
    def nuclearTerm( self ):
        return self.__nuclear

    @nuclearTerm.setter
    def nuclearTerm(self, value):
        if not isinstance(value, nuclearTerm):
            raise TypeError( '%s nuclear term must be instance of %s' % (self.moniker, nuclearTerm.moniker))
        value.setAncestor( self )
        self.__nuclear = value

    @property
    def realInterferenceTerm( self ):
        return self.__realInterference

    @realInterferenceTerm.setter
    def realInterferenceTerm(self, value):
        if not isinstance(value, realInterferenceTerm):
            raise TypeError( '%s real interference term must be instance of %s' %
                             (self.moniker, realInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__realInterference = value

    @property
    def imaginaryInterferenceTerm( self ):
        return self.__imaginaryInterference

    @imaginaryInterferenceTerm.setter
    def imaginaryInterferenceTerm(self, value):
        if not isinstance(value, imaginaryInterferenceTerm):
            raise TypeError( '%s imaginary interference term must be instance of %s' %
                             (self.moniker, imaginaryInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__imaginaryInterference = value

    @property
    def domainMin(self): return self.nuclearTerm.data.domainMin

    @property
    def domainMax(self): return self.nuclearTerm.data.domainMax

    @property
    def domainUnit(self): return self.nuclearTerm.data.domainUnit

    def check( self, info ):
        """
        Check that incident energy domain of all three terms are the same
        FIXME: I'm not sure whether each term is meant to be normalized  (I see L=0 terms != 0)
        There may be more checks we can perform depending on whether particles are identical
        :return:
        """
        return []
        raise NotImplementedError()

    def convertUnits( self, unitMap ):

        self.nuclearTerm.convertUnits( unitMap )
        self.realInterferenceTerm.convertUnits( unitMap )
        self.imaginaryInterferenceTerm.convertUnits( unitMap )

    def evaluate(self, E, mu):
        """
        :param E: incident energy (or numpy array of incident energies)
        :param mu: scattering angle cosine
        :return: differential cross section at E,mu (in b)
        """
        raise NotImplementedError()

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s>' % (indent, self.moniker) ]
        xml += self.nuclearTerm.toXMLList( indent2, **kwargs )
        xml += self.realInterferenceTerm.toXMLList( indent2, **kwargs )
        xml += self.imaginaryInterferenceTerm.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag == nuclearTerm.moniker:
                nuclear = nuclearTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == realInterferenceTerm.moniker:
                real = realInterferenceTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == imaginaryInterferenceTerm.moniker:
                imag = imaginaryInterferenceTerm.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        CoulExp = CoulombExpansion( nuclearTerm=nuclear, realInterference=real, imaginaryInterference=imag )
        xPath.pop()
        return CoulExp


class CoulombDepositionNotSupported( Exception ):
    """ Custom Exception, returned when calculatedDepositionData() called for CoulombExpansion or NuclearPlusCoulomb """
    pass
