# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Containers for unresolved resonance parameters
"""

import fractions
import abc

from PoPs import database as PoPsDatabaseModule
from pqu import PQU as PQUModule
from xData import ancestry as ancestryModule, XYs as XYsModule, regions as regionsModule
from xData.Documentation import documentation as documentationModule

from fudge import suites as suitesModule, abstractClasses as abstractClassesModule

from .common import energyIntervals, resonanceReactions, floatOrint, getAttrs
from .scatteringRadius import scatteringRadius

__metaclass__ = type

class unresolved( abstractClassesModule.component ):

    moniker = 'unresolved'

    def __init__(self, domainMin, domainMax, domainUnit):
        abstractClassesModule.component.__init__(self, allowedClasses=(tabulatedWidths, energyIntervals))
        self.__domainMin = domainMin
        self.__domainMax = domainMax
        self.__domainUnit = domainUnit

    @property
    def domainMin( self ) :

        return( self.__domainMin )

    @property
    def domainMax( self ) :

        return( self.__domainMax )

    @property
    def domainUnit( self ) :

        return( self.__domainUnit )

    def toString( self, simpleString = False ):
        return ("Unresolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        from fudge import warning
        warnings = []
        for L in self.evaluated.Ls:
            for J in L.Js:
                elist = J.levelSpacing.data.convertAxisToUnit(1, self.__domainUnit).domainGrid
                if elist[0] > self.__domainMin or elist[-1] < self.__domainMax:
                    warnings.append( warning.URRdomainMismatch( L.L, J.J, J ) )
                missingPoints = [i1 for i1 in range(1,len(elist)) if elist[i1] > 3*elist[i1-1]]
                for idx in missingPoints:
                    warnings.append( warning.URRinsufficientEnergyGrid( L.L, J.J,
                        PQUModule.PQU(elist[idx-1],self.__domainUnit), PQUModule.PQU(elist[idx],self.__domainUnit), J ) )
        return warnings

    def convertUnits( self, unitMap ):

        if self.__domainUnit in unitMap:
            newUnit = unitMap[self.__domainUnit]
            factor = PQUModule.PQU(1, self.__domainUnit).getValueAs(newUnit)
            self.__domainMin *= factor
            self.__domainMax *= factor
            self.__domainUnit = newUnit
        for form in self:
            form.convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s domainMin="%s" domainMax="%s" domainUnit="%s">' % ( indent, self.moniker, 
                PQUModule.floatToShortestString( self.__domainMin, 12 ), PQUModule.floatToShortestString( self.__domainMax, 12 ), self.__domainUnit ) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )

        URR = cls( **getAttrs(element) )
        for child in element:
            formClass = {
                    tabulatedWidths.moniker: tabulatedWidths,
                    energyIntervals.moniker: energyIntervals,
                    }.get( child.tag )
            if formClass is None: raise Exception("unknown unresolved resonance form '%s'!" % child.tag)
            URR.add(formClass.parseXMLNode( child, xPath, linkData ))

        xPath.pop()
        return URR

class tabulatedWidths(ancestryModule.ancestry):

    moniker = 'tabulatedWidths'

    def __init__(self, label, approximation, resonanceReactions_, Ls_, scatteringRadius=None, PoPs=None,
                 useForSelfShieldingOnly=False):
        """
        Contains unresolved resonance parameters

        :param label: unique label corresponding to a style (i.e. 'eval')
        :param approximation: currently only 'SingleLevelBreitWigner' is supported
        :param resonanceReactions_: resonanceReactions instance
        :param Ls_: URR_Lsections instance
        :param scatteringRadius: optional scatteringRadius instance
        :param PoPs: optional PoPs database
        :param useForSelfShieldingOnly: boolean
        """

        ancestryModule.ancestry.__init__(self)

        self.label = label
        self.approximation = approximation
        self.resonanceReactions = resonanceReactions_
        self.scatteringRadius = scatteringRadius
        self.PoPs = PoPs
        self.Ls = Ls_
        self.useForSelfShieldingOnly = useForSelfShieldingOnly

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

    def __len__(self):
        return len(self.Ls)

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    @property
    def resonanceReactions(self):
        return self.__resonanceReactions

    @resonanceReactions.setter
    def resonanceReactions(self, value):
        if not isinstance(value, resonanceReactions):
            raise TypeError("Must be a resonanceReactions instance")
        value.setAncestor(self)
        self.__resonanceReactions = value

    @property
    def scatteringRadius( self ):
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            from .resonances import resonances
            return self.findClassInAncestry(resonances).scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """ Can be set to None or to a scatteringRadius instance. """
        self.__scatteringRadius = value
        if value is not None:
            if not isinstance(value, scatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            self.__scatteringRadius.setAncestor(self)

    @property
    def Ls(self):
        return self.__Ls

    @Ls.setter
    def Ls(self, value):
        if not isinstance(value, Lsections):
            raise TypeError("Must be a URR_Lsections instance")
        value.setAncestor(self)
        self.__Ls = value

    @property
    def PoPs(self):
        if self.__PoPs is not None:
            return self.__PoPs
        else:
            return self.getRootAncestor().PoPs

    @PoPs.setter
    def PoPs(self, value):
        """ Can be set to None or to a PoPs database instance """
        self.__PoPs = value
        if value is not None:
            if not isinstance(value, PoPsDatabaseModule.database):
                raise TypeError("PoPs can't be set to type '%s'" % type(value))
            self.__PoPs.setAncestor(self)

    def convertUnits( self, unitMap ):
        if self.__scatteringRadius is not None:
            self.__scatteringRadius.convertUnits(unitMap)
        for lval in self.Ls:
            for jval in lval.Js:
                jval.convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ' label="%s" approximation="%s"' % (self.label, self.approximation)
        if self.useForSelfShieldingOnly: attrs += ' useForSelfShieldingOnly="true"'

        xml = ['%s<%s%s>' % (indent, self.moniker, attrs)]

        xml += self.__documentation.toXMLList( indent = indent2, **kwargs )

        if self.__PoPs:
            xml += self.__PoPs.toXMLList( indent2, **kwargs )

        xml += self.resonanceReactions.toXMLList( indent2, **kwargs )
        if self.__scatteringRadius:
            xml += self.__scatteringRadius.toXMLList( indent2, **kwargs )
        xml += self.Ls.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker

        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )

        linkData['conversionTable'] = {'index':int}
        resonanceReactions_ = resonanceReactions()
        resonanceReactions_.parseXMLNode( element.find(resonanceReactions.moniker), xPath, linkData )
        radius = element.find( scatteringRadius.moniker )
        if radius is not None:
            radius = scatteringRadius.parseXMLNode( radius, xPath, linkData )
        pops = element.find( PoPsDatabaseModule.database.moniker )
        if pops is not None:
            pops = PoPsDatabaseModule.database.parseXMLNodeAsClass( pops, xPath, linkData )
        Ls_ = Lsections()
        Ls_.parseXMLNode( element.find( Ls_.moniker ), xPath, linkData )

        result = cls( element.get('label'), element.get('approximation'),
            resonanceReactions_, Ls_=Ls_, scatteringRadius=radius, PoPs=pops,
            useForSelfShieldingOnly=element.get('useForSelfShieldingOnly') == 'true'
        )
        del linkData['conversionTable']

        documentation = element.find( documentationModule.Documentation.moniker )
        if( documentation is not None ) : result.documentation.parseNode( documentation, xPath, linkData )

        xPath.pop()
        return result

class interpolationTable(ancestryModule.ancestry):
    """
    Base class for unresolved level spacing and widths
    Data are stored as XYs1d or regions1d.
    """

    __metaclass__ = abc.ABCMeta

    allowedSubclasses = (XYsModule.XYs1d, regionsModule.regions1d)

    def __init__( self, data ):
        ancestryModule.ancestry.__init__(self)
        self.data = data

    @property
    def data( self ):
        return self.__data

    @data.setter
    def data( self, value ):
        if not isinstance(value, self.allowedSubclasses):
            raise TypeError("levelSpacing must be an XYs1d or regions1d instance")
        value.setAncestor(self)
        self.__data = value

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        child, = element[:]
        data = None
        for subclass in cls.allowedSubclasses:
            if child.tag == subclass.moniker:
                data = subclass.parseXMLNode(child, xPath, linkData)
                break
        result = cls(data)
        xPath.pop()
        return result

class Lsections(suitesModule.suite):
    """
    tabulated widths contains a list of 'L' sections, each of which contains a list of 'J' sections
    """

    moniker = "Ls"

    def __init__( self ):
        suitesModule.suite.__init__(self, (Lsection,))

class Lsection(ancestryModule.ancestry):
    """ unresolved average widths, grouped by L. Contains list of J-values: """
    moniker = 'L'
    def __init__(self, label, L, Js):
        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.L = L
        self.Js = Js

    @property
    def value(self): return self.L

    @property
    def Js(self):
        return self.__Js

    @Js.setter
    def Js(self, value):
        if not isinstance(value, (Jsections)):
            raise TypeError("Must be URR_Jsections instance")
        value.setAncestor(self)
        self.__Js = value

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = ['%s<%s label="%s" value="%d">' % (indent, self.moniker, self.label, self.value)]
        xml += self.Js.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )
        Js = Jsections()
        Js.parseXMLNode( element.find(Js.moniker), xPath, linkData )
        result = cls( element.get('label'), int(element.get('value')), Js )
        xPath.pop()
        return result

class Jsections(suitesModule.suite):

    moniker = "Js"

    def __init__( self ):
        suitesModule.suite.__init__(self, (Jsection,))

class Jsection(ancestryModule.ancestry):
    """
    Unresolved average parameters for a specific L/J.
    Contains interpolation tables for the level spacing and widths, plus degrees of freedom for each open
    channel  (degrees of freedom are typically 1 or 2 indicating how many channel spins contribute)
    """
    moniker = 'J'

    def __init__(self, label, J, levelSpacing_, widths_):
        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.J = J
        self.levelSpacing = levelSpacing_
        self.widths = widths_

    @property
    def value(self): return str(self.J)

    @property
    def levelSpacing(self): return self.__levelSpacing

    @levelSpacing.setter
    def levelSpacing(self, spacing):
        if not isinstance(spacing, levelSpacing):
            raise TypeError("Expected levelSpacing instance, got %s instead" % type(spacing))
        spacing.setAncestor(self)
        self.__levelSpacing = spacing

    @property
    def widths( self ):
        return self.__widths

    @widths.setter
    def widths( self, _widths ):
        if not isinstance(_widths, widths):
            raise TypeError("Expected URR_widths instance, got %s instead" % type(_widths))
        _widths.setAncestor(self)
        self.__widths = _widths

    def convertUnits( self, unitMap ):

        self.levelSpacing.data.convertUnits(unitMap)
        for width in self.widths:
            width.data.convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = ['%s<%s label="%s" value="%s">' % (indent,self.moniker,self.label,self.J)]
        xml += self.levelSpacing.toXMLList(indent2, **kwargs)
        xml += self.widths.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )
        levelSpacing_ = levelSpacing.parseXMLNode( element.find(levelSpacing.moniker), xPath, linkData )
        widths_ = widths()
        widths_.parseXMLNode(element.find(widths.moniker), xPath, linkData)
        Jsec = Jsection(element.get('label'), fractions.Fraction(element.get('value')), levelSpacing_, widths_)
        xPath.pop()
        return Jsec

class levelSpacing(interpolationTable):
    """
    Contains the average level spacing, stored in XYs1d or regions1d.
    """

    moniker = "levelSpacing"

class widths(suitesModule.suite):
    """
    Contains all average channel widths for an L/J combination.
    """

    moniker = "widths"

    def __init__(self):
        suitesModule.suite.__init__(self, allowedClasses=(width,))

class width(interpolationTable):
    """
    Stores the average width and number of degrees of freedom for a single channel

    Degrees of freedom are typically 0, 1 or 2 depending on how many channel spins contribute, but it may be
    non-integer for a reaction that has a threshold somewhere inside the unresolved region.

    Average width is stored as XYs1d or regions1d.
    """

    moniker = "width"

    def __init__(self, label, resonanceReaction, data, degreesOfFreedom = 0):
        interpolationTable.__init__(self, data)
        self.label = label
        self.resonanceReaction = resonanceReaction
        self.degreesOfFreedom = degreesOfFreedom

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attrs = ''
        if self.degreesOfFreedom > 0:
            attrs += ' degreesOfFreedom="%g"' % self.degreesOfFreedom
        xml = ['%s<%s label="%s" resonanceReaction="%s"%s>' % (indent, self.moniker, self.label, self.resonanceReaction, attrs)]
        xml += self.data.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child, = element[:]
        data = None
        for subclass in cls.allowedSubclasses:
            if child.tag == subclass.moniker:
                data = subclass.parseXMLNode(child, xPath, linkData)
                break
        result = cls(element.get('label'), element.get('resonanceReaction'), data, floatOrint(element.get('degreesOfFreedom',0)))
        xPath.pop()
        return result
