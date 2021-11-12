# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Defines classes used in both resolved and unresolved regions
"""

import fractions

from fudge import suites as suitesModule
from fudge.channelData import Q as QModule
from pqu import PQU as PQUModule
from .scatteringRadius import scatteringRadius, hardSphereRadius

from xData import ancestry as ancestryModule, table as tableModule, link as linkModule

__metaclass__ = type

class spin( fractions.Fraction ):
    """
    Store spins for a collection of resonances. Check denominator (must be integer or half-integer)
    """

    def __new__( cls, *args ):
        self = fractions.Fraction.__new__(cls, *args)
        if self.denominator not in (1,2):
            raise ValueError("Illegal spin '%s': must be integer or half-integer" % self)
        return self

    @property
    def value( self ):

        return float(self)

class parity:
    """
    Store parity for a collection of resonances. Allowed values for the parity are +1 or -1.
    """

    def __init__(self, parity):
        self.value = int(parity)
        if self.value not in (1,-1):
            raise ValueError("%d is not a legal value for parity!" % self.value)

    def __str__(self):

        return str(self.value)

    def __int__(self):

        return self.value

    def __copy__( self ) :

        return( parity( self.value ) )

class resonanceReactions( suitesModule.suite ):
    """
    Stores a list of resonanceReaction
    """

    moniker = 'resonanceReactions'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [resonanceReaction])

    def check( self, info ):
        from fudge import warning
        warnings = []
        for c in self:
            warningList = c.check(info)
            if warningList:
                warnings.append(warning.context(str(c.moniker)+' '+str(c.label),warningList))
        return warnings

class resonanceReaction( ancestryModule.ancestry):
    """
    Describes one reaction channel that opens up in the resonance region. In an R-Matrix section,
    all open reaction channels should be described in the list of resonanceReaction elements
    """

    moniker = 'resonanceReaction'

    fission = 'fission'                 # special tokens to support fission reactions
    fissionProduct = 'fissionProduct'

    def __init__( self, label, reactionLink, ejectile, Q=None, scatteringRadius=None, hardSphereRadius=None,
                  computePenetrability=True, computeShiftFactor=False, eliminated=False):
        """
        :param label: unique label for this reaction
        :param reactionLink: href pointing to corresponding reaction in reactionSuite/reactions
        :param ejectile: id for the light particle emitted by this reaction
        :param Q: optional, overloads the linked reaction Q-value
        :param scatteringRadius: optional, overrides resonances.scatteringRadius
        :param hardSphereRadius: optional, overrides resonances.scatteringRadius
        :param computePenetrability: boolean, default=True. If False, use P=1
        :param computeShiftFactor: boolean, default=False
        :param eliminated: boolean, default=False
        """

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.reactionLink = reactionLink
        self.__ejectile = ejectile
        self.__residual = None
        self.Q = Q
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        if( hardSphereRadius is not None ) : self.hardSphereRadius.setAncestor( self )
        self.computePenetrability = computePenetrability
        self.computeShiftFactor = computeShiftFactor
        self.eliminated = eliminated

    @property
    def ejectile(self):
        return self.__ejectile

    @property
    def residual(self):
        if self.__residual is None:
            if self.ejectile == resonanceReaction.fission: return resonanceReaction.fissionProduct

            products = set([p.pid for p in self.reactionLink.link.outputChannel.products])
            products.remove( self.ejectile )
            if len(products) != 1: raise ValueError("Cannot compute resonanceReaction residual!")
            self.__residual = products.pop()
        return self.__residual

    @property
    def scatteringRadius(self):
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            from .resonances import resonances
            return self.findClassInAncestry(resonances).scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a scatteringRadius instance."""

        self.__scatteringRadius = value
        if value is not None:
            if not isinstance(value, scatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            self.__scatteringRadius.setAncestor(self)

    def check( self, info ):
        from fudge import warning
        warnings = []

        # check the reaction link
        theLinkTarget = None
        try:
            theLinkTarget = self.reactionLink.follow(self.getRootAncestor())
            if theLinkTarget is None: warnings.append(warning.unresolvedLink(self.reactionLink))
        except:
            warnings.append(warning.unresolvedLink(self.reactionLink))

        # check the radii
        for thing in [self.scatteringRadius, self.hardSphereRadius]:
            if thing is  None: continue
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.context(thing.moniker,warningList))
        return warnings

    def convertUnits( self, unitMap ):
        for child in ('Q','scatteringRadius','hardSphereRadius'):
            if getattr(self,child) is not None:
                getattr(self,child).convertUnits(unitMap)

    def isFission( self ):

        return self.reactionLink.link.isFission()

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        attrstring = ''
        if self.ejectile is not None: attrstring += ' ejectile="%s"' % self.ejectile
        if not self.computePenetrability: attrstring += ' computePenetrability="false"'
        if self.computeShiftFactor: attrstring += ' computeShiftFactor="true"'
        if self.eliminated: attrstring += ' eliminated="true"'
        xmlString = ['%s<%s label="%s"%s>' % (indent, self.moniker, self.label, attrstring) ]
        xmlString += self.reactionLink.toXMLList( indent=indent2, **kwargs )
        if self.Q is not None: xmlString += self.Q.toXMLList( indent=indent2, **kwargs )
        if self.__scatteringRadius is not None: xmlString += self.__scatteringRadius.toXMLList( indent=indent2, **kwargs)
        if self.hardSphereRadius is not None: xmlString += self.hardSphereRadius.toXMLList( indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        reactionLink = linkModule.link.parseXMLNode( element.find('link'), xPath, linkData )
        Qval, scatRad, hsRad = None, None, None
        if element.find( QModule.component.moniker ):
            Qval = QModule.component()      # FIXME suite.parseXMLNode is not a static method
            Qval.parseXMLNode( element.find( QModule.component.moniker ), xPath, linkData )
        if element.find( scatteringRadius.moniker ):
            scatRad = scatteringRadius.parseXMLNode( element.find( scatteringRadius.moniker ), xPath, linkData )
        if element.find( hardSphereRadius.moniker ):
            hsRad = hardSphereRadius.parseXMLNode( element.find( hardSphereRadius.moniker), xPath, linkData)
        tmp = resonanceReaction( element.get('label'), reactionLink=reactionLink, ejectile=element.get('ejectile'),
                                 Q = Qval, scatteringRadius=scatRad, hardSphereRadius=hsRad,
                                 computePenetrability=getBool( element.get('computePenetrability','true') ),
                                 computeShiftFactor=getBool( element.get('computeShiftFactor','false') ),
                                 eliminated=getBool( element.get('eliminated','false') ) )
        xPath.pop()
        return tmp

class resonanceParameters( ancestryModule.ancestry ):
    """
    Light-weight wrapper around a table.
    """

    moniker = 'resonanceParameters'

    def __init__(self, table):
        ancestryModule.ancestry.__init__(self)
        self.table = table
        self.table.setAncestor(self)

    def check( self, info ):
        warnings=[]
        return warnings

    def convertUnits( self, unitMap ):

        self.table.convertUnits( unitMap )

    def toXMLList(self, indent = '', **kwargs):

        indent2 = indent+'  '
        xmlList = ['%s<%s>' % (indent, self.moniker)]
        xmlList += self.table.toXMLList(indent2, **kwargs)
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        rps = resonanceParameters( tableModule.table.parseXMLNode(
            element.find(tableModule.table.moniker), xPath, linkData ) )
        xPath.pop()
        return rps

class energyIntervals(ancestryModule.ancestry):
    """ Resonance region may be broken up into multiple energy intervals (deprecated) """

    moniker = 'energyIntervals'
    def __init__(self, label):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.__intervals = []

    def __len__(self): return len(self.__intervals)

    def __getitem__(self, item): return self.__intervals[item]

    def append(self, item):
        self.__intervals.append(item)
        item.setAncestor(self, attribute='index')

    def check( self, info ):
        from fudge import warning
        warnings = []
        for idx, interval in enumerate(self):
            info['energyIntervalIndex'] = idx
            warningList = interval.check(info)
            if warningList:
                warnings.append(warning.context('%s[@index="%d"' % (interval.moniker, interval.index) ,warningList))
        del info['energyIntervalIndex']
        return warnings

    def convertUnits(self, unitMap):
        for interval in self:
            interval.convertUnits(unitMap)

    @property
    def useForSelfShieldingOnly(self):
        for interval in self:
            if interval.evaluated.useForSelfShieldingOnly: return True
        return False

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlString = [ '%s<%s label="%s">' % (indent, self.moniker, self.label) ]
        for interval in self:
            xmlString += interval.toXMLList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append( element.tag )
        EIs = cls(element.get('label'))
        for child in element:
            EIs.append( energyInterval.parseXMLNode(child,xPath,linkData) )
        xPath.pop()
        return EIs

class energyInterval(ancestryModule.ancestry):
    """ single energy interval, for use inside energyIntervals """

    moniker = 'energyInterval'
    def __init__(self, index, data, domainMin, domainMax, domainUnit):
        ancestryModule.ancestry.__init__(self)
        self.index = index
        data.setAncestor(self)
        self.evaluated = data
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

    def check( self, info ):
        from fudge import warning
        warnings = []
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.context( self.evaluated.moniker, warningList ))
        return warnings

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQUModule.PQU(1, self.domainUnit).getValueAs(newUnit)
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit
        self.evaluated.convertUnits(unitMap)

    def toString(self, simpleString = False):
        return ("%s resonances, %s to %s. Contains %i resonances" %
                (self.evaluated, self.domainMin, self.domainMax, len(self.evaluated) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = ['%s<%s index="%s" domainMin="%s" domainMax="%s" domainUnit="%s">' % ( indent, self.moniker, self.index, 
                PQUModule.floatToShortestString( self.domainMin, 12 ), PQUModule.floatToShortestString( self.domainMax, 12 ), self.domainUnit ) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @property
    def BreitWigner(self):
        from .resolved import BreitWigner
        if isinstance(self.evaluated, BreitWigner): return self.evaluated
        return None

    @property
    def RMatrix(self):
        from .resolved import RMatrix
        if isinstance(self.evaluated, RMatrix): return self.evaluated
        return None

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        from .resolved import BreitWigner, RMatrix
        from .unresolved import tabulatedWidths

        xPath.append('%s[@label="%s"]' % (cls.moniker, element.get('label')))

        formClass = {
            BreitWigner.moniker: BreitWigner,
            RMatrix.moniker: RMatrix,
            tabulatedWidths.moniker: tabulatedWidths,
        }.get( element[0].tag )
        if formClass is None: raise Exception("unknown unresolved resonance form '%s'!" % element[0].tag)
        data = formClass.parseXMLNode(element[0], xPath, linkData)
        EI = cls( data=data, **getAttrs(element) )

        xPath.pop()
        return EI

# helper functions for reading in from xml:
def getBool( value ):
    return {'true':True, '1':True, 'false':False, '0':False}[value]

def floatOrint( value ):
    if float( value ).is_integer(): return int( value )
    return float( value )

def getAttrs(element, exclude=(), required=()):
    """
    Convert attributes to proper type (float, bool, int, etc), returning a dictionary.
    Anything in the 'exclude' list is omitted, anything in the 'required' list is automatically set to False
    if not present.
    """
    conversionTable = {
            'domainMin':float, 'domainMax':float, 'value':float, 'channelSpin':spin,
            'calculateChannelRadius':getBool, 'useForSelfShieldingOnly':getBool,
            'computeShiftFactor':getBool,'computePenetrability':getBool, 'columnIndex':int,
            'ENDF_MT':int, 'index':int, 'L':int, 'reducedWidthAmplitudes':getBool,
            'spin':spin, 'parity':parity, 'boundaryConditionValue':float,
            'supportsAngularReconstruction':getBool,
            }
    attrs = dict( element.items() )
    for key in attrs.keys():
        if key in exclude: attrs.pop(key)
        elif key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
    for val in required:
        if val not in attrs: attrs[val] = False
    return attrs

