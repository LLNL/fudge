# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Containers for resolved resonance parameters
"""

from pqu import PQU as PQUModule
from xData import ancestry as ancestryModule
from xData.Documentation import documentation as documentationModule

from PoPs import database as PoPsDatabaseModule
from fudge import abstractClasses as abstractClassesModule, suites as suitesModule
from fudge.resonances.resonances import resonances
from fudge.resonances import externalRMatrix as externalRMatrixModule
from fudge.resonances.scatteringRadius import scatteringRadius, hardSphereRadius
from fudge.resonances.common import spin, resonanceParameters, getAttrs, resonanceReactions, energyIntervals

__metaclass__ = type


class resolved( abstractClassesModule.component ):
    """ class for resolved resonances """

    moniker = 'resolved'

    def __init__(self, domainMin, domainMax, domainUnit):

        abstractClassesModule.component.__init__(self,
                                                 allowedClasses=(energyIntervals, BreitWigner, RMatrix))
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

    def toString(self, simpleString = False):
        if isinstance(self.evaluated, energyIntervals):
            return ("Resolved region with DEPRECATED multiple regions\n")
        else:
            return ("Resolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        from fudge import warning
        warnings = []
        if isinstance(self.evaluated, energyIntervals):
            warnings.append( warning.RRmultipleRegions() )
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.context(self.evaluated.moniker,warningList))
        return warnings

    def convertUnits( self, unitMap ):

        if self.__domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQUModule.PQU(1, self.domainUnit).getValueAs( newUnit )
            self.__domainMin *= factor
            self.__domainMax *= factor
            self.__domainUnit = newUnit
        for form in self:
            form.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s domainMin="%s" domainMax="%s" domainUnit="%s">' % ( indent, self.moniker,
                PQUModule.floatToShortestString( self.__domainMin, 12 ), PQUModule.floatToShortestString( self.__domainMax, 12 ), self.__domainUnit) ]
        for form in self : xmlString += form.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )

        RRR = cls(**getAttrs(element))
        for child in element:
            formClass = {
                    BreitWigner.moniker: BreitWigner,
                    RMatrix.moniker: RMatrix,
                    energyIntervals.moniker: energyIntervals,
                    }.get( child.tag )
            if formClass is None: raise Exception("unknown resolved resonance form '%s'!" % child.tag)
            RRR.add(formClass.parseXMLNode( child, xPath, linkData ))

        xPath.pop()
        return RRR

class BreitWigner( ancestryModule.ancestry ) :

    moniker = 'BreitWigner'
    singleLevel = 'singleLevel'
    multiLevel = 'multiLevel'
    ancestryMembers = ( 'scatteringRadius', 'resonanceParameters' )
    optAttrList = ( 'calculateChannelRadius', 'computeAngularDistribution', 'LvaluesNeededForConvergence',
                    'useForSelfShieldingOnly')

    def __init__(self, label, approximation, resonanceParameters=None, scatteringRadius=None, PoPs=None, **kwargs):
        """
        Container for resonance parameters using Single-Level or Multi-Level Breit-Wigner approximation

        :param label: corresponds to a style (i.e. 'eval')
        :param approximation: 'singleLevel' or 'multiLevel'
        :param resonanceParameters: optional ResonanceParameters instance
        :param scatteringRadius: optional ScatteringRadius instance
        :param PoPs: optional PoPs.database instance
        :param kwargs:
        """

        ancestryModule.ancestry.__init__( self )

        for attr in self.optAttrList:
            setattr( self, attr, kwargs.get(attr) )
        if self.computeAngularDistribution:
            self.computeAngularDistribution = bool(self.computeAngularDistribution)

        self.label = label
        if approximation not in (BreitWigner.singleLevel, BreitWigner.multiLevel):
            raise TypeError("Unknown approximation '%s' for BreitWigner resonance section" % approximation)
        self.approximation = approximation
        self.resonanceParameters = resonanceParameters or []
        if self.resonanceParameters:
            self.resonanceParameters.setAncestor( self )
        self.scatteringRadius = scatteringRadius
        self.PoPs = PoPs

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    def check( self, info ):
        return _resonance_checker(self, info, [self.scatteringRadius,self.resonanceParameters])

    def convertUnits( self, unitMap ):

        if self.scatteringRadius is not None:
            self.scatteringRadius.convertUnits(unitMap)
        self.resonanceParameters.convertUnits(unitMap)

    def __getitem__(self, idx):
        return self.resonanceParameters.table[idx]

    def __len__(self):
        return len(self.resonanceParameters.table)

    def addResonance(self, resonance):
        """ insert a new resonance in the resonance parameter table """
        #resonance = (energy, J, l, ... )
        self.resonanceParameters.table.addRow( resonance )

    @property
    def scatteringRadius(self):
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            return self.findClassInAncestry(resonances).scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a scatteringRadius instance."""

        self.__scatteringRadius = value
        if value is not None:
            if not isinstance(value, scatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            self.__scatteringRadius.setAncestor(self)

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

    def toXMLList( self, indent = '', **kwargs ):

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        if not self.moniker :
            raise NotImplementedError ("Please use specific formalisms (MLBW, RM, etc) instead of the BaseClass")
        xmlString = '%s<%s label="%s" approximation="%s"' % ( indent, self.moniker, self.label, self.approximation )
        for attr in self.optAttrList:
            if getattr(self,attr,None):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal )
        xmlString = [xmlString+'>']

        xmlString += self.__documentation.toXMLList( indent = indent2, **kwargs )

        if self.__PoPs is not None:
            xmlString += self.__PoPs.toXMLList( indent2, **kwargs )
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXMLList( indent2, **kwargs )
        if self.resonanceParameters:
            xmlString.extend( self.resonanceParameters.toXMLList( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )

        pops = element.find( PoPsDatabaseModule.database.moniker )
        if pops is not None:
            pops = PoPsDatabaseModule.database.parseXMLNodeAsClass( pops, xPath, linkData )

        radius = element.find( scatteringRadius.moniker )
        if radius is not None:
            radius = scatteringRadius.parseXMLNode( radius, xPath, linkData )
        linkData['conversionTable'] = {'index':int, 'L':int}    # inform table class how to treat columns
        parameters = resonanceParameters.parseXMLNode( element.find( resonanceParameters.moniker ),
                xPath, linkData )
        attrs = getAttrs( element )
        label = attrs.pop("label")
        approximation = attrs.pop("approximation")
        resonanceData = cls( label, approximation, parameters, radius, pops, **attrs )
        del linkData['conversionTable']

        documentation = element.find( documentationModule.Documentation.moniker )
        if( documentation is not None ) : resonanceData.documentation.parseNode( documentation, xPath, linkData )

        xPath.pop()
        return resonanceData

class RMatrix( ancestryModule.ancestry ):
    """
    RMatrix is a general container for resolved resonance parameters.
    It can handle standard Lane & Thomas R-Matrix, but can also use various approximations
    (Single- and Multi-level Breit Wigner plus Reich-Moore).

    Internally, resonances are sorted into spin groups, each with a conserved total angular momentum and parity.
    """

    moniker = 'RMatrix'
    ancestryMembers = ( 'resonanceReactions', 'spinGroups', 'PoPs', 'documentation' )
    ReichMooreToken = "Reich_Moore"
    RMatrixToken = "Full R-Matrix"

    optAttrList = ('approximation', 'boundaryCondition', 'boundaryConditionValue', 'relativisticKinematics',
                   'supportsAngularReconstruction', 'reducedWidthAmplitudes', 'calculateShift',
                   'calculateChannelRadius', 'useForSelfShieldingOnly')

    def __init__(self, label, resonanceReactions_, spinGroups_, PoPs=None, **kwargs):
        """
        Container for R-Matrix resonance parameters

        :param label: unique label corresponding to a style (i.e., 'eval')
        :param resonanceReactions_: resonanceReactions instance
        :param spinGroups_: spinGroups instance
        :param PoPs: optional PoPs database instance
        :param kwargs: see RMatrix.optAttrList for supported options
        """

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.resonanceReactions = resonanceReactions_
        self.spinGroups = spinGroups_
        self.PoPs = PoPs
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

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
    def spinGroups(self):
        return self.__spinGroups

    @spinGroups.setter
    def spinGroups(self, value):
        if not isinstance(value, spinGroups):
            raise TypeError("Must be a spinGroups instance")
        value.setAncestor(self)
        self.__spinGroups = value

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

    def check( self, info ):
        return _resonance_checker(self, info, [self.resonanceReactions] + list(self.spinGroups))

    def convertUnits( self, unitMap ):

        for reac in self.resonanceReactions:
            reac.convertUnits( unitMap )
        for sg in self.spinGroups:
            sg.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + '  '
        xmlString = ['%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        for attr in self.optAttrList:
            if getattr(self,attr):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString[0] += ' %s="%s"' % (attr,attrVal)
        xmlString[0] += '>'

        xmlString += self.__documentation.toXMLList( indent = indent2, **kwargs )

        if self.__PoPs is not None:
            xmlString += self.__PoPs.toXMLList( indent=indent2, **kwargs )
        xmlString += self.resonanceReactions.toXMLList( indent=indent2, **kwargs )
        xmlString += self.spinGroups.toXMLList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )

        pops = element.find( PoPsDatabaseModule.database.moniker )
        if pops is not None:
            pops = PoPsDatabaseModule.database.parseXMLNodeAsClass( pops, xPath, linkData )

        RRs = resonanceReactions()
        RRs.parseXMLNode(element.find(resonanceReactions.moniker), xPath, linkData)
        linkData['conversionTable'] = {'index':int, 'L':int, 'channelSpin':spin}
        SGs = spinGroups()
        SGs.parseXMLNode( element.find(spinGroups.moniker), xPath, linkData )
        attrs = getAttrs(element, required=RMatrix.optAttrList)
        label = attrs.pop("label")

        tmp = RMatrix( label, RRs, SGs, PoPs=pops, **attrs )
        del linkData['conversionTable']

        documentation = element.find( documentationModule.Documentation.moniker )
        if( documentation is not None ) : tmp.documentation.parseNode( documentation, xPath, linkData )

        xPath.pop()
        return tmp

class BoundaryCondition:
    """
    Defines allowed values for the 'boundaryCondition' attribute.
    """
    EliminateShiftFunction = "EliminateShiftFunction"
    NegativeOrbitalMomentum = "NegativeOrbitalMomentum"
    Brune = "Brune"
    Given = "Given"

class channels( suitesModule.suite ):
    """
    Stores a list of channels (used to override global definitions for a single spinGroup)
    """

    moniker = 'channels'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [channel])

class channel( ancestryModule.ancestry ):
    """
    Defines an open channel for a single R-Matrix spin group.
    May be used to override the scattering radius or hard-sphere radius, etc. for this channel.
    """

    moniker = 'channel'
    ancestryMembers = ( 'scatteringRadius', 'hardSphereRadius' )

    def __init__(self, label, resonanceReaction, columnIndex, L=None, channelSpin=None,
                 scatteringRadius=None, hardSphereRadius=None, penetrability=None,
                 shiftFactor=None, phaseShift=None, boundaryConditionValue=None,
                 externalRMatrix=None):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.resonanceReaction = resonanceReaction
        self.columnIndex = columnIndex
        self.L = L
        self.channelSpin = channelSpin
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        if( self.hardSphereRadius is not None ) : self.hardSphereRadius.setAncestor( self )
        self.penetrability = penetrability
        self.shiftFactor = shiftFactor
        self.phaseShift = phaseShift
        self.boundaryConditionValue = boundaryConditionValue
        self.externalRMatrix = externalRMatrix
        if( self.externalRMatrix is not None ) : self.externalRMatrix.setAncestor( self )

    @property
    def scatteringRadius(self):
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            return self.findClassInAncestry(RMatrix).resonanceReactions[self.resonanceReaction].scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a scatteringRadius instance."""

        if( value is not None ) :
            if not isinstance( value, scatteringRadius ) :
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__scatteringRadius = value

    def convertUnits( self, unitMap ):
        for child in ('scatteringRadius', 'hardSphereRadius', 'externalRMatrix'):
            if getattr(self,child) is not None:
                getattr(self,child).convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        attrs = ""
        if self.L is not None: attrs += ' L="%d"' % self.L
        if self.channelSpin is not None: attrs += ' channelSpin="%s"' % self.channelSpin
        if self.boundaryConditionValue is not None:
            attrs += ' boundaryConditionValue="%s"' % self.boundaryConditionValue
        xmlString = ['%s<%s label="%s" resonanceReaction="%s"%s columnIndex="%d">' %
                     (indent, self.moniker, self.label, self.resonanceReaction, attrs, self.columnIndex) ]
        if self.externalRMatrix is not None:
            xmlString += self.externalRMatrix.toXMLList(indent=indent2, **kwargs)
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXMLList(indent=indent2, **kwargs)
        for attr in ('hardSphereRadius', 'penetrability', 'shiftFactor', 'phaseShift'):
            if getattr(self, attr) is not None:
                xmlString += getattr(self, attr).toXMLList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        kwargs = getAttrs(element)
        for child in element:
            childClass = {'externalRMatrix': externalRMatrixModule.externalRMatrix,
                          'scatteringRadius': scatteringRadius,
                          'hardSphereRadius': hardSphereRadius,
                          'penetrability': None,    # FIXME last three not yet implemented
                          'shiftFactor': None,
                          'phaseShift': None}.get( child.tag )
            kwargs[ child.tag ] = childClass.parseXMLNode( child, xPath, linkData )
        tmp = channel( **kwargs )
        xPath.pop()
        return tmp

class spinGroups( suitesModule.suite ):
    """
    Contains a list of spinGroup nodes
    """

    moniker = 'spinGroups'

    def __init__(self):

        suitesModule.suite.__init__(self, [spinGroup])

class spinGroup( ancestryModule.ancestry ):
    """
    Single group with same Jpi (conserved). Each spin group contains an AP (scattering radius),
    along with 1 or more resonance widths.
    """

    moniker = 'spinGroup'
    ancestryMembers = ( 'channels', )       # FIXME, this is not complete.

    def __init__(self, label, spin, parity, channels, resonanceParameters):
        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.spin = spin
        self.parity = parity
        self.channels = channels
        self.resonanceParameters = resonanceParameters

    def __getitem__(self, idx):
        return self.resonanceParameters[idx]

    def __len__(self):
        return len(self.resonanceParameters.table)

    def __lt__(self, other):
        """ for sorting spin groups by Jpi. group J values together """
        return (self.spin, self.parity) < (other.spin, other.parity)

    @property
    def channels(self):
        return self.__channels

    @channels.setter
    def channels(self, value):
        if not isinstance(value, channels):
            raise TypeError("Must be a channels instance")
        value.setAncestor(self)
        self.__channels = value

    @property
    def resonanceParameters(self):
        return self.__resonanceParameters

    @resonanceParameters.setter
    def resonanceParameters(self, value):
        if not isinstance(value, resonanceParameters):
            raise TypeError("Must be a resonanceParameters instance")
        value.setAncestor(self)
        self.__resonanceParameters = value

    def check( self, info ):
        from fudge import warning
        warnings = []
        for thing in []:
            if thing is  None: continue
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.context(thing.moniker,warningList))
        return warnings

    def convertUnits( self, unitMap ):
        for chan in self.channels:
            chan.convertUnits(unitMap)
        self.resonanceParameters.convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        xml = ['%s<%s label="%s" spin="%s" parity="%s">' % (indent, self.moniker, self.label, self.spin, self.parity)]
        xml += self.channels.toXMLList(indent=indent2, **kwargs)
        if self.resonanceParameters.table.columns:    # need not contain any data
            xml.extend( self.resonanceParameters.toXMLList( indent2, **kwargs ) )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        parameters = resonanceParameters.parseXMLNode(element.find(resonanceParameters.moniker), xPath, linkData)
        chs = channels()
        chs.parseXMLNode(element.find(channels.moniker), xPath, linkData)
        SG = spinGroup(channels=chs, resonanceParameters=parameters, **getAttrs(element))
        xPath.pop()

        return SG


def _resonance_checker(self, info, things):
    """
    Common resolved resonance checks

    :param self: a reference to the resonanceFormalismBaseClass (or derived) or RMatrix class
    :param info: the dictionary of checking options, shared by all .check() member functions.
                 For the potential scattering convergence test, the 'potentialScatteringCoverganceCriteria'
                 value is the limit on the ratio of the Lth cross section term to the 0th one.  It is printed
                 as a percent but should be entered as a fraction.
    :param things: list of things to check the obvious way (with a .check() function)
    :return:
    """
    from fudge import warning
    from collections import OrderedDict
    warnings = []

    # Check the member data
    for thing in things:
        warningList = thing.check(info)
        if warningList:
            warnings.append(warning.context(thing.moniker, warningList))

    # Setup for the physics checks
    import fudge.processing.resonances.reconstructResonances as rrReconstructModule
    rrReconstructor = rrReconstructModule.getResonanceReconstructionClass(self)(self)

    # Check for mismatched spins
    warnings += rrReconstructor.setResonanceParametersByChannel(warnOnly=True)  # multipleSScheme='NJOY', 'ENDF' or None

    # Check for missing channels via sum rule of statistical weights
    sumgJ = {}
    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if (c.reaction, c.l) not in sumgJ:
                sumgJ[(c.reaction, c.l)] = 0.0
            sumgJ[(c.reaction, c.l)] += c.gfact
    keys = sorted(sumgJ.keys())
    for rxn, L in keys:
        if abs(abs(sumgJ[(rxn, L)]) - abs(2.0*L+1.0)) > 1e-6:
            warnings.append(warning.badSpinStatisticalWeights(L, sumgJ[(rxn, L)], 2.*L+1, rxn))

    # setup for checking allowed angular momentum:
    def getSList(Ia, Ib):
        """Get possible spins"""
        smin = abs(Ia-Ib)
        smax = Ia+Ib
        nS = int(smax-smin)+1
        return [iS+smin for iS in range(nS)]
    LMax = 0
    allSs = []
    rxnList = []

    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if c.eliminated:
                continue
            LMax = max(c.l, LMax)
            if c.reaction not in rxnList:
                rxnList.append(c.reaction)
            if c.channelClass not in [rrReconstructModule.FISSIONCHANNEL, rrReconstructModule.COMPETITIVECHANNEL]:
                try:
                    spinList = getSList(*rrReconstructor.getParticleSpins(c.reaction))
                    for s in spinList:
                        if s not in allSs:
                            allSs.append(s)
                except:
                    warnings.append(warning.unknownSpinParity(c.reaction))

    # determine the min & max J allowed
    Jmin = min(allSs)
    Jmax = LMax + max(allSs)
    nJ = int(Jmax - Jmin)

    # Check the allowed angular momenta
    for L in range(0, LMax+1):
        for iJ in range(0, nJ+1):
            J = iJ+Jmin
            for rxn in rxnList:
                if 'ission' not in rxn:
                    try:
                        spinList = getSList(*rrReconstructor.getParticleSpins(c.reaction))
                    except:
                        warnings.append(warning.unknownSpinParity(c.reaction))
                        continue

                    for S in spinList:
                        if J not in rrReconstructModule.getAllowedTotalSpins(L, S, useFactor2Trick=False):
                            continue
                        gotIt = False
                        for cc in rrReconstructor.channels:
                            if isinstance(cc, OrderedDict):
                                iterateThroughThis = cc  # for SLBW
                            else:
                                iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
                            for c in iterateThroughThis:
                                if rxn == c.reaction and \
                                        rrReconstructModule.spins_equal(c.l, L) and \
                                        rrReconstructModule.spins_equal(c.J, J) and \
                                        rrReconstructModule.spins_equal(c.s, S):
                                    gotIt = True
                            if gotIt:
                                break
                        if not (gotIt or 'apture' in rxn):
                            theWarning = warning.missingResonanceChannel(L, S, J, rxn)
                            if str(theWarning) not in [str(w) for w in warnings]:
                                warnings.append(theWarning)

    # Check for convergence in L
    import numpy
    almostXS = {}
    for cc in rrReconstructor.channels:
        if isinstance(cc, OrderedDict):
            iterateThroughThis = cc  # for SLBW
        else:
            iterateThroughThis = [cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            egrid = numpy.array([max(rrReconstructor.lowerBound, rrReconstructor.lowerBound + c.Xi),
                                 rrReconstructor.upperBound])
            phis = rrReconstructor.phiByChannel(c, egrid)
            almostXS.setdefault(c.l, []).extend(list(pow(numpy.sin(phis)/rrReconstructor.k(egrid), 2.0)))
    fom = max(almostXS[max(almostXS.keys())])/max(almostXS[min(almostXS.keys())])
    fomTarget = info.get('potentialScatteringCoverganceCriteria', 0.001)
    if fom > fomTarget:
        warnings.append(warning.potentialScatteringNotConverged(c.l, rrReconstructor.upperBound, fom, fomTarget))

    return warnings
