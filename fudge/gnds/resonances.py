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
This module contains the resonances, resolved and unresolved classes.

class hierarchy:

    * resonances contains one or more of scatteringRadius, resolved and unresolved

    * resolved and unresolved contain a list of energy regions. In each region,

        - resolved may have SLBW, MLBW, RM, RML formats
        - unresolved has L-dependant or E-dependant format

      each of these subsections has a list of resonances and attributes
"""

# TODO: move to a 'resonances' module

import fractions
from pqu import PQU
from xData import XYs as XYsModule, constant as constantModule, link as linkModule, table as tableModule, regions as regionsModule
from fudge.gnds import suites as suitesModule
from fudge.gnds import abstractClasses as abstractClassesModule
from fudge.gnds.channelData import Q as QModule

import xData.ancestry as ancestryModule

__metaclass__ = type


class resonances( ancestryModule.ancestry ) :
    """ 
    This is the top-level class for storing resonance parameters.
    For light targets it may contain only a scattering radius. For heavier targets it typically
    contains a resolved and/or unresolved section.
    
    resonances also has a boolean flag 'reconstructCrossSection'. If False, either the cross section
    has already been reconstructed, or the parameters are given for information only and no reconstruction
    should be performed.
    """

    moniker = 'resonances'
    children = ('scatteringRadius', 'resolved', 'unresolved')

    def __init__(self, scatteringRadius_, resolved_=None, unresolved_=None):

        ancestryModule.ancestry.__init__( self )

        self.scatteringRadius = scatteringRadius_
        self.resolved = resolved_
        self.unresolved = unresolved_

        for child in (self.scatteringRadius, self.resolved, self.unresolved):
            if child is not None: child.setAncestor( self )
        
    def  __str__( self ) :
        """ string representation """
        return( self.toString( simpleString = False ) )

    def convertUnits( self, unitMap ):
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        for child in self.children:
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits( unitMap )

    def check( self, info ):
        from fudge.gnds import warning
        warnings = []
        for child in self.children:
            section = getattr(self,child)
            if section is not None:
                warningList = section.check(info)
                if warningList:
                    warnings.append( warning.context( section.moniker, warningList ) )
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
    
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for child in self.children:
            section = getattr(self, child)
            if section is not None:
                xmlString += section.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def domain( self, unitTo = None, asPQU = False ):
        """ Return resonance region domain as a tuple of floats: (lowest edge, highest edge).

        options:
          unitTo: convert output to specified unit (given as a string).
          asPQU = True: return a tuple of PhysicalQuantityWithUncertainty instances instead of floats.
        """

        bounds = [ (self.scatteringRadius.domainMin, self.scatteringRadius.domainMax) ]
        if self.resolved:
            if self.resolved.multipleRegions:
                bounds += [(reg.domainMin, reg.domainMax) for reg in self.resolved.regions]
            else: bounds.append( (self.resolved.domainMin, self.resolved.domainMax) )
        if self.unresolved:
            bounds.append( (self.unresolved.domainMin, self.unresolved.domainMax) )

        for idx in range(len(bounds)-1):
            assert bounds[idx][1] == bounds[idx+1][0], "Resonance region boundaries don't match!"

        if( asPQU ):
            return (bounds[0][0], bounds[-1][1])
        elif unitTo:
            return (bounds[0][0].getValue(unitTo), bounds[-1][1].getValue(unitTo))
        else:
            return (bounds[0][0].value, bounds[-1][1].value)

    @property
    def reconstructCrossSection( self ):
        if self.resolved and not self.resolved.evaluated.useForSelfShieldingOnly: return True
        if self.unresolved and not self.unresolved.evaluated.useForSelfShieldingOnly: return True
        return False

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )

        scatRadius, RRR, URR = None,None,None
        for child in element:
            if child.tag==scatteringRadius.moniker:
                scatRadius = scatteringRadius.parseXMLNode( child, xPath, linkData )
            elif child.tag==resolved.moniker:
                RRR = resolved.parseXMLNode( child, xPath, linkData )
            elif child.tag==unresolved.moniker:
                URR = unresolved.parseXMLNode( child, xPath, linkData )
            else:
                raise Exception("unknown element '%s' encountered in resonances!" % child.tag)

        res = cls( scatteringRadius_ = scatRadius, resolved_=RRR, unresolved_=URR )
        xPath.pop()
        return res
    
    def toString( self, simpleString = False ) :
        """Returns a string representation of self. If simpleString is True, 
        the string contains only an overview without listing resonances"""
        s = 'Resonances:\n'
        if self.scatteringRadius:
            s += self.scatteringRadius.toString( simpleString = simpleString )
        if self.resolved:
            s += self.resolved.toString( simpleString = simpleString )
        if self.unresolved:
            s += self.unresolved.toString( simpleString = simpleString )
        return( s )


class scatteringRadius( ancestryModule.ancestry ):
    """
    The scatteringRadius determines scattering length. May be constant or energy-dependent.
    Each resonances section contains a scatteringRadius. The resolved / unresolved regions may also
    provide their own copy that overrides the top-level definition.
    The RMatrix class also provides a way of defining channel-specific radii (deprecated but necessary to handle
    ENDF evaluations)
    """

    moniker = 'scatteringRadius'

    def __init__( self, form=None ) :

        ancestryModule.ancestry.__init__( self )
        form.setAncestor( self )
        self.form = form

    def __eq__(self, other):
        if not isinstance(other, scatteringRadius): return False
        return( self.form==other.form )

    def __str__(self): return self.form.moniker

    def __nonzero__(self): return bool(self.form)

    def toString( self, simpleString = False ) :
        """Returns a string representation of self. If simpleString is True,
        the string contains only an overview without listing resonances"""
        if simpleString: return str(self)
        return str(self)

    def check( self, info ):
        """
        Checks that the scattering radius is within a factor of 3 of the expected scattering radius
        of AP = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08

        This default AP is roughly the nuclear radius, so getting within a factor of 3 of the default
        shouldn't be a problem.

        A similar test is included in PSYCHE, but PSYCHE's cannot handle LRF=7
        :param info:
        :return:
        """
        from fudge.gnds import warning
        warnings = []
        target = info['reactionSuite'].PoPs[ info['reactionSuite'].target ]
        if target.id in info['reactionSuite'].PoPs.aliases:
            target = info['reactionSuite'].PoPs[ target.pid ]
        expectedAP = 10.0*( 0.123 * target.getMass('amu')**(1./3.) + 0.08 ) # expected radius in fm
        factor=3.0
        if self.isEnergyDependent():
            egrid=self.form.domainGrid
            APs = self.getValueAs( 'fm', energy_grid=egrid )
            for iE,AP in enumerate(APs):
                if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                    warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP, E=egrid[iE])
        else:
            AP = self.form.constant
            if self.form.rangeUnit != 'fm':
                AP *= PQU.PQU(1, self.form.rangeUnit).getValueAs('fm')
            if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP)
        return warnings

    def convertUnits( self, unitMap ):

        self.form.convertUnits( unitMap )

    def isEnergyDependent(self):
        return isinstance( self.form, XYsModule.XYs1d )

    def getValueAs( self, unit, energy_grid=None ):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise NameError("Missing: energy_grid to evaluate E-dependent scattering radius")
            energy_unit = self.form.axes[-1].unit
            xScale = self.form.domainUnitConversionFactor( energy_unit )
            yScale = self.form.rangeUnitConversionFactor( unit )
            return [ yScale * self.form.evaluate( xScale * e ) for e in energy_grid ]
        else:
            oldUnit = self.form.rangeUnit
            factor = PQU.PQU(1, oldUnit).getValueAs( unit )
            return self.form.constant * factor

    def toXMLList( self, indent = '', **kwargs ):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.form.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )
        form = {
            constantModule.constant1d.moniker: constantModule.constant1d,
            XYsModule.XYs1d.moniker: XYsModule.XYs1d,
        }[ element[0].tag ].parseXMLNode( element[0], xPath, linkData )
        SR = cls( form )
        xPath.pop()
        return SR


class hardSphereRadius( scatteringRadius):

    moniker = 'hardSphereRadius'


class spin( fractions.Fraction ):
    """
    Store spins for a collection of resonances. Check denominator (must be integer or half-integer)
    """

    def __init__( self, *args ):
        fractions.Fraction.__init__( self, *args )
        if self.denominator not in (1,2):
            raise ValueError("Illegal spin '%s': must be integer or half-integer" % self)


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


class resonanceParameters( ancestryModule.ancestry ):   # FIXME: still needed? It's an extra level of nesting...
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


# resolved/unresolved regions:
class resolved( abstractClassesModule.component ):
    """ class for resolved resonances """

    moniker = 'resolved'

    def __init__(self, domainMin, domainMax, domainUnit):

        abstractClassesModule.component.__init__( self,
            allowedClasses=(energyIntervals, BreitWigner, RMatrix) )
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

    def toString(self, simpleString = False):
        if isinstance( self.evaluated, energyIntervals ):
            return ("Resolved region with DEPRECATED multiple regions\n")
        else:
            return ("Resolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        import warning
        warnings = []
        if isinstance(self.evaluated, energyIntervals):
            warnings.append( warning.RRmultipleRegions() )
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.context(self.evaluated.moniker,warningList))
        return warnings

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQU.PQU(1, self.domainUnit).getValueAs( newUnit )
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit
        for form in self:
            form.convertUnits( unitMap )
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s domainMin="%r" domainMax="%r" domainUnit="%s">' % ( indent, self.moniker,
                self.domainMin, self.domainMax, self.domainUnit) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
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


class unresolved( abstractClassesModule.component ):

    moniker = 'unresolved'

    def __init__(self, domainMin, domainMax, domainUnit):
        abstractClassesModule.component.__init__( self, allowedClasses=(unresolvedTabulatedWidths, energyIntervals) )
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit
    
    def toString( self, simpleString = False ):
        return ("Unresolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        from fudge.gnds import warning
        warnings = []
        for L in self.evaluated.Ls:
            for J in L.Js:
                elist = J.levelSpacing.data.convertAxisToUnit(1, self.domainUnit).domainGrid
                if elist[0] > self.domainMin or elist[-1] < self.domainMax:
                    warnings.append( warning.URRdomainMismatch( L.L, J.J, J ) )
                missingPoints = [i1 for i1 in range(1,len(elist)) if elist[i1] > 3*elist[i1-1]]
                for idx in missingPoints:
                    warnings.append( warning.URRinsufficientEnergyGrid( L.L, J.J,
                        PQU.PQU(elist[idx-1],self.domainUnit), PQU.PQU(elist[idx],self.domainUnit), J ) )
        return warnings

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQU.PQU(1, self.domainUnit).getValueAs(newUnit)
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit
        for form in self:
            form.convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s domainMin="%s" domainMax="%s" domainUnit="%s">' %
                ( indent, self.moniker, self.domainMin, self.domainMax, self.domainUnit ) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )

        URR = cls( **getAttrs(element) )
        for child in element:
            formClass = {
                    unresolvedTabulatedWidths.moniker: unresolvedTabulatedWidths,
                    energyIntervals.moniker: energyIntervals,
                    }.get( child.tag )
            if formClass is None: raise Exception("unknown unresolved resonance form '%s'!" % child.tag)
            URR.add(formClass.parseXMLNode( child, xPath, linkData ))

        xPath.pop()
        return URR


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
        import warning
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
        import warning
        warnings = []
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.context( self.evaluated.moniker, warningList ))
        return warnings

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQU.PQU(1, self.domainUnit).getValueAs(newUnit)
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit
        self.evaluated.convertUnits(unitMap)

    def toString(self, simpleString = False):
        return ("%s resonances, %s to %s. Contains %i resonances" % 
                (self.evaluated, self.domainMin, self.domainMax, len(self.evaluated) ) )
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = ['%s<%s index="%s" domainMin="%r" domainMax="%r" domainUnit="%s">'
                % ( indent, self.moniker, self.index, self.domainMin, self.domainMax, self.domainUnit ) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @property
    def BreitWigner(self):
        if isinstance(self.evaluated, BreitWigner): return self.evaluated
        return None

    @property
    def RMatrix(self):
        if isinstance(self.evaluated, RMatrix): return self.evaluated
        return None

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):
        xPath.append('%s[@label="%s"]' % (cls.moniker, element.get('label')))

        formClass = {
            BreitWigner.moniker: BreitWigner,
            RMatrix.moniker: RMatrix,
            unresolvedTabulatedWidths.moniker: unresolvedTabulatedWidths,
        }.get( element[0].tag )
        if formClass is None: raise Exception("unknown unresolved resonance form '%s'!" % element[0].tag)
        data = formClass.parseXMLNode(element[0], xPath, linkData)
        EI = cls( data=data, **getAttrs(element) )

        xPath.pop()
        return EI


def _resonance_checker(self,info,things):
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
    from fudge.gnds import warning
    from collections import OrderedDict
    warnings=[]

    # Check the member data
    for thing in things:
        warningList = thing.check(info)
        if warningList:
            warnings.append(warning.context(thing.moniker,warningList))

    # Setup for the physics checks
    import fudge.processing.resonances.reconstructResonances as rrReconstructModule
    rrReconstructor = rrReconstructModule.getResonanceReconstructionClass(self)(
            info['reactionSuite'], sectionIndex = info.get('energyIntervalIndex') )

    # Check for mismatched spins
    warnings+=rrReconstructor.setResonanceParametersByChannel(warnOnly=True) #multipleSScheme='NJOY', 'ENDF' or None

    # Check for missing channels via sum rule of statistical weights
    sumgJ={}
    for cc in rrReconstructor.channels:
        if isinstance(cc,OrderedDict): iterateThroughThis=cc # for SLBW
        else: iterateThroughThis=[cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if (c.reaction,c.l) not in sumgJ: sumgJ[(c.reaction,c.l)]=0.0
            sumgJ[(c.reaction,c.l)]+=c.gfact
    keys = sumgJ.keys()
    keys.sort()
    for rxn,L in keys:
        if abs(abs(sumgJ[(rxn,L)]) - abs(2.0*L+1.0)) > 1e-6:
            warnings.append(warning.badSpinStatisticalWeights(L,sumgJ[(rxn,L)],2.*L+1,rxn))

    # setup for checking allowed angular momentum:
    def getSList(Ia,Ib):
        '''Get possible spins'''
        smin=abs(Ia-Ib)
        smax=Ia+Ib
        nS = int(smax-smin)+1
        return [iS+smin for iS in range(nS)]
    LMax = 0
    allSs = []
    rxnList = []

    for cc in rrReconstructor.channels:
        if isinstance(cc,OrderedDict): iterateThroughThis=cc # for SLBW
        else: iterateThroughThis=[cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            if c.eliminated: continue
            LMax = max(c.l,LMax)
            if c.reaction not in rxnList: rxnList.append(c.reaction)
            if c.channelClass not in [rrReconstructModule.FISSIONCHANNEL,rrReconstructModule.COMPETATIVECHANNEL]:
                try:
                    spinList = getSList(*rrReconstructor.getParticleSpins(c.reaction))
                    for s in spinList:
                        if s not in allSs: allSs.append(s)
                except:
                    warnings.append(warning.unknownSpinParity(c.reaction))
    # determine the min & max J allowed
    Jmin = min(allSs)
    Jmax = LMax + max(allSs)
    nJ = int(Jmax - Jmin)

    # Check the allowed angular momenta
    for L in range(0,LMax+1):
        for iJ in range(0,nJ+1):
            J=iJ+Jmin
            for rxn in rxnList:
                if 'ission' not in rxn:
                    try:
                        spinList = getSList(*rrReconstructor.getParticleSpins(c.reaction))
                    except:
                        warnings.append(warning.unknownSpinParity(c.reaction))
                        continue

                    for S in spinList:
                        if not J in rrReconstructModule.getAllowedTotalSpins(L,S,useFactor2Trick=False): continue
                        for cc in rrReconstructor.channels:
                            if isinstance(cc, OrderedDict): iterateThroughThis = cc # for SLBW
                            else: iterateThroughThis = [cc] # so everyone else can iterate like SLBW
                            for c in iterateThroughThis:
                                gotIt = False
                                if rxn==c.reaction and \
                                        rrReconstructModule.spins_equal(c.l,L) and \
                                        rrReconstructModule.spins_equal(c.J,J) and \
                                        rrReconstructModule.spins_equal(c.s,S):
                                    gotIt=True
                                    break
                            if not gotIt and 'apture' not in rxn:
                                theWarning=warning.missingResonanceChannel(L,S,J,rxn)
                                if str(theWarning) not in map(str,warnings):
                                    warnings.append(theWarning)

    # Check for convergence in L
    import numpy
    almostXS={}
    for cc in rrReconstructor.channels:
        if isinstance(cc,OrderedDict): iterateThroughThis=cc # for SLBW
        else: iterateThroughThis=[cc]  # so everyone else can iterate like SLBW
        for c in iterateThroughThis:
            egrid = numpy.array([max(rrReconstructor.lowerBound, rrReconstructor.lowerBound + c.Xi), rrReconstructor.upperBound])
            phis=rrReconstructor.phiByChannel(c,egrid)
            almostXS.setdefault(c.l,[]).extend( list(pow(numpy.sin(phis)/rrReconstructor.k(egrid),2.0)) )
    fom=max(almostXS[max(almostXS.keys())])/max(almostXS[min(almostXS.keys())])
    fomTarget=info.get('potentialScatteringCoverganceCriteria',0.001)
    if fom > fomTarget:
        warnings.append(warning.potentialScatteringNotConverged(c.l, rrReconstructor.upperBound, fom, fomTarget))

    return warnings


################################################################
# now we come to specific implementations for resolved region: #
################################################################
class BreitWigner( ancestryModule.ancestry ) :

    moniker = 'BreitWigner'
    singleLevel = 'singleLevel'
    multiLevel = 'multiLevel'
    optAttrList = ( 'calculateChannelRadius', 'computeAngularDistribution', 'LvaluesNeededForConvergence',
                    'useForSelfShieldingOnly')

    def __init__(self, label, approximation, resonanceParameters=None, scatteringRadius=None, **kwargs):
        """LdependentScatteringRadii is a list of dictionaries, one for each L with APL != AP
        only used in ENDF for Reich_Moore form
        
        resonanceParameters should be a gnds.core.math.table.table instance."""

        index = 0
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
        self.__scatteringRadius = scatteringRadius
        if self.__scatteringRadius: self.__scatteringRadius.setAncestor( self )
        ancestryModule.ancestry.__init__( self )

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
        """ Can be set to None or to a scatteringRadius instance. """
        self.__scatteringRadius = value
        if value is not None:
            if not isinstance(value, scatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            self.__scatteringRadius.setAncestor(self)

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
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXMLList( indent2, **kwargs )
        if self.resonanceParameters:
            xmlString.extend( self.resonanceParameters.toXMLList( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        radius = element.find( scatteringRadius.moniker )
        if radius is not None:
            radius = scatteringRadius.parseXMLNode( radius, xPath, linkData )
        linkData['conversionTable'] = {'index':int, 'L':int}    # inform table class how to treat columns
        parameters = resonanceParameters.parseXMLNode( element.find( resonanceParameters.moniker ),
                xPath, linkData )
        attrs = getAttrs( element )
        label = attrs.pop("label")
        approximation = attrs.pop("approximation")
        resonanceData = cls( label, approximation, parameters, radius, **attrs )
        del linkData['conversionTable']
        xPath.pop()
        return resonanceData
    
#############################################
# R-matrix (LRF=7 in ENDF) is more complex: #
#############################################
class RMatrix( ancestryModule.ancestry ):
    """
    RMatrix is a general container for resolved resonance parameters.
    It can handle standard Lane & Thomas R-Matrix, but can also use various approximations
    (Single- and Multi-level Breit Wigner plus Reich-Moore).
    
    Internally, resonances are sorted into spin groups, each with a conserved total angular momentum and parity.
    """

    moniker = 'RMatrix'

    optAttrList = ('approximation','boundaryCondition','relativisticKinematics','supportsAngularReconstruction',
            'reducedWidthAmplitudes','calculatePenetrability','calculateShift','calculateChannelRadius',
            'useForSelfShieldingOnly')

    def __init__(self, label, resonanceReactions_, spinGroups_, **kwargs):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.resonanceReactions = resonanceReactions_
        self.spinGroups = spinGroups_
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

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

        xmlString += self.resonanceReactions.toXMLList( indent=indent2, **kwargs )
        xmlString += self.spinGroups.toXMLList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        RRs = resonanceReactions()
        RRs.parseXMLNode( element.find(resonanceReactions.moniker), xPath, linkData )
        linkData['conversionTable'] = {'index':int, 'L':int, 'channelSpin':spin}
        SGs = spinGroups()
        SGs.parseXMLNode( element.find(spinGroups.moniker), xPath, linkData )
        attrs = getAttrs(element, required=RMatrix.optAttrList)
        label = attrs.pop("label")
        tmp = RMatrix( label, RRs, SGs, **attrs )
        del linkData['conversionTable']
        xPath.pop()
        return tmp

class resonanceReactions( suitesModule.suite ):

    moniker = 'resonanceReactions'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [resonanceReaction])

    def check( self, info ):
        from fudge.gnds import warning
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
                  computeShiftFactor=False, eliminated=False):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.reactionLink = reactionLink
        self.__ejectile = ejectile
        self.__residual = None
        self.Q = Q
        self.__scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
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
            return self.findClassInAncestry(resonances).scatteringRadius

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """ Can be set to None or to a scatteringRadius instance. """
        self.__scatteringRadius = value
        if value is not None:
            if not isinstance(value, scatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            self.__scatteringRadius.setAncestor(self)

    def check( self, info ):
        from fudge.gnds import warning
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

        return self.reactionLink.link.outputChannel.isFission()

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        attrstring = ''
        if self.ejectile is not None: attrstring += ' ejectile="%s"' % self.ejectile
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
                computeShiftFactor=getBool( element.get('computeShiftFactor','false') ),
                eliminated=getBool( element.get('eliminated','false') ) )
        xPath.pop()
        return tmp

class channels( suitesModule.suite ):

    moniker = 'channels'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [channel])

class channel( ancestryModule.ancestry ):
    """
    Use if necessary to override the true scattering radius, effective radius, etc. for a single
    spin-group / channel combination.
    """

    moniker = 'channel'

    def __init__(self, label, resonanceReaction, columnIndex, L=None, channelSpin=None,
                 scatteringRadius=None, hardSphereRadius=None, penetrability=None,
                 shiftFactor=None, phaseShift=None, boundaryConditionOverride=None):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.resonanceReaction = resonanceReaction
        self.columnIndex = columnIndex
        self.L = L
        self.channelSpin = channelSpin
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.penetrability = penetrability
        self.shiftFactor = shiftFactor
        self.phaseShift = phaseShift
        self.boundaryConditionOverride = boundaryConditionOverride

    def convertUnits( self, unitMap ):
        for child in ('scatteringRadius','hardSphereRadius'):
            if getattr(self,child) is not None:
                getattr(self,child).convertUnits(unitMap)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        attrs = ""
        if self.L is not None: attrs += ' L="%d"' % self.L
        if self.channelSpin is not None: attrs += ' channelSpin="%s"' % self.channelSpin
        if self.boundaryConditionOverride is not None:
            attrs += ' boundaryConditionOverride="%s"' % self.boundaryConditionOverride
        xmlString = ['%s<%s label="%s" resonanceReaction="%s"%s columnIndex="%d">' %
                     (indent, self.moniker, self.label, self.resonanceReaction, attrs, self.columnIndex) ]
        for attr in ('scatteringRadius', 'hardSphereRadius', 'penetrability', 'shiftFactor', 'phaseShift'):
            if getattr(self, attr) is not None:
                xmlString += getattr(self, attr).toXMLList( indent=indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        kwargs = getAttrs(element)
        for child in element:
            childClass = {'scatteringRadius': scatteringRadius,
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
        # for sorting spin groups by Jpi. group J values together
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
        from fudge.gnds import warning
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


#############################################
#    end of R-Matrix (LRF=7) classes        #
#############################################


####################################################
# remaining classes are for unresolved resonances: #
####################################################
class unresolvedTabulatedWidths(ancestryModule.ancestry):

    moniker = 'tabulatedWidths'

    def __init__(self, label, approximation, resonanceReactions_, Ls_, scatteringRadius=None,
                 useForSelfShieldingOnly=False):
        """
        L_sections is a list of URR_Lsections, which in turn contains a list of URR_Jsections.
        The average widths are given in the J_sections.
        """
        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.approximation = approximation
        self.resonanceReactions = resonanceReactions_
        self.scatteringRadius = scatteringRadius
        self.Ls = Ls_
        self.useForSelfShieldingOnly = useForSelfShieldingOnly

    def __len__(self):
        return len(self.Ls)

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
        if not isinstance(value, URR_Lsections):
            raise TypeError("Must be a URR_Lsections instance")
        value.setAncestor(self)
        self.__Ls = value

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
        Ls_ = URR_Lsections()
        Ls_.parseXMLNode( element.find( Ls_.moniker ), xPath, linkData )

        result = cls( element.get('label'), element.get('approximation'),
            resonanceReactions_, Ls_=Ls_, scatteringRadius=radius,
            useForSelfShieldingOnly=element.get('useForSelfShieldingOnly') == 'true'
        )
        del linkData['conversionTable']
        xPath.pop()
        return result


class URR_interpolationTable(ancestryModule.ancestry):
    """
    Base class for unresolved level spacing and widths
    Data are stored as XYs1d or regions1d.
    """
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


class URR_Lsections( suitesModule.suite ):
    """
    tabulated widths contains a list of 'L' sections, each of which contains a list of 'J' sections
    """

    moniker = "Ls"

    def __init__( self ):
        suitesModule.suite.__init__(self, (URR_Lsection,))


class URR_Lsection(ancestryModule.ancestry):
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
        if not isinstance(value, (URR_Jsections)):
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
        Js = URR_Jsections()
        Js.parseXMLNode( element.find(Js.moniker), xPath, linkData )
        result = cls( element.get('label'), int(element.get('value')), Js )
        xPath.pop()
        return result


class URR_Jsections( suitesModule.suite ):

    moniker = "Js"

    def __init__( self ):
        suitesModule.suite.__init__(self, (URR_Jsection,))


class URR_Jsection(ancestryModule.ancestry):
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
    def widths( self, widths ):
        if not isinstance(widths, URR_widths):
            raise TypeError("Expected URR_widths instance, got %s instead" % type(widths))
        widths.setAncestor(self)
        self.__widths = widths

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
        widths_ = URR_widths()
        widths_.parseXMLNode( element.find(URR_widths.moniker), xPath, linkData )
        Jsec = URR_Jsection( element.get('label'), fractions.Fraction(element.get('value')), levelSpacing_, widths_ )
        xPath.pop()
        return Jsec


class levelSpacing(URR_interpolationTable):
    """
    Contains the average level spacing, stored in XYs1d or regions1d.
    """

    moniker = "levelSpacing"


class URR_widths( suitesModule.suite ):
    """
    Contains all average channel widths for an L/J combination.
    """

    moniker = "widths"

    def __init__(self):
        suitesModule.suite.__init__(self, allowedClasses=(URR_width,))


class URR_width(URR_interpolationTable):
    """
    Stores the average width and number of degrees of freedom for a single channel

    Degrees of freedom are typically 0, 1 or 2 depending on how many channel spins contribute, but it may be
    non-integer for a reaction that has a threshold somewhere inside the unresolved region.

    Average width is stored as XYs1d or regions1d.
    """

    moniker = "width"

    def __init__(self, resonanceReaction, data, degreesOfFreedom = 0):
        URR_interpolationTable.__init__(self, data)
        self.resonanceReaction = resonanceReaction
        self.degreesOfFreedom = degreesOfFreedom

    @property
    def label(self): return self.resonanceReaction

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attrs = ''
        if self.degreesOfFreedom > 0:
            attrs += ' degreesOfFreedom="%g"' % self.degreesOfFreedom
        xml = ['%s<%s resonanceReaction="%s"%s>' % (indent, self.moniker, self.resonanceReaction, attrs)]
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
        result = cls(element.get('resonanceReaction'), data, floatOrint(element.get('degreesOfFreedom',0)))
        xPath.pop()
        return result


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
            'reconstructCrossSection':getBool, 'supportsAngularReconstruction':getBool, 'multipleRegions': getBool,
            'calculateChannelRadius':getBool, 'computeAngularDistribution':getBool, 'useForSelfShieldingOnly':getBool,
            'calculateShift':getBool,'calculatePenetrability':getBool, 'channelIndex':int, 'columnIndex':int,
            'LvaluesNeededForConvergence':int, 'ENDF_MT':int, 'index':int, 'L':int,
            'spin':spin, 'parity':parity, 'boundaryConditionOverride':float,
            }
    attrs = dict( element.items() )
    for key in attrs.keys():
        if key in exclude: attrs.pop(key)
        elif key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
    for val in required:
        if val not in attrs: attrs[val] = False
    return attrs
