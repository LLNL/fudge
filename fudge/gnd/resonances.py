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
from xData import XYs as XYsModule, constant as constantModule, link as linkModule, table as tableModule
from fudge.gnd import suites as suitesModule
from fudge.gnd import abstractClasses as abstractClassesModule

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

    def __init__(self, scatteringRadius=None, resolved=None, unresolved=None):

        ancestryModule.ancestry.__init__( self )

        if (scatteringRadius is not None) and (resolved is not None or unresolved is not None):
            raise TypeError, ("Resonances should contain scattering radius OR resonance parameters, not both")
        self.scatteringRadius = scatteringRadius
        self.resolved = resolved
        self.unresolved = unresolved

        for child in (self.scatteringRadius, self.resolved, self.unresolved):
            if child is not None: child.setAncestor( self )
        
    def  __str__( self ) :
        """ string representation """
        return( self.toString( simpleString = False ) )

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for section in (self.scatteringRadius, self.resolved, self.unresolved):
            if section is not None:
                warningList = section.check(info)
                if warningList:
                    warnings.append( warning.context( section.moniker, warningList ) )
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
    
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        if self.scatteringRadius is not None:
            xmlString += self.scatteringRadius.toXMLList( indent2, **kwargs )
        if self.resolved is not None:
            xmlString += self.resolved.toXMLList( indent2, **kwargs )
        if self.unresolved is not None:
            xmlString += self.unresolved.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def domain( self, unitTo = None, asPQU = False ):
        """ Return resonance region domain as a tuple of floats: (lowest edge, highest edge).

        options:
          unitTo: convert output to specified unit (given as a string).
          asPQU = True: return a tuple of PhysicalQuantityWithUncertainty instances instead of floats.
        """

        bounds = []
        if self.scatteringRadius:
            bounds.append( (self.scatteringRadius.domainMin, self.scatteringRadius.domainMax) )
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
        if self.resolved and self.resolved.reconstructCrossSection: return True
        if self.unresolved and self.unresolved.reconstructCrossSection: return True
        return False

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
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

        res = resonances( scatteringRadius = scatRadius, resolved=RRR, unresolved=URR )
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


# Q-value in resonance region may not be consistent with other parts of the evaluation
class Q( ancestryModule.ancestry ):

    moniker = "Q"

    def __init__( self, value, unit ):
        ancestryModule.ancestry.__init__(self)
        self.value = value
        self.unit = unit

    def getValueAs(self, unit):
        factor = PQU.PQU('1. ' + self.unit).getValueAs(unit)
        return self.value * factor

    def toXMLList( self, indent='', **kwargs ):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml.append( '%s<constant value="%r" unit="%s"/>' % (indent2, self.value, self.unit) )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )

        tmp = Q( float(element.getchildren()[0].get('value')), element.getchildren()[0].get('unit') )
        xPath.pop()
        return tmp


# scatteringRadius can be constant, energy-dependent or L-dependent:
class scatteringRadius( ancestryModule.ancestry ):

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
        '''
        Checks that the scattering radius is within a factor of 3 of the expected scattering radius
        of AP = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08

        This default AP is roughly the nuclear radius, so getting within a factor of 3 of the default
        shouldn't be a problem.

        A similar test is included in PSYCHE, but PSYCHE's cannot handle LRF=7
        :param info:
        :return:
        '''
        from fudge.gnd import warning
        warnings = []
        target = info['reactionSuite'].PoPs[ info['reactionSuite'].target ]
        expectedAP = 10.0*( 0.123 * target.getMass('amu')**(1./3.) + 0.08 ) # expected radius in fm
        factor=3.0
        if self.isEnergyDependent():
            egrid=self.form.domainGrid
# BRB6 hardwired
            APs = self.getValueAs( 'fm', energy_grid=egrid )
            for iE,AP in enumerate(APs):
                if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                    warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP, E=egrid[iE])
        elif isinstance(self.form, constantScatteringRadius):    # FIXME goes away once these are all converted to constant1d
            AP = self.form.value
            if self.form.unit != 'fm':
                AP *= PQU.PQU('1. '+self.form.unit).getValueAs('fm')
            if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP)
        else:
            AP = self.form.constant
            if self.form.rangeUnit != 'fm':
                AP *= PQU.PQU(1, self.form.rangeUnit).getValueAs('fm')
            if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP)
        return warnings

    def isEnergyDependent(self):
        return isinstance( self.form, XYsModule.XYs1d )

    def getValueAs( self, unit, energy_grid=None, L=None ):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise NameError("Missing: energy_grid to evaluate E-dependent scattering radius")
            energy_unit = self.form.axes[-1].unit
            xScale = self.form.domainUnitConversionFactor( energy_unit )
            yScale = self.form.rangeUnitConversionFactor( unit )
            return [ yScale * self.form.evaluate( xScale * e ) for e in energy_grid ]
        elif isinstance(self.form, constantScatteringRadius):   # FIXME goes away once these are converted to constant1d
            oldUnit  = self.form.unit
            factor = PQU.PQU(1,oldUnit).getValueAs( unit )
            return self.form.value * factor
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
            constantScatteringRadius.moniker: constantScatteringRadius,
            XYsModule.XYs1d.moniker: XYsModule.XYs1d,
        }[ element[0].tag ].parseXMLNode( element[0], xPath, linkData )
        SR = cls( form )
        xPath.pop()
        return SR


class hardSphereRadius( scatteringRadius):

    moniker = 'hardSphereRadius'


class baseScatteringRadiusForm( abstractClassesModule.form ) :

    pass

class constantScatteringRadius( baseScatteringRadiusForm ):

    moniker = 'constant'

    def __init__( self, value, unit, domainMin = None, domainMax = None, domainUnit = None, label = None ) :

        baseScatteringRadiusForm.__init__(self)
        if isinstance( value, str ):
            value = PQU.PQU( value )
        self.value = value
        self.unit = unit
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

        if( label is not None ) :
            if not isinstance( label, str ):
                raise TypeError( 'label must be a string' )
        self.__label = label

    @property
    def label( self ) :

        return( self.__label )

    def toXMLList( self, indent = '', **kwargs ):

        attributeStr, boundsString = '', ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if self.domainMin is not None:
            boundsString = ' domainMin="%s" domainMax="%s" domainUnit="%s"' % (
                self.domainMin, self.domainMax, self.domainUnit )
        return [ '%s<%s%s value="%s" unit="%s" %s/>' %
                 ( indent, self.moniker, attributeStr, self.value, self.unit, boundsString )]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        CSR = constantScatteringRadius( **getAttrs(element) )
        xPath.pop()
        return CSR


class spin( fractions.Fraction ):
    """
    Store spins for a collection of resonances. Check denominator (must be integer or half-integer)
    """

    def __init__( self, *args ):
        fractions.Fraction.__init__( self, *args )
        if self.denominator not in (1,2):
            raise ValueError("Illegal spin '%s': must be integer or half-integer" % self)


class parity:
    """Store parity for a collection of resonances. Allowed values for the parity are +1 or -1."""

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
class resolved( ancestryModule.ancestry ):
    """ class for resolved resonances """

    moniker = 'resolved'

    def __init__(self, formalism, domainMin, domainMax, domainUnit, reconstructCrossSection=True,
                 reconstructAngularDistributions=False, multipleRegions=False):

        ancestryModule.ancestry.__init__( self )
        if not multipleRegions:
            setattr(self, formalism.moniker, formalism)
        self.evaluated = formalism
        if self.evaluated is not None: self.evaluated.setAncestor( self )
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit
        self.reconstructCrossSection = reconstructCrossSection
        self.reconstructAngularDistributions = reconstructAngularDistributions
        self.multipleRegions = multipleRegions
        self.regions = []  # contains multiple energyInterval (deprecated)

    def toString(self, simpleString = False):
        if self.regions:
            return ("Resolved region with DEPRECATED multiple regions\n")
        else:
            return ("Resolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        import warning
        warnings = []
        if self.multipleRegions:
            warnings.append( warning.RRmultipleRegions() )
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.context(self.evaluated.moniker,warningList))
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ' reconstructCrossSection="%s"' % str(self.reconstructCrossSection).lower()
        if self.reconstructAngularDistributions:
            attrs += ' reconstructAngularDistributions="true"'
        xmlString = [ '%s<%s' % ( indent, self.moniker ) ]
        if self.multipleRegions:
            xmlString[0] += ' multipleRegions="true" domainMin="%r" domainMax="%r" domainUnit="%s"%s>' % (
                self.domainMin, self.domainMax, self.domainUnit, attrs )
            for region in self.regions:
                xmlString += region.toXMLList( indent2, **kwargs )
        else:
            xmlString[0] += ' domainMin="%r" domainMax="%r" domainUnit="%s" formalism="%s"%s>' % (
                    self.domainMin, self.domainMax, self.domainUnit, self.evaluated.moniker, attrs )
            xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        def readResolved( child ):
            formClass = {SLBW.moniker: SLBW,
                    MLBW.moniker: MLBW,
                    RM.moniker: RM,
                    RMatrix.moniker: RMatrix
                    }.get( child.tag )
            if formClass is None: raise Exception("unknown resolved resonance form '%s'!" % child.tag)
            return formClass.parseXMLNode( child, xPath, linkData )

        regions = getBool( element.get('multipleRegions','false') )
        if regions:
            regions = []
            for region in element.findall('region'):
                resonanceSection = readResolved( region[0] )
                thisRegion = energyInterval( formalism = resonanceSection, **getAttrs( region, exclude=('formalism') ) )
                regions.append( thisRegion )
            RRR = resolved( formalism=None, **getAttrs( element ) )
            RRR.regions = regions
        else:
            formalism = readResolved(element[0])
            RRR = resolved( formalism, **getAttrs( element, exclude=('formalism',) ) )
        xPath.pop()
        return RRR


class unresolved( ancestryModule.ancestry ):

    moniker = 'unresolved'

    def __init__(self, formalism, domainMin, domainMax, domainUnit, reconstructCrossSection=True):
        ancestryModule.ancestry.__init__( self )
        setattr(self, formalism.moniker, formalism)
        self.evaluated = formalism
        if isinstance( self.evaluated, ancestryModule.ancestry ): self.evaluated.setAncestor( self )
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit
        self.reconstructCrossSection = reconstructCrossSection

    
    def toString( self, simpleString = False ):
        return ("Unresolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for L in self.tabulatedWidths.L_values:
            for J in L.J_values:
# BRB6 hardwired
                elist = J.energyDependentWidths.getColumn('energy',self.domainUnit)
                if elist is None:
                    warnings.append( warning.URRmissingEnergyList( L.L, J.J, J ) )
                else:
                    if elist[0] > self.domainMin or elist[-1] < self.domainMax:
                        warnings.append( warning.URRdomainMismatch( L.L, J.J, J ) )
                    missingPoints = [i1 for i1 in range(1,len(elist)) if elist[i1] > 3*elist[i1-1]]
                    for idx in missingPoints:
                        warnings.append( warning.URRinsufficientEnergyGrid( L.L, J.J,
                            PQU.PQU(elist[idx-1],self.domainUnit), PQU.PQU(elist[idx],self.domainUnit), J ) )
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        reconstructString = ' reconstructCrossSection="false"' if not self.reconstructCrossSection else ''
        xmlString = [ '%s<%s domainMin="%s" domainMax="%s" domainUnit="%s" formalism="%s"%s>' %
                ( indent, self.moniker, self.domainMin, self.domainMax, self.domainUnit, self.evaluated.moniker, reconstructString ) ]
        xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        table = unresolvedTabulatedWidths.parseXMLNode( element.find('tabulatedWidths'), xPath, linkData )
        URR = unresolved( formalism = table, **getAttrs( element, exclude=('formalism',) ) )
        xPath.pop()
        return URR


class energyInterval:
    """ resolved region may be made up of multiple energy intervals (deprecated) """
    def __init__(self,index, formalism, domainMin, domainMax, domainUnit):
        self.index = index
        self.evaluated = formalism
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit
    
    def setIndex(self, index):
        self.index = index
    
    def toString(self, simpleString = False):
        return ("%s resonances, %s to %s. Contains %i resonances" % 
                (self.evaluated, self.domainMin, self.domainMax, len(self.evaluated) ) )
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [indent+
                '<region index="%s" domainMin="%r" domainMax="%r" domainUnit="%s" formalism="%s">'
                % ( self.index, self.domainMin, self.domainMax, self.domainUnit, self.evaluated.moniker ) ]
        if self.evaluated:
            xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        else:
            raise Exception( "Resonance section contains no data!" )
        xmlString[-1] += '</region>'
        return xmlString


def _resonance_checker(self,info,things):
    '''
    Common resolved resonance checks

    :param self: a reference to the resonanceFormalismBaseClass (or derived) or RMatrix class
    :param info: the dictionary of checking options, shared by all .check() member functions.
                 For the potential scattering convergence test, the 'potentialScatteringCoverganceCriteria'
                 value is the limit on the ratio of the Lth cross section term to the 0th one.  It is printed
                 as a percent but should be entered as a fraction.
    :param things: list of things to check the obvious way (with a .check() function)
    :return:
    '''
    from fudge.gnd import warning
    warnings=[]

    # check the member data
    for thing in things:
        warningList = thing.check(info)
        if warningList:
            warnings.append(warning.context(thing.moniker,warningList))

    # now for the physics checks
    import fudge.processing.resonances.reconstructResonances as rrReconstructModule
    rrReconstructor = rrReconstructModule.getResonanceReconstructionClass(self.moniker)(info['reactionSuite'])

    # check for mismatched spins
    warnings+=rrReconstructor.setResonanceParametersByChannel(\
        warnOnly=True,\
        useReichMooreApproximation=self.moniker==RM.moniker) #multipleSScheme='NJOY', 'ENDF' or None

    # check for missing channels via sum rule of statistical weights
    sumgJ={}
    for c in rrReconstructor.channels:
        if (c.reaction,c.l) not in sumgJ: sumgJ[(c.reaction,c.l)]=0.0
        sumgJ[(c.reaction,c.l)]+=c.gfact
    keys = sumgJ.keys()
    keys.sort()
    for rxn,L in keys:
        if abs(abs(sumgJ[(rxn,L)]) - abs(2.0*L+1.0)) > 1e-6:
            warnings.append(warning.badSpinStatisticalWeights(L,sumgJ[(rxn,L)],2.*L+1,rxn))

    # check for missing channels by looping through all the allowed angular momenta
    def getSList(Ia,Ib):
        '''Get possible spins'''
        smin=abs(Ia-Ib)
        smax=Ia+Ib
        nS = int(smax-smin)+1
        return [iS+smin for iS in range(nS)]
    LMax = 0
    allSs = []
    rxnList = []
    for c in rrReconstructor.channels:
        if c.eliminated: continue
        LMax = max(c.l,LMax)
        if c.reaction not in rxnList: rxnList.append(c.reaction)
        if c.channelClass not in [rrReconstructModule.FISSIONCHANNEL,rrReconstructModule.COMPETATIVECHANNEL]:
            for s in getSList(*rrReconstructor.getParticleSpins(c.reaction)):
                if s not in allSs: allSs.append(s)
    for L in range(0,LMax+1):
        Jmin=min(allSs)
        Jmax=LMax+max(allSs)
        nJ=int(Jmax-Jmin)
        for iJ in range(0,nJ+1):
            J=iJ+Jmin
            for rxn in rxnList:
                if 'ission' not in rxn:
                    for S in getSList(*rrReconstructor.getParticleSpins(rxn)):
                        if not J in rrReconstructModule.getAllowedTotalSpins(L,S,useFactor2Trick=False): continue
                        for c in rrReconstructor.channels:
                            gotIt = False
                            if rxn==c.reaction and \
                                    rrReconstructModule.spins_equal(c.l,L) and \
                                    rrReconstructModule.spins_equal(c.J,J) and \
                                    rrReconstructModule.spins_equal(c.s,S):
                                gotIt=True
                                break
                        if not gotIt:
                            warnings.append(warning.missingResonanceChannel(L,S,J,rxn))

    # check for convergence in L
    import numpy
    almostXS={}
    for c in rrReconstructor.channels:
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
class resonanceFormalismBaseClass( ancestryModule.ancestry ) :

    optAttrList = ( 'calculateChannelRadius', 'computeAngularDistribution', 'LvaluesNeededForConvergence' )

    def __init__(self, resonanceParameters=None, scatteringRadius=None, **kwargs):
        """LdependentScatteringRadii is a list of dictionaries, one for each L with APL != AP
        only used in ENDF for Reich_Moore form
        
        resonanceParameters should be a gnd.core.math.table.table instance."""

        index = 0
        for attr in self.optAttrList:
            setattr( self, attr, kwargs.get(attr) )
        if self.computeAngularDistribution:
            self.computeAngularDistribution = bool(self.computeAngularDistribution)

        self.resonanceParameters = resonanceParameters or []
        if self.resonanceParameters:
            self.resonanceParameters.setAncestor( self )
        self.scatteringRadius = scatteringRadius
        if self.scatteringRadius: self.scatteringRadius.setAncestor( self )
        ancestryModule.ancestry.__init__( self )

    def check( self, info ):
        return _resonance_checker(self, info, [self.scatteringRadius,self.resonanceParameters])

    def __getitem__(self, idx):
        return self.resonanceParameters.table[idx]
    
    def __len__(self):
        return len(self.resonanceParameters.table)
    
    def addResonance(self, resonance):
        """ insert a new resonance in the resonance parameter table """
        #resonance = (energy, J, l, ... )
        self.resonanceParameters.table.addRow( resonance )

    def toXMLList( self, indent = '', **kwargs ):

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        if not self.moniker :
            raise NotImplementedError ("Please use specific formalisms (MLBW, RM, etc) instead of the BaseClass")
        xmlString = '%s<%s' % ( indent, self.moniker )
        for attr in self.optAttrList:
            if getattr(self,attr,None):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal )
        xmlString = [xmlString+'>']
        xmlString += self.scatteringRadius.toXMLList( indent2, **kwargs )
        if self.resonanceParameters:
            xmlString.extend( self.resonanceParameters.toXMLList( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        radius = scatteringRadius.parseXMLNode( element.find('scatteringRadius'), xPath, linkData )
        linkData['conversionTable'] = {'index':int, 'L':int}    # inform table class how to treat columns
        parameters = resonanceParameters.parseXMLNode( element.find( resonanceParameters.moniker ),
                xPath, linkData )
        attrs = getAttrs( element )
        resonanceData = cls( parameters, radius, **attrs )
        del linkData['conversionTable']
        xPath.pop()
        return resonanceData
    
class SLBW( resonanceFormalismBaseClass ) :
    """Resonance region in Single-Level Breit-Wigner form."""

    moniker = 'SingleLevel_BreitWigner'

    def check(self,info):
        '''
        Common resolved resonance checks

        :param self: a reference to the resonanceFormalismBaseClass (or derived) or RMatrix class
        :param info: the dictionary of checking options, shared by all .check() member functions.
                     For the potential scattering convergence test, the 'potentialScatteringCoverganceCriteria'
                     value is the limit on the ratio of the Lth cross section term to the 0th one.  It is printed
                     as a percent but should be entered as a fraction.
        :return:
        '''
        from fudge.gnd import warning
        warnings=[]

        # check the member data
        for thing in  [self.scatteringRadius,self.resonanceParameters]:
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.context(thing.moniker,warningList))
        from fudge.gnd import warning
        warnings=[]

        # now for the physics checks
        import fudge.processing.resonances.reconstructResonances as rrReconstructModule
        rrReconstructor = rrReconstructModule.getResonanceReconstructionClass(self.moniker)(info['reactionSuite'])
        warnings+=rrReconstructor.setResonanceParametersByChannel(warnOnly=True) #multipleSScheme='NJOY', 'ENDF' or None

        for i,channels in enumerate(rrReconstructor.channels):
            # check for missing channels
            sumgJ={}
            for c in channels:
                if not c.isElastic: continue
                #FIXME: mixing gamma, neutron and fission channels!!!!
                #FIXME: I think gJ should use particle pair spins, not neutron+target!!!!
                if c.l not in sumgJ: sumgJ[c.l]=0.0
                sumgJ[c.l]+=c.gfact
            Ls = sumgJ.keys()
            Ls.sort()
            for L in Ls:
                if abs(abs(sumgJ[L]) - abs(2.0*L+1.0)) > 1e-6:
                    warnings.append(warning.badSpinStatisticalWeights(L,sumgJ[L],2.*L+1))

            # check for convergence in L
            import numpy
            almostXS={}
            for c in channels:
                if c.l not in almostXS:
                    egrid=numpy.array([rrReconstructor.lowerBound, rrReconstructor.upperBound])
                    phis=rrReconstructor.phiByChannel(c,egrid)
                    almostXS[c.l]=pow(numpy.sin(phis)/rrReconstructor.k(egrid),2.0)
            fom=max(almostXS[max(almostXS.keys())])/max(almostXS[min(almostXS.keys())])
            fomTarget=info.get('potentialScatteringCoverganceCriteria',0.001)
            if fom > fomTarget:
                warnings.append(warning.potentialScatteringNotConverged(c.l, rrReconstructor.upperBound, fom, fomTarget))

        return warnings

class MLBW( resonanceFormalismBaseClass ) :
    """Resonance region in Multi-Level Breit-Wigner form."""

    moniker = 'MultiLevel_BreitWigner'

class RM( resonanceFormalismBaseClass ) :
    """Reich_Moore formalism."""

    moniker = 'Reich_Moore'

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

    optAttrList = ('approximation','boundaryCondition','relativisticKinematics','reconstructAngular',
            'reducedWidthAmplitudes','calculatePenetrability','calculateShift','calculateChannelRadius',
            'ENDFconversionFlag')

    def __init__(self, resonanceReactions, spinGroups, **kwargs):

        ancestryModule.ancestry.__init__(self)
        self.resonanceReactions = resonanceReactions
        self.spinGroups = spinGroups
        for sg in self.spinGroups: sg.setAncestor(self, attribute='index')
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

    def check( self, info ):
        return _resonance_checker(self, info, [self.resonanceReactions]+self.spinGroups)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + '  '
        xmlString = ['%s<%s' % ( indent, self.moniker ) ]
        for attr in self.optAttrList:
            if getattr(self,attr):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString[0] += ' %s="%s"' % (attr,attrVal)
        xmlString[0] += '>'
        xmlString += self.resonanceReactions.toXMLList( indent=indent2, **kwargs )
        for gr in self.spinGroups:
            xmlString.append( gr.toXMLList( indent, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        RRs = resonanceReactions()
        RRs.parseXMLNode( element.find(resonanceReactions.moniker), xPath, linkData)
        spinGroups = []
        linkData['conversionTable'] = {'index':int, 'L':int, 'channelSpin':spin}
        for sg in [ch for ch in element if ch.tag=='spinGroup']:
            parameters = resonanceParameters.parseXMLNode( sg.find(resonanceParameters.moniker), xPath, linkData )
            chs = channels()
            chs.parseXMLNode( sg.find(channels.moniker), xPath, linkData )
            spinGroups.append( spinGroup( channels=chs, resonanceParameters=parameters, **getAttrs(sg, required=('background',) ) ) )
        tmp = RMatrix( RRs, spinGroups, **getAttrs(element,
            required=RMatrix.optAttrList ) )
        del linkData['conversionTable']
        xPath.pop()
        return tmp

class resonanceReactions( suitesModule.suite ):

    moniker = 'resonanceReactions'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [resonanceReaction])

    def check( self, info ):
        from fudge.gnd import warning
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
        self.scatteringRadius = scatteringRadius
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

            products = set([p.name for p in self.reactionLink.link.outputChannel.products])
            products.remove( self.ejectile )
            if len(products) != 1: raise ValueError("Cannot compute resonanceReaction residual!")
            self.__residual = products.pop()
        return self.__residual

    def check( self, info ):
        from fudge.gnd import warning
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

    def isFission( self ):

        return self.reactionLink.link.outputChannel.isFission()

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        attrstring = ' computeShiftFactor="true"' if self.computeShiftFactor else ''
        if self.eliminated: attrstring += ' eliminated="true"'
        xmlString = ['%s<%s label="%s" ejectile="%s"%s>' % (indent, self.moniker, self.label, self.ejectile, attrstring) ]
        xmlString += self.reactionLink.toXMLList( indent=indent2, **kwargs )
        if self.Q is not None: xmlString += self.Q.toXMLList( indent=indent2, **kwargs )
        if self.scatteringRadius is not None: xmlString += self.scatteringRadius.toXMLList( indent=indent2, **kwargs)
        if self.hardSphereRadius is not None: xmlString += self.hardSphereRadius.toXMLList( indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        reactionLink = linkModule.link.parseXMLNode( element.find('link'), xPath, linkData )
        Qval, scatRad, hsRad = None, None, None
        if element.find( Q.moniker ):
            Qval = Q.parseXMLNode( element.find( Q.moniker ), xPath, linkData )
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

class spinGroup( ancestryModule.ancestry ):
    """
    Single group with same Jpi (conserved). Each spin group contains an AP (scattering radius),
    along with 1 or more resonance widths.
    """

    moniker = 'spinGroup'

    def __init__(self, index, spin, parity, channels, resonanceParameters, background=None, applyPhaseShift=False):
        ancestryModule.ancestry.__init__(self)
        self.index = index
        self.spin = spin
        self.parity = parity
        self.channels = channels
        self.channels.setAncestor( self )
        self.resonanceParameters = resonanceParameters
        self.resonanceParameters.setAncestor( self )
        self.background = background
        self.applyPhaseShift = applyPhaseShift

    def __getitem__(self, idx):
        return self.resonanceParameters[idx]

    def __len__(self):
        return len(self.resonanceParameters.table)
    
    def __lt__(self, other):
        # for sorting spin groups by Jpi. group J values together
        return (self.spin, self.parity) < (other.spin, other.parity)

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for thing in []:
            if thing is  None: continue
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.context(thing.moniker,warningList))
        return warnings

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        xmlString = '%s<%s index="%i" spin="%s" parity="%s"' % (indent, self.moniker, self.index, self.spin, self.parity)
        for opt in 'background','applyPhaseShift':
            if getattr(self, opt):
                attrVal = getattr(self,opt)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (opt, attrVal)
        xmlString = [xmlString + '>']
        xmlString += self.channels.toXMLList(indent=indent2, **kwargs)
        if self.resonanceParameters.table.columns:    # need not contain any data
            xmlString.extend( self.resonanceParameters.toXMLList( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return '\n'.join(xmlString)

#############################################
#    end of R-Matrix (LRF=7) classes        #
#############################################


####################################################
# remaining classes are for unresolved resonances: #
####################################################
class unresolvedTabulatedWidths:

    moniker = 'tabulatedWidths'
    optAttrList = ('interpolation','ENDFconversionFlag')

    def __init__(self, L_values=None, scatteringRadius=None, **kwargs):
        """ L_values is a list of URR_Lsections, which in turn contains a list of J_sections. 
        The average widths are given in the J_sections """
        index = 0
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))
        self.L_values = L_values or []
        self.scatteringRadius = scatteringRadius
        if self.scatteringRadius: self.scatteringRadius.setAncestor( self )
    
    def __len__(self):
        return len(self.L_values)
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = '%s<%s' % (indent, self.moniker)
        for attr in self.optAttrList:
            if getattr(self,attr,None):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal)
        xmlString = [xmlString+'>']
        xmlString += self.scatteringRadius.toXMLList( indent2, **kwargs )
        for L in self.L_values:
            xmlString += L.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        L_values = []
        scatRadius = scatteringRadius.parseXMLNode( element.find('scatteringRadius'), xPath, linkData )
        linkData['conversionTable'] = {'index':int}
        for lval in element:
            if lval.tag != URR_Lsection.moniker: continue
            J_values = []
            for jval in lval:
                consWid, edepWid = {}, None
                for width in jval:
                    if width.tag=='constantWidths':
                        consWid = getAttrs( width )
                    else: edepWid = tableModule.table.parseXMLNode( width[0], xPath, linkData )
                if edepWid is None: edepWid = tableModule.table()
                for quant in ('energy','levelSpacing','neutronWidth','captureWidth','fissionWidthA','competitiveWidth'):
                    if ( quant not in consWid and quant not in [e.name for e in edepWid.columns] ):
                        consWid[quant] = PQU.PQU(0,'eV')
                J_values.append( URR_Jsection( spin(jval.get('J')), edepWid, consWid,
                    **getAttrs( jval, exclude=("J"), required=("neutronDOF","gammaDOF","fissionDOF","competitiveDOF") ) ) )
            L_values.append( URR_Lsection( int(lval.get("L")), J_values ) )

        table = unresolvedTabulatedWidths( L_values, scatRadius, **getAttrs( element ) )
        del linkData['conversionTable']
        xPath.pop()
        return table

    
# each of the above sections contains one or more 'L' sections, which each contain one or more 'J's:
class URR_Lsection:
    """ unresolved average widths, grouped by L. Contains list of J-values: """
    moniker = 'L_section'
    def __init__(self, L, J_values):
        self.L = L
        self.J_values = J_values
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = ['%s<%s L="%s">' % (indent, self.moniker, self.L)]
        for J in self.J_values:
            xmlString += J.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString


class URR_Jsection:
    """ unresolved average widths, grouped by J. Should contain *either*
     - energyDependentWidths, with a table of widths / level densities as function of energy, or
     - constantWidths + energyGrid (the grid is necessary to tell where to do the reconstruction) """
    moniker = 'J_section'
    optAttrList = ('competitiveDOF','neutronDOF','gammaDOF','fissionDOF')
    def __init__(self, J, energyDependentWidths=None, constantWidths={}, **kwargs):
        self.J = J
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))
        self.constantWidths = constantWidthSection( **constantWidths )
        self.energyDependentWidths = energyDependentWidths
    
    def eliminateRedundantInfo(self):
        """ ENDF often lists more info than necessary. 
        If we encounter section where everything (except energy) is constant, 
        move data from 'energyDependentWidths' to 'constantWidths' """

        allEliminated = False
        edep = self.energyDependentWidths
        for colId in range(edep.nColumns)[::-1]:
            column = edep.columns[colId]
            columnData = edep.getColumn( column.name, unit='eV' )
            if len(set( columnData ) ) == 1:
                setattr( self.constantWidths, column.name, PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( columnData[0] ), column.unit ) )
                [d.pop(colId) for d in edep.data]
                edep.columns.pop(colId)
        for idx, col in enumerate( edep.columns ): col.index = idx  #re-number
        #if edep.nColumns == 1 and edep.columns[0].name == 'energy':
        #    edep.columns, edep.data = [],[] # all widths are constant
        #    allEliminated = True
        return allEliminated
    
    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xmlString = ['%s<%s J="%s"' % (indent,self.moniker,self.J) ]
        for attr in URR_Jsection.optAttrList:
            if getattr(self, attr):
                xmlString[-1] += ' %s="%s"' % (attr, getattr(self,attr))
        xmlString[-1] += '>'
        constwidth = self.constantWidths.toXMLList( indent2, **kwargs )
        if constwidth: xmlString.append(constwidth)
        if self.energyDependentWidths:
# BRB6 Hardwried
            xmlString.extend( ['%s<energyDependentWidths>' % (indent2)] +
                    self.energyDependentWidths.toXMLList( indent3, **kwargs ) )
            xmlString[-1] += '</energyDependentWidths>'
# BRB6 Hardwried
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString


# single region (specified by L, J and optionally energy) for which we provide average widths in URR:
class constantWidthSection:
    """
    for widths or level spacings that depend only on L and J (not on energy)
    """
    optAttrList = ('levelSpacing','neutronWidth','captureWidth','competitiveWidth',
            'fissionWidthA','fissionWidthB')
    def __init__(self, **kwargs):
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr, None))
    
    def toXMLList( self, indent = '', **kwargs ) :

        if not any( [ getattr(self, attr) for attr in self.optAttrList ] ):
            return ''
        xmlString = '%s<constantWidths' % indent
        for attr in self.optAttrList:
            if getattr(self, attr):
                xmlString += ' %s="%s"' % ( attr, getattr(self, attr).toString( keepPeriod = False ) )
        xmlString += '/>'
        return xmlString


# helper functions for reading in from xml:
def getBool( value ):
    return {'true':True, '1':True, 'false':False, '0':False}[value]

def floatOrint( value ):
    if float( value ).is_integer(): return int( value )
    return float( value )

def getAttrs(element, exclude=(), required=()):
    """ read attributes from one xml element. Convert to PhysicalQuantityWithUncertainty (or bool, int, etc)
    where appropriate, and skip anything in the 'exclude' list: """
    conversionTable = {'domainMin':float, 'domainMax':float, 'value':float, 'energy':PQU.PQU,
            'neutronWidth':PQU.PQU, 'captureWidth':PQU.PQU, 'fissionWidthA':PQU.PQU, 'fissionWidthB':PQU.PQU, 'competitiveWidth':PQU.PQU,
            'levelSpacing':PQU.PQU, 'radius':PQU.PQU, 'hardSphereRadius':PQU.PQU, 'channelSpin':spin,
            'reconstructCrossSection':getBool, 'reconstructAngular':getBool, 'multipleRegions': getBool,
            'calculateChannelRadius':getBool, 'computeAngularDistribution':getBool,
            'calculateShift':getBool,'calculatePenetrability':getBool, 'channelIndex':int, 'columnIndex':int,
            'LvaluesNeededForConvergence':int, 'ENDF_MT':int, 'index':int, 'L':int,
            'neutronDOF':floatOrint, 'gammaDOF':floatOrint, 'competitiveDOF':floatOrint, 'fissionDOF':floatOrint,
            'spin':spin, 'parity':parity, 'boundaryConditionOverride':float,
            }
    attrs = dict( element.items() )
    for key in attrs.keys():
        if key in exclude: attrs.pop(key)
        elif key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
    for val in required:
        if val not in attrs: attrs[val] = False
    return attrs
