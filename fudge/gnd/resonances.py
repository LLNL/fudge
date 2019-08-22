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

import math
from pqu import PQU
from xData import XYs as XYsModule, link as linkModule, table as tableModule
from fudge.gnd import xParticle, baseClasses
from fudge.gnd import suites as suitesModule

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

    def __init__(self, scatteringRadius=None, resolved=None, unresolved=None, reconstructCrossSection=True):

        ancestryModule.ancestry.__init__( self )

        if (scatteringRadius is not None) and (resolved is not None or unresolved is not None):
            raise TypeError, ("Resonances should contain scattering radius OR resonance parameters, not both")
        self.scatteringRadius = scatteringRadius
        self.resolved = resolved
        self.unresolved = unresolved

        for child in (self.scatteringRadius, self.resolved, self.unresolved):
            if child is not None: child.setAncestor( self )
        self.reconstructCrossSection = reconstructCrossSection
        
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
    
        xmlString = [ '%s<%s reconstructCrossSection="%s">' % 
                ( indent, self.moniker, str(self.reconstructCrossSection).lower() ) ]
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
            bounds.append( (self.scatteringRadius.lowerBound, self.scatteringRadius.upperBound) )
        if self.resolved:
            if self.resolved.multipleRegions:
                bounds += [(reg.lowerBound, reg.upperBound) for reg in self.resolved.regions]
            else: bounds.append( (self.resolved.lowerBound, self.resolved.upperBound) )
        if self.unresolved:
            bounds.append( (self.unresolved.lowerBound, self.unresolved.upperBound) )

        for idx in range(len(bounds)-1):
            assert bounds[idx][1] == bounds[idx+1][0], "Resonance region boundaries don't match!"

        if( asPQU ):
            return (bounds[0][0], bounds[-1][1])
        elif unitTo:
            return (bounds[0][0].getValue(unitTo), bounds[-1][1].getValue(unitTo))
        else:
            return (bounds[0][0].value, bounds[-1][1].value)

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )

        scatRadius, RRR, URR = None,None,None
        kReconstruct = (element.get("reconstructCrossSection") == "true")
        for child in element:
            if child.tag==scatteringRadius.moniker:
                scatRadius = scatteringRadius.parseXMLNode( child, xPath, linkData )
            elif child.tag==resolved.moniker:
                RRR = resolved.parseXMLNode( child, xPath, linkData )
            elif child.tag==unresolved.moniker:
                URR = unresolved.parseXMLNode( child, xPath, linkData )
            else:
                raise Exception("unknown element '%s' encountered in resonances!" % child.tag)

        res = resonances( scatteringRadius = scatRadius, resolved=RRR, unresolved=URR,
                reconstructCrossSection=kReconstruct )
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

    def check( self, info ):
        return []

    def isEnergyDependent(self):
        return isinstance( self.form, XYsModule.XYs1d )

    def isLdependent(self):
        return isinstance(self.form, LdependentScatteringRadii)

    def getValueAs( self, unit, energy_grid=None, L=None ):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise Exception("Missing: energy grid to evaluate E-dependent scattering radius")
            energy_unit = self.form.axes[-1].unit
            return [self.form.evaluate( PQU.PQU( e, energy_unit ), unitTo = unit ) for e in energy_grid]
        elif self.isLdependent():
            if L in self.form.lvals:
                return self.form.lvals[L].getValueAs( unit )
            return self.form.default.getValueAs( unit )
        else:
            return self.form.value.getValueAs( unit )

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
            constantScatteringRadius.moniker: constantScatteringRadius,
            LdependentScatteringRadii.moniker: LdependentScatteringRadii,
            XYsModule.XYs1d.moniker: XYsModule.XYs1d,
        }[ element[0].tag ].parseXMLNode( element[0], xPath, linkData )
        SR = cls( form )
        xPath.pop()
        return SR


class effectiveRadius( scatteringRadius ):

    moniker = 'effectiveRadius'


class constantScatteringRadius( baseClasses.formBase ):

    moniker = 'constant'

    def __init__( self, value, bounds = None, label = None ) :
        if isinstance( value, str ):
            value = PQU.PQU( value )
        self.value = value
        self.bounds = bounds

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
        if self.bounds is not None:
            boundsString = ' lowerBound="%s" upperBound="%s"' % self.bounds
        return [ '%s<%s%s value="%s"%s/>' % ( indent, self.moniker, attributeStr, self.value, boundsString )]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        value = PQU.PQU(element.get('value'))
        bounds = None
        if 'lowerBound' in element.keys():
            bounds = tuple([PQU.PQU(element.get( bound )) for bound in ('lowerBound','upperBound') ])
        CSR = constantScatteringRadius( value, bounds=bounds, label = element.get( 'label' ) )
        xPath.pop()
        return CSR


class LdependentScatteringRadii( baseClasses.formBase ):

    moniker = 'Ldependent'

    def __init__( self, default, lvals, label = None ) :
        self.default = default
        self.lvals = lvals

        if( label is not None ) :
            if not isinstance( label, str ):
                raise TypeError( 'label must be a string' )
        self.__label = label

    @property
    def label( self ) :

        return( self.__label )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        xml = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        xml.append( indent2 + '<default value="%s"/>' % self.default )
        for lval in sorted(self.lvals.keys()):
            xml.append( indent2 + '<Lspecific L="%d" value="%s"/>' % (lval, self.lvals[lval]) )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        default = PQU.PQU( element.find('default').get('value') )
        lvals = {}
        for child in element:
            if child.tag == 'Lspecific':
                lval = int(child.get('L'))
                lvals[lval] = PQU.PQU( child.get('value') )
        result = LdependentScatteringRadii( default, lvals, element.get( 'label' ) )
        xPath.pop()
        return result


class resonanceParameters( ancestryModule.ancestry ):   # FIXME: still needed? It's an extra level of nesting...
    """
    Light-weight wrapper around a table.
    """

    moniker = 'resonanceParameters'

    def __init__(self, table):
        self.table = table
        self.table.setAncestor(self)

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

    def __init__(self, formalism=None, lowerBound='', upperBound='', multipleRegions=False):

        ancestryModule.ancestry.__init__( self )
        if not multipleRegions:
            setattr(self, formalism.moniker, formalism)
        self.evaluated = formalism
        if self.evaluated is not None: self.evaluated.setAncestor( self )
        self.lowerBound = lowerBound
        self.upperBound = upperBound
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
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s' % ( indent, self.moniker ) ]
        if self.multipleRegions:
            xmlString[0] += ' multipleRegions="true">'
            for region in self.regions:
                xmlString += region.toXMLList( indent2, **kwargs )
        else:
            lowerBound = self.lowerBound.toString( keepPeriod = False )
            upperBound = self.upperBound.toString( keepPeriod = False )
            xmlString[0] += ' lowerBound="%s" upperBound="%s" formalism="%s">' % (
                    lowerBound, upperBound, self.evaluated.moniker )
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
            RRR = resolved( multipleRegions = True )
            RRR.regions = regions
        else:
            formalism = readResolved(element[0])
            RRR = resolved( formalism, **getAttrs( element, exclude=('formalism',) ) )
        xPath.pop()
        return RRR


class unresolved( ancestryModule.ancestry ):

    moniker = 'unresolved'

    def __init__(self, formalism, lowerBound, upperBound):
        ancestryModule.ancestry.__init__( self )
        setattr(self, formalism.moniker, formalism)
        self.evaluated = formalism
        if isinstance( self.evaluated, ancestryModule.ancestry ): self.evaluated.setAncestor( self )
        self.lowerBound = lowerBound
        self.upperBound = upperBound
    
    def toString( self, simpleString = False ):
        return ("Unresolved resonances in %s form\n" % self.evaluated.moniker )

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for L in self.tabulatedWidths.L_values:
            for J in L.J_values:
                elist = J.energyDependentWidths.getColumn('energy','eV')
                if elist is None:
                    warnings.append( warning.URRmissingEnergyList( L.L, J.J.value, J ) )
                else:
                    if elist[0] > self.lowerBound.getValueAs('eV') or elist[-1] < self.upperBound.getValueAs('eV'):
                        warnings.append( warning.URRdomainMismatch( L.L, J.J.value, J ) )
                    missingPoints = [i1 for i1 in range(1,len(elist)) if elist[i1] > 3*elist[i1-1]]
                    for idx in missingPoints:
                        warnings.append( warning.URRinsufficientEnergyGrid( L.L, J.J.value, PQU.PQU(elist[idx-1],'eV'),
                            PQU.PQU(elist[idx],'eV'), J ) )
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        lowerBound = self.lowerBound.toString( keepPeriod = False )
        upperBound = self.upperBound.toString( keepPeriod = False )
        xmlString = [ '%s<%s lowerBound="%s" upperBound="%s" formalism="%s">' %
                ( indent, self.moniker, lowerBound, upperBound, self.evaluated.moniker ) ]
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
    def __init__(self,index, formalism, lowerBound, upperBound):
        self.index = index
        self.evaluated = formalism
        self.lowerBound = lowerBound
        self.upperBound = upperBound
    
    def setIndex(self, index):
        self.index = index
    
    def toString(self, simpleString = False):
        return ("%s resonances, %s to %s. Contains %i resonances" % 
                (self.evaluated, self.lowerBound, self.upperBound, len(self.evaluated) ) )
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        lowerBound = self.lowerBound.toString( keepPeriod = False )
        upperBound = self.upperBound.toString( keepPeriod = False )
        xmlString = [indent+
                '<region index="%s" lowerBound="%s" upperBound="%s" formalism="%s">'
                % ( self.index, lowerBound, upperBound, self.evaluated.moniker ) ]
        if self.evaluated:
            xmlString += self.evaluated.toXMLList( indent2, **kwargs )
        else:
            raise Exception( "Resonance section contains no data!" )
        xmlString[-1] += '</region>'
        return xmlString


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
    """R-Matrix Limited formalism. """

    moniker = 'R_Matrix_Limited'

    optAttrList = ('approximation','boundaryCondition','relativisticKinematics',
            'reducedWidthAmplitudes','calculatePenetrability','calculateShift','calculateChannelRadius')

    def __init__(self, openChannels, spinGroups, **kwargs):
        """R-matrix is sorted by 'spin groups', each with same Jpi (conserved)."""

        ancestryModule.ancestry.__init__(self)
        self.channels = openChannels
        self.spinGroups = spinGroups
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + '  '
        xmlString = ['%s<%s' % ( indent, self.moniker ) ]
        for attr in self.optAttrList:
            if getattr(self,attr):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString[0] += ' %s="%s"' % (attr,attrVal)
        xmlString[0] += '>'
        xmlString += self.channels.toXMLList( indent=indent2, **kwargs )
        for gr in self.spinGroups:
            xmlString.append( gr.toXMLList( indent, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        channels = openChannels()
        channels.parseXMLNode( element.find(openChannels.moniker), xPath, linkData )
        spinGroups = []
        linkData['conversionTable'] = {'index':int, 'L':int, 'channelSpin':float}
        for sg in [ch for ch in element if ch.tag=='spinGroup']:
            parameters = resonanceParameters.parseXMLNode( sg.find(resonanceParameters.moniker), xPath, linkData )
            spinGroups.append( spinGroup( resonanceParameters=parameters, **getAttrs(sg, required=('background',) ) ) )
            if sg.find( overrides.moniker ):
                spinGroups[-1].overrides.parseXMLNode( sg.find(overrides.moniker), xPath, linkData )
        tmp = RMatrix( channels, spinGroups, **getAttrs(element,
            required=RMatrix.optAttrList ) )
        del linkData['conversionTable']
        xPath.pop()
        return tmp

class openChannels( suitesModule.suite ):

    moniker = 'channels'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [openChannel] )

class openChannel( ancestryModule.ancestry ):
    """
    Describes one reaction channel that opens up in the resonance region. In an R-Matrix section,
    all open reaction channels should be described in the list of openChannel elements
    """

    moniker = 'channel'

    optAttrList = ('Q', 'boundaryCondition', 'calculatePenetrability','calculateShift')
    def __init__(self, label, reactionLink, name, ENDF_MT, scatteringRadius=None, effectiveRadius=None, **kwargs):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.reactionLink = reactionLink
        self.name = name
        self.ENDF_MT = ENDF_MT
        self.scatteringRadius = scatteringRadius
        self.effectiveRadius = effectiveRadius
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr,None))
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        xmlString = ['%s<%s label="%s" name="%s" ENDF_MT="%s"' % (indent, self.moniker, self.label, self.name, self.ENDF_MT) ]
        for attr in self.optAttrList:
            if getattr(self,attr):
                if( attr in [ 'Q' ] ) :
                    xmlString[-1] += ' %s="%s"' % ( attr, getattr( self, attr ).toString( keepPeriod = False ) )
                else :
                    attrVal = getattr(self,attr)
                    if type(attrVal) is bool: attrVal = str(attrVal).lower()
                    xmlString[-1] += ' %s="%s"' % (attr, attrVal)
        xmlString[-1] += '>'
        xmlString += self.reactionLink.toXMLList( indent=indent2, **kwargs )
        if self.scatteringRadius is not None: xmlString += self.scatteringRadius.toXMLList( indent=indent2, **kwargs)
        if self.effectiveRadius is not None: xmlString += self.effectiveRadius.toXMLList( indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        reactionLink = linkModule.link.parseXMLNode( element.find('link'), xPath, linkData )
        attrs = getAttrs( element )
        if 'Q' not in attrs : attrs['Q'] = PQU.PQU(0,'eV')
        scatRad, effRad = None, None
        if element.find( scatteringRadius.moniker ):
            scatRad = scatteringRadius.parseXMLNode( element.find( scatteringRadius.moniker ), xPath, linkData )
        if element.find( effectiveRadius.moniker ):
            effRad = effectiveRadius.parseXMLNode( element.find( effectiveRadius.moniker ), xPath, linkData )
        tmp = openChannel( attrs.pop('label'), reactionLink=reactionLink, scatteringRadius=scatRad,
                effectiveRadius=effRad, **attrs )
        xPath.pop()
        return tmp

class overrides( suitesModule.suite ):

    moniker = 'overrides'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [channelOverride] )

class channelOverride( ancestryModule.ancestry ):
    """
    Use if necessary to override the true scattering radius, effective radius, etc. for a single
    spin-group / channel combination.
    """

    moniker = 'channelOverride'

    def __init__(self, label, scatteringRadius=None, effectiveRadius=None, penetrability=None,
                 shiftFactor=None, phaseShift=None):

        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.scatteringRadius = scatteringRadius
        self.effectiveRadius = effectiveRadius
        self.penetrability = penetrability
        self.shiftFactor = shiftFactor
        self.phaseShift = phaseShift

    def __nonzero__(self):
        return (self.scatteringRadius is not None or self.effectiveRadius is not None or self.penetrability is not None
                or self.shiftFactor is not None or self.phaseShift is not None)

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent+'  '
        xmlString = ['%s<%s label="%s">' % (indent, self.moniker, self.label) ]
        for attr in ('scatteringRadius', 'effectiveRadius', 'penetrability', 'shiftFactor', 'phaseShift'):
            if getattr(self, attr) is not None:
                xmlString += getattr(self, attr).toXMLList( indent=indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        xPath.append( element.tag )
        data = {}
        for child in element:
            childClass = {'scatteringRadius': scatteringRadius,
                          'effectiveRadius': effectiveRadius,
                          'penetrability': None,    # FIXME last three not yet implemented
                          'shiftFactor': None,
                          'phaseShift': None}.get( child.tag )
            data[ child.tag ] = childClass.parseXMLNode( child, xPath, linkData )
        tmp = channelOverride( element.get('label'), **data )
        xPath.pop()
        return tmp

class spinGroup( ancestryModule.ancestry ):
    """
    Single group with same Jpi (conserved). Each spin group contains an AP (scattering radius),
    along with 1 or more resonance widths.
    """

    moniker = 'spinGroup'

    def __init__(self, index=None, spin=None, parity=None, background=None, applyPhaseShift=False,
            resonanceParameters=None):
        ancestryModule.ancestry.__init__(self)
        self.index = index
        self.spin = spin
        self.parity = parity
        self.background = background
        self.applyPhaseShift = applyPhaseShift
        self.overrides = overrides()
        self.resonanceParameters = resonanceParameters
    
    def __getitem__(self, idx):
        return self.resonanceParameters[idx]
    
    def __len__(self):
        return len(self.resonanceParameters.table)
    
    def __lt__(self, other):
        # for sorting spin groups by Jpi. group J values together
        return (self.spin,self.parity) < (
                other.spin,other.parity)
    
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
        #for ch in self.channelInfo:
        #    xmlString.append( ch.toXMLList( indent2, **kwargs ) )
        xmlString += self.overrides.toXMLList(indent=indent2, **kwargs)
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
    optAttrList = ('interpolation','forSelfShieldingOnly','ENDFconversionFlag')

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
                J_values.append( URR_Jsection( xParticle.spin(jval.get('J')), edepWid, consWid,
                    **getAttrs( jval, exclude=("J"), required=("neutronDOF","gammaDOF","fissionDOF","competitiveDOF") ) ) )
            L_values.append( URR_Lsection( int(lval.get("L")), J_values ) )

        table = unresolvedTabulatedWidths( L_values, scatRadius, **getAttrs( element, required=('forSelfShieldingOnly',) ) )
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
            columnData = edep.getColumn( column.name, units='eV' )
            if len(set( columnData ) ) == 1:
                setattr( self.constantWidths, column.name, PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( columnData[0] ), column.units ) )
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

        xmlString = ['%s<J_section J="%s"' % (indent,self.J) ]
        for attr in URR_Jsection.optAttrList:
            if getattr(self, attr):
                xmlString[-1] += ' %s="%s"' % (attr, getattr(self,attr))
        xmlString[-1] += '>'
        constwidth = self.constantWidths.toXMLList( indent2, **kwargs )
        if constwidth: xmlString.append(constwidth)
        if self.energyDependentWidths:
            xmlString.extend( ['%s<energyDependentWidths>' % (indent2)] +
                    self.energyDependentWidths.toXMLList( indent3, **kwargs ) )
            xmlString[-1] += '</energyDependentWidths>'
        xmlString[-1] += '</J_section>'
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
    conversionTable = {'lowerBound':PQU.PQU, 'upperBound':PQU.PQU, 'value':PQU.PQU, 'energy':PQU.PQU,
            'neutronWidth':PQU.PQU, 'captureWidth':PQU.PQU, 'fissionWidthA':PQU.PQU, 'fissionWidthB':PQU.PQU, 'competitiveWidth':PQU.PQU,
            'levelSpacing':PQU.PQU, 'Q':PQU.PQU, 'radius':PQU.PQU, 'effectiveRadius':PQU.PQU,
            'reconstructCrossSection':getBool, 'multipleRegions': getBool, 'LdependentScatteringRadii': getBool,
            'calculateChannelRadius':getBool, 'computeAngularDistribution':getBool, 'forSelfShieldingOnly': getBool,
            'calculateShift':getBool,'calculatePenetrability':getBool,
            'LvaluesNeededForConvergence':int, 'ENDF_MT':int, 'index':int, 'L':int,
            'neutronDOF':floatOrint, 'gammaDOF':floatOrint, 'competitiveDOF':floatOrint, 'fissionDOF':floatOrint,
            'spin':xParticle.spin, 'parity':xParticle.parity,
            'scatteringRadius':(lambda foo: scatteringRadius(PQU.PQU(foo)) if foo!='energyDependent' else foo),
            }
    attrs = dict( element.items() )
    for key in attrs.keys():
        if key in exclude: attrs.pop(key)
        elif key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
    for val in required:
        if val not in attrs: attrs[val] = False
    return attrs

