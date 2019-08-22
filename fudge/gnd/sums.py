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
This module contains the 'sum', 'summand' and 'listOfSummands' classes
"""

import xData.link as linkModule
import xData.ancestry as ancestryModule
from pqu import PQU as PQUModule

from fudge.gnd import warning as warningModule

from channelData import Q as QModule
from fudge.gnd import suites as suitesModule
from fudge.gnd.reactionData import crossSection as crossSectionModule
from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd import documentation as documentationModule

__metaclass__ = type

class sums( ancestryModule.ancestry ):
    """
    Contains all summed quantities. Currently supports summed cross sections and multiplicities,
    could extend to other types of sums later.
    """

    moniker = 'sums'
    ancestryMembers = ('crossSections','multiplicities')

    def __init__(self):
        ancestryModule.ancestry.__init__(self)
        self.crossSections = suitesModule.crossSections()
        self.multiplicities = suitesModule.multiplicities()

    @property
    def crossSections(self): return self.__crossSections

    @crossSections.setter
    def crossSections(self, value):
        if not isinstance(value, suitesModule.crossSections): raise TypeError
        value.setAncestor(self)
        self.__crossSections = value

    @property
    def multiplicities( self ):
        return self.__multiplicities

    @multiplicities.setter
    def multiplicities( self, value ):
        if not isinstance(value, suitesModule.multiplicities): raise TypeError
        value.setAncestor(self)
        self.__multiplicities = value

    def convertUnits(self, unitMap):
        for child in self.ancestryMembers:
            getattr(self, child).convertUnits( unitMap )

    def toXMLList(self, indent='', **kwargs):

        if not any([getattr(self,val) for val in self.ancestryMembers]): return []
        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        for child in self.ancestryMembers:
            val = getattr(self,child)
            if val: xml += val.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag in self.ancestryMembers:
                getattr(self, child.tag).parseXMLNode( child, xPath, linkData )
            else:
                print("Warning: unexpected child '%s' in sums" % child.tag)
        xPath.pop()


class crossSectionSum( ancestryModule.ancestry ):
    """
    Stores a summed quantity (cross section, multiplicity, etc.) along with the list of what's being summed over
    """

    moniker = 'crossSectionSum'
    ancestryMembers = ( 'summands', 'Q', 'crossSection' )

    def __init__( self, label, ENDF_MT, summands = None, crossSection = None, documentation = None ):
        """
        Create a 'sum' instance
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :param summands: list showing what we are summing over
        :param crossSection: should equal the sum of summands' cross sections
        :param documentation:
        :return:
        """

        ancestryModule.ancestry.__init__( self )
        self.label = label
        self.summands = summands or listOfSummands()
        self.summands.setAncestor( self )

        self.Q = QModule.component( )
        self.Q.setAncestor( self )

        self.crossSection = crossSectionModule.component( )
        self.crossSection.setAncestor( self )
        if( crossSection is not None ) : self.crossSection.add( crossSection )

        self.documentation = documentation
        self.ENDF_MT = int( ENDF_MT )

    def __str__( self ):
        return self.label

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.Q.convertUnits( unitMap )
        self.crossSection.convertUnits( unitMap )

    def getENDL_CS_ENDF_MT( self ) :
        """
        Returns the reaction's ENDL C, S and ENDF's MT values as integers in a python dictionary with 
        keys 'C', 'S' and 'MT' (e.g., { 'C' : 11, 'S' : 1, 'MT' : 53 }).
        """

        from fudge.legacy.converting import endf_endl

        MT = self.ENDF_MT
        C, S = endf_endl.getCSFromMT( MT )
        return( { 'C' : C, 'S' : S, 'MT' : MT } )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False,
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False ) :

        return( self.crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
                heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero, addToSuite = True ) )

    def sumSummands(self, label=None):
        """
        Calculate the sum by looping over summands. If a style is specified, use that style if available.
        Otherwise sum the evaluated style.
        :param style: optional string identifying what style of data to sum ('eval', 'recon', etc.)
        :return:
        """
        def getPointwiseLinearForm( component ):
            if label is not None and label in component: return component[ label ]
            else: return component.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)

        sum_ = getPointwiseLinearForm( self.summands[0].link )
        for idx in range(1,len(self.summands)):
            sum_, current = sum_.mutualify( 1e-8, 1e-8, 0,
                    getPointwiseLinearForm( self.summands[idx].link ), 1e-8, 1e-8, 0 )
            sum_ += current
        return sum_

    def scaleSummandsToMatchSum(self, label=None):
        """
        rescale the cross sections in the summands list so that their sum matches the cross section in self
        :return:
        """
        if label is not None and label in self.crossSection:
            correctXsec =self.crossSection[ label ]
        else:
            correctXsec = self.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        summedXsec = self.sumSummands(label=label)
        if correctXsec.domain() != summedXsec.domain():
            raise warningModule.summedCrossSectionDomainMismatch( obj=self )
        scaleFactor=correctXsec/summedXsec

        for s in self.summands:
            label = s.link.evaluated.label
            cutHere=s.link.evaluated.domainMin
            try:
                cutHere=max(cutHere,s.link.evaluated.getAncestor().getAncestor().getThreshold(unit=s.link.evaluated.domainUnit))
            except: # if getting threshold fails, cut at the domainMin
                pass
            answer = s.link.evaluated.returnAsClass(s.link.evaluated,s.link.evaluated*scaleFactor).domainSlice(domainMin=cutHere)
            answer.label=label
            s.link.remove(label)
            s.link.add(answer)


    def scaleSumToMatchSummands(self, label=None):
        """
        Sum the summands cross section list, then set the sum to the cross section in self.
        If a style is specified, the sum will be stored under that style. Otherwise it replaces the
        evaluated style.
        :return:
        """
        resummed = self.sumSummands()
        if label is None: label = self.crossSection.evaluated.label
        resummed.label = resummed.style = label
        if label == self.crossSection.evaluated.label:
            self.crossSection.remove( label )
        self.crossSection.add( resummed )

    def check( self, info ) :

        warnings = self.crossSection.check( info )

        # does self.crossSection equal the sum over summand reaction cross sections?
        sum_ = self.sumSummands(label=info['reconstructedStyleName'])
        if info['reconstructedStyleName'] in self.crossSection:
            quotedXsec =self.crossSection[ info['reconstructedStyleName'] ]
        else:
            quotedXsec = self.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sum_.domain() != quotedXsec.domain():
            warnings.append( warningModule.summedCrossSectionDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sum_
            diff.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax with 'inf'
                if relativeDiff.rangeMax > info['crossSectionMaxDiff']:
                    warnings.append( warningModule.summedCrossSectionMismatch( relativeDiff.rangeMax * 100, obj=self ) )
            except TypeError:
                warnings.append( warningModule.summedCrossSectionZeroDivision( obj=self ) )

        return warnings

    def getLabel(self):

        return self.label

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ""
        xmlString = [ '%s<%s label="%s" ENDF_MT="%s"%s>' % ( indent, self.moniker, self.label, self.ENDF_MT, attrStr ) ]
        xmlString += self.summands.toXMLList( indent2, **kwargs )
        xmlString += self.Q.toXMLList( indent2, **kwargs )
        xmlString += self.crossSection.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )
        attributes = {}
        for key in element.keys( ) :
            if( key in ( 'ENDF_MT', 'label' ) ) :
                attributes[key] = element.get( key )
            elif( key == 'Q' ) :
                attributes[ key ] = PQUModule.PQU( element.get(key) )
            else:
                raise ValueError( 'Unsupported attribute "%s"' % key )
        MT = int( attributes.pop( 'ENDF_MT' ) )
        label = attributes.pop( 'label' )

        summands_, documentation = None, None
        for child in element:
            if child.tag == listOfSummands.moniker:
                summands_ = listOfSummands.parseXMLNode( child, xPath, linkData )
            elif child.tag == documentationModule.documentation.moniker:
                documentation = documentationModule.documentation.parseXMLNode( child, xPath, linkData )
            elif child.tag in ( QModule.component.moniker, crossSectionModule.component.moniker ):
                pass    # these are read below
            else:
                raise TypeError( 'Unexpected element "%s" encountered' % child.tag )
        sum_ = crossSectionSum( label, MT, summands=summands_, documentation=documentation )
        sum_.crossSection.parseXMLNode( element.find(crossSectionModule.component.moniker), xPath, linkData )
        sum_.Q.parseXMLNode( element.find(QModule.component.moniker), xPath, linkData )
        xPath.pop()
        return sum_


class multiplicitySum( ancestryModule.ancestry ):
    """
    Stores a summed multiplicity along with the list of what's being summed over
    """

    moniker = 'multiplicitySum'
    ancestryMembers = ( 'summands', 'multiplicity' )

    def __init__( self, label, ENDF_MT, summands = None, documentation = None, attributes = None ):
        """
        Create a 'sum' instance
        :param name: descriptive string
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :param summands: list showing what we are summing over
        :param documentation:
        :param attributes: additional info
        :return:
        """

        ancestryModule.ancestry.__init__( self )
        self.label = label

        self.summands = summands or listOfSummands()
        self.summands.setAncestor( self )

        self.multiplicity = multiplicityModule.component( )
        self.multiplicity.setAncestor( self )

        self.documentation = documentation
        self.attributes = attributes or {}
        self.ENDF_MT = int( ENDF_MT )

    def __str__( self ):
        return self.label

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.multiplicity.convertUnits( unitMap )

    def check( self, info ) :

        warnings = []

        if self.multiplicity.isConstant() and self.multiplicity.getConstant() < 1:
            warnings.append( warningModule.negativeMultiplicity( self.getConstant(), self ) )
        else:   # energy-dependent mult.
            for form in self.multiplicity :
                if hasattr(form, 'rangeMin') and form.rangeMin < 0:
                    warnings.append( warningModule.negativeMultiplicity( form.rangeMin, obj=form ) )

        # does multiplicity equal the sum over its summand multiplicities?
        sum_ = self.summands[0].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        for idx in range(1,len(self.summands)):
            sum_, current = sum_.mutualify( 1e-8, 1e-8, 0,
                    self.summands[idx].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8), 1e-8, 1e-8, 0 )
            sum_ += current
        quotedXsec = self.multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sum_.domain() != quotedXsec.domain():
            warnings.append( warningModule.summedMultiplicityDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sum_
            diff.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax with 'inf'
                if relativeDiff.rangeMax > info['multiplicityMaxDiff']:
                    warnings.append( warningModule.summedMultiplicityMismatch( relativeDiff.rangeMax * 100, obj=self ) )
            except TypeError:
                warnings.append( warningModule.summedMultiplicityZeroDivision( obj=self ) )

        return warnings

    def getLabel(self):

        return self.label

    def setAttribute( self, name, value ):

        self.attributes[name] = value

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ""
        for attr in sorted(self.attributes): attrStr += ' %s="%s"' % (attr, self.attributes[attr])
        xmlString = [ '%s<%s label="%s" ENDF_MT="%s"%s>' % ( indent, self.moniker, self.label, self.ENDF_MT, attrStr ) ]
        xmlString += self.summands.toXMLList( indent2, **kwargs )
        xmlString += self.multiplicity.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )
        attributes = {}
        for key in element.keys( ) :
            if( key in ( 'ENDF_MT', 'label' ) ) :
                attributes[key] = element.get( key )
            elif( key == 'Q' ) :
                attributes[ key ] = PQUModule.PQU( element.get(key) )
            else:
                raise ValueError( 'Unsupported attribute "%s"' % key )
        MT = int( attributes.pop( 'ENDF_MT' ) )
        label = attributes.pop( 'label' )

        summands_, documentation = None, None
        for child in element:
            if child.tag == listOfSummands.moniker:
                summands_ = listOfSummands.parseXMLNode( child, xPath, linkData )
            elif child.tag == documentationModule.documentation.moniker:
                documentation = documentationModule.documentation.parseXMLNode( child, xPath, linkData )
            elif child.tag == multiplicityModule.component.moniker:
                pass    # to be read below
            else:
                raise TypeError( 'Unexpected element "%s" encountered' % child.tag )

        sum_ = multiplicitySum( label, MT, summands=summands_, documentation=documentation,
                attributes = attributes )
        sum_.multiplicity.parseXMLNode( element.find( multiplicityModule.component.moniker ), xPath, linkData )

        xPath.pop()
        return sum_

class listOfSummands( ancestryModule.ancestry ):

    moniker = 'summands'
    ancestryMembers = ( '', )

    def __init__(self, summandList=None):

        ancestryModule.ancestry.__init__( self )
        self.__summandList = summandList or []

    def __getitem__(self, index):

        return self.summandList[index]

    def __len__(self):

        return len(self.summandList)

    @property
    def summandList(self): return self.__summandList

    def toXMLList( self, indent = '', **kwargs ) :

        xmlList = ['%s<%s>' % (indent, self.moniker)]
        for summand_ in self.summandList:
            xmlList.append( '  %s%s' % (indent, summand_.toXML() ) )
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        summandList = []
        for summand_ in element:
            if summand_.tag == add.moniker:
                summandList.append( add.parseXMLNode(summand_, xPath, linkData) )
            elif summand_.tag == subtract.moniker:
                summandList.append( subtract.parseXMLNode(summand_, xPath, linkData) )
        newList = listOfSummands( summandList=summandList )
        xPath.pop()
        return newList


class add( linkModule.link ):
    """ link representing one of the quantities that is added to the sum """

    moniker = 'add'

class subtract( linkModule.link ):
    """ link representing one of the quantities that is subtracted from the sum """

    moniker = 'subtract'
