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

import fudge
import xData.link as linkModule
import xData.ancestry as ancestryModule
from pqu import PQU as PQUModule

from channelData import Q as QModule
from fudge.gnd.reactionData import crossSection as crossSectionModule
from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd import channelData as channelDataModule
from fudge.gnd import documentation as documentationModule

__metaclass__ = type

class crossSectionSum( ancestryModule.ancestry ):
    """
    Stores a summed quantity (cross section, multiplicity, etc.) along with the list of what's being summed over
    """

    moniker = 'crossSectionSum'

    def __init__( self, name, label, ENDF_MT, summands = None, crossSection = None, documentation = None, date = None ):
        """
        Create a 'sum' instance
        :param name: descriptive string
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :param summands: list showing what we are summing over
        :param crossSection: should equal the sum of summands' cross sections
        :param documentation:
        :return:
        """

        ancestryModule.ancestry.__init__( self )
        self.name = name
        self.label = label
        self.summands = summands or listOfSummands()
        self.summands.setAncestor( self )

        self.Q = QModule.component( )
        self.Q.setAncestor( self )

        self.crossSection = crossSectionModule.component( )
        self.crossSection.setAncestor( self )
        if( crossSection is not None ) : self.crossSection.add( crossSection )

        self.documentation = documentation
        self.date = date
        self.ENDF_MT = int( ENDF_MT )

    def __str__( self ):
        return self.name

    def getENDL_CS_ENDF_MT( self ) :
        """
        Returns the reaction's ENDL C, S and ENDF's MT values as integers in a python dictionary with 
        keys 'C', 'S' and 'MT' (e.g., { 'C' : 11, 'S' : 1, 'MT' : 53 }).
        """

        from fudge.legacy.converting import endf_endl

        MT = self.ENDF_MT
        C, S = endf_endl.getCSFromMT( MT )
        return( { 'C' : C, 'S' : S, 'MT' : MT } )

    def check( self, info ) :

        from fudge.gnd import warning
        warnings = self.crossSection.check( info )

        # does self.crossSection equal the sum over summand reaction cross sections?
        def getPointwiseLinearForm( component ):

            if info['reconstructedStyle'] in component: return component[ info['reconstructedStyle'] ]
            else: return component.toPointwise_withLinearXYs()

        sum_ = getPointwiseLinearForm( self.summands[0].link )
        for idx in range(1,len(self.summands)):
            sum_, current = sum_.mutualify( 1e-8, 1e-8, 0,
                    getPointwiseLinearForm( self.summands[idx].link ), 1e-8, 1e-8, 0 )
            sum_ += current
        quotedXsec = getPointwiseLinearForm( self.crossSection )
        if sum_.domain() != quotedXsec.domain():
            warnings.append( warning.summedCrossSectionDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sum_
            diff.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax() with 'inf'
                if relativeDiff.rangeMax() > info['crossSectionMaxDiff']:
                    warnings.append( warning.summedCrossSectionMismatch( relativeDiff.rangeMax()*100, obj=self ) )
            except TypeError:
                warnings.append( warning.summedCrossSectionZeroDivision( obj=self ) )

        return warnings

    def getLabel(self):

        return self.label

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ""
        if( self.date is not None ) : attrStr += ' date="%s"' % self.date
        xmlString = [ '%s<%s label="%s" name="%s" ENDF_MT="%s"%s>' % ( indent, self.moniker, self.label, self.name, self.ENDF_MT, attrStr ) ]
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
            if( key in ( 'name', 'date', 'ENDF_MT', 'label' ) ) :
                attributes[key] = element.get( key )
            elif( key == 'Q' ) :
                attributes[ key ] = PQUModule.PQU( element.get(key) )
            else:
                raise ValueError( 'Unsupported attribute "%s"' % key )
        MT = int( attributes.pop( 'ENDF_MT' ) )
        name = attributes.pop('name')
        label = attributes.pop( 'label' )
        date = attributes.get( 'date', None )

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
        sum_ = crossSectionSum( name, label, MT, summands=summands_, documentation=documentation, date = date )
        sum_.crossSection.parseXMLNode( element.find(crossSectionModule.component.moniker), xPath, linkData )
        sum_.Q.parseXMLNode( element.find(QModule.component.moniker), xPath, linkData )
        xPath.pop()
        return sum_


class multiplicitySum( ancestryModule.ancestry ):
    """
    Stores a summed multiplicity along with the list of what's being summed over
    """

    moniker = 'multiplicitySum'

    def __init__( self, name, label, ENDF_MT, summands = None, documentation = None, attributes = None ):
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
        self.name = name
        self.label = label

        self.summands = summands or listOfSummands()
        self.summands.setAncestor( self )

        self.multiplicity = multiplicityModule.component( )
        self.multiplicity.setAncestor( self )

        self.documentation = documentation
        self.attributes = attributes or {}
        self.ENDF_MT = int( ENDF_MT )

    def __str__( self ):
        return self.name

    def check( self, info ) :

        from fudge.gnd import warning
        warnings = []

        if self.multiplicity.isConstant() and self.multiplicity.getConstant() < 1:
            warnings.append( warning.negativeMultiplicity( self.getConstant(), self ) )
        else:   # energy-dependent mult.
            for form in self.multiplicity :
                if hasattr(form, 'rangeMin') and form.rangeMin() < 0:
                    warnings.append( warning.negativeMultiplicity( form.rangeMin(), obj=form ) )

        # does multiplicity equal the sum over its summand multiplicities?
        sum_ = self.summands[0].link.toPointwise_withLinearXYs()
        for idx in range(1,len(self.summands)):
            sum_, current = sum_.mutualify( 1e-8, 1e-8, 0,
                    self.summands[idx].link.toPointwise_withLinearXYs(), 1e-8, 1e-8, 0 )
            sum_ += current
        quotedXsec = self.multiplicity.toPointwise_withLinearXYs()
        if sum_.domain() != quotedXsec.domain():
            warnings.append( warning.summedMultiplicityDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sum_
            diff.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax() with 'inf'
                if relativeDiff.rangeMax() > info['multiplicityMaxDiff']:
                    warnings.append( warning.summedMultiplicityMismatch( relativeDiff.rangeMax()*100, obj=self ) )
            except TypeError:
                warnings.append( warning.summedMultiplicityZeroDivision( obj=self ) )

        return warnings

    def getLabel(self):

        return self.label

    def setAttribute( self, name, value ):

        self.attributes[name] = value

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ""
        for attr in sorted(self.attributes): attrStr += ' %s="%s"' % (attr, self.attributes[attr])
        xmlString = [ '%s<%s label="%s" name="%s" ENDF_MT="%s"%s>' % ( indent, self.moniker, self.label, self.name, self.ENDF_MT, attrStr ) ]
        xmlString += self.summands.toXMLList( indent2, **kwargs )
        xmlString += self.multiplicity.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )
        attributes = {}
        for key in element.keys( ) :
            if( key in ( 'name', 'date', 'ENDF_MT', 'label' ) ) :
                attributes[key] = element.get( key )
            elif( key == 'Q' ) :
                attributes[ key ] = PQUModule.PQU( element.get(key) )
            else:
                raise ValueError( 'Unsupported attribute "%s"' % key )
        MT = int( attributes.pop( 'ENDF_MT' ) )
        name = attributes.pop('name')
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

        sum_ = multiplicitySum( name, label, MT, summands=summands_, documentation=documentation,
                attributes = attributes )
        sum_.multiplicity.parseXMLNode( element.find( multiplicityModule.component.moniker ), xPath, linkData )

        xPath.pop()
        return sum_

class listOfSummands( ancestryModule.ancestry ):

    moniker = 'summands'

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
