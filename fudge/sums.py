# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the 'sum', 'summand' and 'Summands' classes
"""

from xData import link as linkModule
from xData import ancestry as ancestryModule
from xData import formatVersion as formatVersionModule

from fudge import warning as warningModule
from fudge import suites as suitesModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.channelData import Q as QModule
from fudge.productData import multiplicity as multiplicityModule

__metaclass__ = type

class sums( ancestryModule.ancestry ):
    """
    Contains all summed quantities. Currently supports summed cross sections and multiplicities,
    could extend to other types of sums later.
    """

    moniker = 'sums'
    ancestryMembers = ( 'crossSectionSums', 'multiplicitySums' )
    legacyMemberNameMapping = { 'crossSections': 'crossSectionSums', 'multiplicities': 'multiplicitySums' }

    def __init__(self):

        ancestryModule.ancestry.__init__( self )

        self.__crossSectionSums = suitesModule.crossSectionSums( )
        self.__crossSectionSums.setAncestor( self )

        self.__multiplicitySums = suitesModule.multiplicitySums( )
        self.__multiplicitySums.setAncestor( self )

    @property
    def crossSectionSums( self ) :

        return self.__crossSectionSums

    @property
    def multiplicitySums( self ) :

        return self.__multiplicitySums

    def convertUnits( self, unitMap ) :

        for child in self.ancestryMembers :
            getattr( self, child ).convertUnits( unitMap )

    def cullStyles( self, styleList ) :

        self.__crossSectionSums.cullStyles( styleList )
        self.__multiplicitySums.cullStyles( styleList )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers : item = getattr( self, ancestryMember ).findLinks( links )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__crossSectionSums.removeStyles( styleLabels )
        self.__multiplicitySums.removeStyles( styleLabels )

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
            elif linkData['format'] == formatVersionModule.version_1_10 and child.tag in self.legacyMemberNameMapping:
                getattr(self, self.legacyMemberNameMapping[child.tag]).parseXMLNode( child, xPath, linkData )
            else:
                print("Warning: unexpected child '%s' in sums" % child.tag)
        xPath.pop()

class crossSectionSum( ancestryModule.ancestry ):
    """
    Stores a summed quantity (cross section, multiplicity, etc.) along with the list of what's being summed over
    """

    moniker = 'crossSectionSum'
    ancestryMembers = ( 'summands', 'Q', 'crossSection' )

    def __init__( self, label, ENDF_MT ):
        """
        Create a 'sum' instance
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :return:
        """

        ancestryModule.ancestry.__init__( self )

        self.label = label
        self.ENDF_MT = int( ENDF_MT )

        self.__summands = Summands( )
        self.__summands.setAncestor( self )

        self.__Q = QModule.component( )
        self.__Q.setAncestor( self )

        self.__crossSection = crossSectionModule.component( )
        self.__crossSection.setAncestor( self )

    def __str__( self ):
        return self.label

    @property
    def summands( self ) :

        return( self.__summands )

    @property
    def crossSection( self ) :

        return( self.__crossSection )

    @property
    def Q( self ) :

        return( self.__Q )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__Q.convertUnits( unitMap )
        self.__crossSection.convertUnits( unitMap )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers : item = getattr( self, ancestryMember ).findLinks( links )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001, heatAllPoints = False,
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False, verbose = 0 ) :

        return( self.__crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
                heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero, addToSuite = True ) )

    def cullStyles( self, styleList ) :

        self.__Q.cullStyles( styleList )
        self.__crossSection.cullStyles( styleList )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__Q.removeStyles( styleLabels )
        self.__crossSection.removeStyles( styleLabels )

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ' ) :

        self.__crossSection.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        self.__crossSection.processMultiGroup( style, tempInfo, indent )

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

        sums = getPointwiseLinearForm( self.__summands[0].link )
        for idx in range( 1, len( self.__summands ) ) :
            pointwise = getPointwiseLinearForm( self.__summands[idx].link )
            sums, current = sums.mutualify( 1e-8, 1e-8, 0, pointwise, 1e-8, 1e-8, 0 )
            sums += current
        return sums

    def scaleSummandsToMatchSum(self, label=None):
        """
        rescale the cross sections in the summands list so that their sum matches the cross section in self
        :return:
        """

        if label is not None and label in self.__crossSection:
            correctXsec = self.__crossSection[ label ]
        else:
            correctXsec = self.__crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        summedXsec = self.sumSummands(label=label)
        if correctXsec.domain() != summedXsec.domain():
            raise warningModule.summedCrossSectionDomainMismatch( obj=self )
        scaleFactor=correctXsec/summedXsec

        for s in self.__summands:
            label = s.link.evaluated.label
            cutHere=s.link.evaluated.domainMin
            try:
                cutHere=max(cutHere,s.link.evaluated.ancestor.ancestor.getThreshold(unit=s.link.evaluated.domainUnit))
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
        if label is None: label = self.__crossSection.evaluated.label
        resummed.label = resummed.style = label
        if label == self.__crossSection.evaluated.label:
            self.__crossSection.remove( label )
        self.__crossSection.add( resummed )

    def check( self, info ) :

        warnings = self.__crossSection.check( info )

        # does self.__crossSection equal the sum over summand reaction cross sections?
        sums = self.sumSummands(label=info['reconstructedStyleName'])
        if info['reconstructedStyleName'] in self.__crossSection:
            quotedXsec = self.__crossSection[ info['reconstructedStyleName'] ]
        else:
            quotedXsec = self.__crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sums.domain() != quotedXsec.domain():
            warnings.append( warningModule.summedCrossSectionDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sums
            diff.nf_pointwiseXY.setSafeDivide(True)
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

        xmlString = [ '%s<%s label="%s" ENDF_MT="%s">' % ( indent, self.moniker, self.label, self.ENDF_MT ) ]
        xmlString += self.__summands.toXMLList( indent2, **kwargs )
        xmlString += self.__Q.toXMLList( indent2, **kwargs )
        xmlString += self.__crossSection.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )

        for child in element:
            if( child.tag not in ( Summands.moniker, QModule.component.moniker, crossSectionModule.component.moniker ) ) :
                raise TypeError( 'Unexpected element "%s" encountered' % child.tag )

        sums = crossSectionSum( element.get( 'label' ), int( element.get( 'ENDF_MT' ) ) )
        sums.__summands.parseXMLNode( element.find( Summands.moniker ), xPath, linkData )
        sums.__crossSection.parseXMLNode( element.find( crossSectionModule.component.moniker ), xPath, linkData )
        sums.__Q.parseXMLNode( element.find( QModule.component.moniker ), xPath, linkData )

        xPath.pop()
        return sums

class multiplicitySum( ancestryModule.ancestry ):
    """
    Stores a summed multiplicity along with the list of what's being summed over
    """

    moniker = 'multiplicitySum'
    ancestryMembers = ( 'summands', 'multiplicity' )

    def __init__( self, label, ENDF_MT ):
        """
        Create a 'sum' instance
        :param name: descriptive string
        :param label: unique string identifying this sum
        :param ENDF_MT: for converting back to ENDF
        :return:
        """

        ancestryModule.ancestry.__init__( self )
        self.label = label

        self.ENDF_MT = int( ENDF_MT )

        self.__summands = Summands( )
        self.__summands.setAncestor( self )

        self.__multiplicity = multiplicityModule.component( )
        self.__multiplicity.setAncestor( self )

    def __str__( self ):
        return self.label

    @property
    def summands( self ) :

        return( self.__summands )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__multiplicity.convertUnits( unitMap )

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers : item = getattr( self, ancestryMember ).findLinks( links )

    def check( self, info ) :

        warnings = []

        if self.__multiplicity.isConstant() and self.__multiplicity.getConstant() < 1:
            warnings.append( warningModule.negativeMultiplicity( self.getConstant(), self ) )
        else:   # energy-dependent mult.
            for form in self.__multiplicity :
                if hasattr(form, 'rangeMin') and form.rangeMin < 0:
                    warnings.append( warningModule.negativeMultiplicity( form.rangeMin, obj=form ) )

        # does multiplicity equal the sum over its summand multiplicities?
        sums = self.__summands[0].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        for idx in range( 1, len( self.__summands ) ) :
            sums, current = sums.mutualify( 1e-8, 1e-8, 0,
                    self.__summands[idx].link.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8), 1e-8, 1e-8, 0 )
            sums += current
        quotedXsec = self.__multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
        if sums.domain() != quotedXsec.domain():
            warnings.append( warningModule.summedMultiplicityDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sums
            diff.nf_pointwiseXY.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            try:    # FIXME instead of try/catch, should really compare relativeDiff.rangeMax with 'inf'
                if relativeDiff.rangeMax > info['multiplicityMaxDiff']:
                    warnings.append( warningModule.summedMultiplicityMismatch( relativeDiff.rangeMax * 100, obj=self ) )
            except TypeError:
                warnings.append( warningModule.summedMultiplicityZeroDivision( obj=self ) )

        return warnings

    def getLabel(self):

        return self.label

    def cullStyles( self, styleList ) :

        self.__multiplicity.cullStyles( styleList )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__multiplicity.removeStyles( styleLabels )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ""
        xmlString = [ '%s<%s label="%s" ENDF_MT="%s"%s>' % ( indent, self.moniker, self.label, self.ENDF_MT, attrStr ) ]
        xmlString += self.__summands.toXMLList( indent2, **kwargs )
        xmlString += self.__multiplicity.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )

        for key in element.keys( ) :
            if( key == 'ENDF_MT' ) :
                MT = int( element.get( key ) )
            elif( key == 'label' ) :
                label = element.get( key )
            else:
                raise ValueError( 'Unsupported attribute "%s"' % key )

        summands = None
        for child in element:
            if( child.tag not in ( Summands.moniker, multiplicityModule.component.moniker ) ) :
                raise TypeError( 'Unexpected element "%s" encountered' % child.tag )

        sums = multiplicitySum( label, MT )
        sums.__summands.parseXMLNode( element.find( Summands.moniker ), xPath, linkData )
        sums.__multiplicity.parseXMLNode( element.find( multiplicityModule.component.moniker ), xPath, linkData )

        xPath.pop()
        return sums

class Summands( ancestryModule.ancestry ):

    moniker = 'summands'
    ancestryMembers = ( '', )

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )
        self.__summands = []

    def __getitem__(self, index):

        return self.__summands[index]

    def __len__(self):

        return len( self.__summands )

    def append( self, item ) :

#        if( not( isinstance( iter, ( add, subtract ) ) ) ) : raise TypeError( 'Invalid item to add: type = %s.' % type( item ) )

        self.__summands.append( item )
        item.setAncestor( self )

    def findLinks( self, links ) :

        for item in self :
            if( isinstance( item, ( linkModule.link, linkModule.link2 ) ) ) : links.append( [ item, item.link, item.path ] )

    def index( self, item ) :
        """Returns the index of item in self."""

        return( self.__summands.index( item ) )

    def pop( self, index = -1 ) :
        """Removes and return item at index. Default is the last item."""

        return( self.__summands.pop( index ) )

    def toXMLList( self, indent = '', **kwargs ) :

        xmlList = ['%s<%s>' % (indent, self.moniker)]
        for summand in self.__summands :
            xmlList.append( '  %s%s' % (indent, summand.toXML() ) )
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if child.tag == add.moniker :
                self.append( add.parseXMLNode( child, xPath, linkData ) )
            elif child.tag == subtract.moniker :
                self.append( subtract.parseXMLNode( child, xPath, linkData ) )

        xPath.pop( )

class add( linkModule.link ):
    """Link representing one of the quantities that is added to the sum."""

    moniker = 'add'

class subtract( linkModule.link ):
    """Link representing one of the quantities that is subtracted from the sum."""

    moniker = 'subtract'
