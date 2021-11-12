# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from PoPs import IDs as IDsPoPsModule
from PoPs.fissionFragmentData import rate as rateModule

from fudge import suites as suitesModule
from fudge import product as productModule

class delayedNeutron( ancestryModule.ancestry ) :

    moniker = 'delayedNeutron'

    def __init__( self, label, product ) :

        ancestryModule.ancestry.__init__( self )

        self.__label = label

        self.__rate = rateModule.suite( )
        self.__rate.setAncestor( self )

        self.__product = product
        self.__product.setAncestor( self )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def rate( self ) :

        return( self.__rate )

    @property
    def product( self ) :

        return( self.__product )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        self.__product.calculateAverageProductData( style, indent = indent, **kwargs )

    def check( self, info ):
        """
        Check delayed neutron rate, multiplicity and distribution

        :param info: dict
        :return: list of warnings
        """
        from fudge import warning
        warnings = []

        if self.rate[0].value < 0:
            warnings.append( warning.negativeDelayedRate(self.rate.time, self) )

        warnings += self.product.check( info )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__rate.convertUnits( unitMap )
        self.__product.convertUnits( unitMap )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        self.__product.processMC_cdf( style, tempInfo, indent = indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productName'] = self.__product.id
        tempInfo['productLabel'] = self.__product.label
        self.__product.processMultiGroup( style, tempInfo, indent )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

#        self.__rate.removeStyles( styleLabels )
        self.__product.removeStyles( styleLabels )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        XMLStringList += self.__rate.toXMLList( indent2, **kwargs )
        if( self.__product is not None ) : XMLStringList += self.__product.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 2 ) : return( [] )
        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        product = None
        for child in element :
            if( child.tag == productModule.product.moniker ) :
                product = productModule.product.parseXMLNode( child, xPath, linkData )
                break

        delayedNeutron1 = cls( element.get( 'label' ), product )

        for child in element :
            if( child.tag == delayedNeutron1.__rate.moniker ) :
                delayedNeutron1.__rate.parseXMLNode( child, xPath, linkData )
            elif( child.tag == productModule.product.moniker ) :
                pass
            else :
                raise TypeError( "Encountered unknown node '%s' in %s" % ( child.tag, element.tag ) )

        xPath.pop( )

        return( delayedNeutron1 )

class delayedNeutrons( suitesModule.suite ) :

    moniker = 'delayedNeutrons'

    def __init__( self ) :

        suitesModule.suite.__init__( self, allowedClasses = ( delayedNeutron, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        for delayedNeutron1 in self : delayedNeutron1.calculateAverageProductData( style, indent = indent, **kwargs )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        for delayedNeutron1 in self : delayedNeutron1.processMC_cdf( style, tempInfo, indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        for delayedNeutron1 in self : delayedNeutron1.processMultiGroup( style, tempInfo, indent )
