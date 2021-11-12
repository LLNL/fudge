# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from xData import ancestry as ancestryModule

from fudge.channelData.fissionFragmentData import fissionEnergyRelease as fissionEnergyReleaseModule

from . import delayedNeutron as delayedNeutronModule
from . import productYield as productYieldModule

class fissionFragmentData( ancestryModule.ancestry ) :

    moniker = 'fissionFragmentData'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__delayedNeutrons = delayedNeutronModule.delayedNeutrons( )
        self.__delayedNeutrons.setAncestor( self )

        self.__fissionEnergyReleases = fissionEnergyReleaseModule.fissionEnergyReleases( )
        self.__fissionEnergyReleases.setAncestor( self )

        self.__productYields = productYieldModule.suite( )
        self.__productYields.setAncestor( self )
        
    @property
    def delayedNeutrons( self ) :

        return( self.__delayedNeutrons )

    @property
    def fissionEnergyReleases( self ) :

        return( self.__fissionEnergyReleases )

    @property
    def productYields( self ) :

        return( self.__productYields )

    def amendForPatch( self, fromLabel, toLabel ) :

        self.__delayedNeutrons.amendForPatch( fromLabel, toLabel )
        self.__fissionEnergyReleases.amendForPatch( fromLabel, toLabel )
        self.__productYields.amendForPatch( fromLabel, toLabel )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        self.__delayedNeutrons.calculateAverageProductData( style, indent = indent, **kwargs )

    def check( self, info ):
        """
        Check for problems in delayed fission products, product yields and energy release data

        :param info:
        :return:
        """
        from fudge import warning
        warnings = []

        for term in ('productYields', 'delayedNeutrons', 'fissionEnergyReleases'):
            for data in getattr(self, term):
                dwarnings = data.check( info )
                if dwarnings:
                    warnings.append( warning.context('%s: %s' % (term, data.label), dwarnings) )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.__delayedNeutrons.convertUnits( unitMap )
        self.__fissionEnergyReleases.convertUnits( unitMap )
        self.__productYields.convertUnits( unitMap )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        self.__delayedNeutrons.processMC_cdf( style, tempInfo, indent = indent )

    def processMultiGroup( self, style, tempInfo, indent ) :

        self.__delayedNeutrons.processMultiGroup( style, tempInfo, indent )
        self.__fissionEnergyReleases.processMultiGroup( style, tempInfo, indent )

    def replicate( self, other ) :

        self.__delayedNeutrons = other.delayedNeutrons
        self.__delayedNeutrons.setAncestor( self )

        self.__fissionEnergyReleases = other.fissionEnergyReleases
        self.__fissionEnergyReleases.setAncestor( self )

        self.__productYields = other.productYields
        self.__productYields.setAncestor( self )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__delayedNeutrons.removeStyles( styleLabels )
        self.__fissionEnergyReleases.removeStyles( styleLabels )
        self.__productYields.removeStyles( styleLabels )

    def saveToOpenedFile( self, fOut, **kwargs ) :

        xmlString = self.toXMLList( **kwargs )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, **kwargs ) :
        """
        Save the fissionFragmentData in GNDS/xml structure to specified file.
        To suppress extra messages, change the 'verbosity' flag:
        >>>self.saveToFile( "output.xml", flags = { 'verbosity' : 0 } )
        """

        dirname = os.path.dirname( fileName )
        if( ( len( dirname ) > 0 ) and not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        with open( fileName, "w" ) as fout :
            self.saveToOpenedFile( fout, **kwargs )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__delayedNeutrons.toXMLList( indent2, **kwargs )
        XMLStringList += self.__fissionEnergyReleases.toXMLList( indent2, **kwargs )
        XMLStringList += self.__productYields.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : XMLStringList = []

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if( child.tag == self.__delayedNeutrons.moniker ) :
                self.__delayedNeutrons.parseXMLNode( child, xPath, linkData )
            elif( child.tag == self.__fissionEnergyReleases.moniker ) :
                self.__fissionEnergyReleases.parseXMLNode( child, xPath, linkData )
            elif( child.tag == self.__productYields.moniker ) :
                self.__productYields.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError( "Encountered unknown node '%s' in %s" % ( child.tag, element.tag ) )

        xPath.pop( )
