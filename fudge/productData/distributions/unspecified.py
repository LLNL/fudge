# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains class and functions needed for storing and processing a GNDS *unspecified* distribution. It
contains the following classes:
            
    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | Form          | Class representing GNDS unspecified distribution.                     |
    +---------------+-----------------------------------------------------------------------+
"""

from xData import enums as xDataEnumsModule
from . import base as baseModule

class Form( baseModule.Form ) :
    """
    Class representing a GNDS unspecified distribution form.

    :param label:                   The label for this form.
    :param productFrame:            The frame the product data are specified in.
    """

    moniker = 'unspecified'
    subformAttributes = []

    def __init__(self, label, productFrame=xDataEnumsModule.Frame.lab):

        if productFrame is None:
            productFrame = xDataEnumsModule.Frame.lab      # This is a kludge and needs to be fixed. Ergo, None should never be allowed as a productFrame value.
        baseModule.Form.__init__( self, label, productFrame, [] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """     
        Parses *element* into an instance *cls*.

        :param cls:                 Form class to return.
        :param element:             Node to parse.
        :param xPath:               List containing xPath to current node, useful mostly for debugging.
        :param linkData:            dict that collects unresolved links.    
        :param kwargs:              A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
                        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        UC = cls( element.get( 'label' ), element.get( 'productFrame' ) )
        xPath.pop( )
        return( UC )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        return( self.rootAncestor.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        return( [ None, None ] )

    def copy( self ):
        """
        Returns a copy of *self*.

        :return:            :py:class`Form` instance.
        """

        return Form( self.label, self.productFrame )

    def fixDomains(self, energyMin, energyMax, domainsToFix):
        """
        The method does nothing and returns 0.

        :return:                0.
        """

        return 0

    def isSpecified( self ) :
        """
        Since this is the GNDS unspecified distribution, this method returns False.

        :return:                False.
        """

        return( False )

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods does nothing and returns None.
    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                None.
        """     

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        This methods does nothing and returns None.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                None.
        """

        return( None )
