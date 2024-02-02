# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""        
This module contains classes and functions supporting coherent and incoherent photon scattering distribution.
                
    This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | BaseForm                          | Base class for :py:class:`CoherentPhotonScattering.From` and          |
    |                                   | :py:class:`IncoherentPhotonScattering.From` classes.                  |
    +-----------------------------------+-----------------------------------------------------------------------+
    | CoherentPhotonScattering          | Filler class for storing its Form class.                              |
    +-----------------------------------+-----------------------------------------------------------------------+
    | CoherentPhotonScattering.From     | A distribution form class with a link to an coherent photon           |
    |                                   | scattering instance in the double differential cross section node.    |
    +-----------------------------------+-----------------------------------------------------------------------+
    | IncoherentPhotonScattering        | Filler class for storing its Form class.                              |
    +-----------------------------------+-----------------------------------------------------------------------+
    | IncoherentPhotonScattering.From   | A distribution form class with a link to an incoherent photon         |
    |                                   | scattering instance in the double differential cross section node.    |
    +-----------------------------------+-----------------------------------------------------------------------+
"""     

from xData import enums as xDataEnumsModule
from xData import link as linkModule

from . import base as baseModule

class BaseForm( baseModule.Form, linkModule.Link ) :
    """
    Base class for :py:class:`CoherentPhotonScattering.From` and :py:class:`IncoherentPhotonScattering.From` classes.

    :param link:        The linked instance.
    :param root:        ?.
    :param path:        The xLink string to the linked instance.
    :param relative:    If True, link is relative; otherwise it is an absolute link.
    :param label:    
    """

    def __init__( self, link = None, root = None, path = None, relative = False, label = None ) :

        linkModule.Link.__init__( self, link = link, root = root, path = path, relative = relative, label = label )
        baseModule.Form.__init__( self, label, xDataEnumsModule.Frame.none, [] )

    @property
    def productFrame( self ) :
        """Returns the frame the product data are specified in."""

        return( self.link.productFrame )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        return( [] )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        return( self.link.calculateAverageProductData( style, indent, **kwargs ) )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """This method does nothing."""

        return 0

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an :py:class:`energyAngularMCModule.Form` instance representing *self*.

    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        return( self.link.processMC_cdf( style, tempInfo, indent ) )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        return( self.link.processMultiGroup( style, tempInfo, indent ) )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        return( linkModule.Link.toXML_strList( self, indent = indent, **kwargs ) )

class CoherentPhotonScattering :

    class Form(BaseForm):

        moniker = 'coherentPhotonScattering'

class IncoherentPhotonScattering :

    class Form(BaseForm):

        moniker = 'incoherentPhotonScattering'
