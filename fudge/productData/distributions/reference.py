# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains distribution forms that a references to other distribution forms.

This module contains the following classes:
        
    +-------------------------------+-----------------------------------------------------------------------+
    | Class                         | Description                                                           |
    +===============================+=======================================================================+
    | Form                          | Class for representing GNDS a reference distribution form.            |
    +-------------------------------+-----------------------------------------------------------------------+
    | CoulombPlusNuclearElastic     | Class for :py:class:`CoulombPlusNuclearElastic`.                      |
    +-------------------------------+-----------------------------------------------------------------------+
    | ThermalNeutronScatteringLaw   | Class for :py:class:`ThermalNeutronScatteringLaw`.                    |
    +-------------------------------+-----------------------------------------------------------------------+
""" 

from xData import enums as xDataEnumsModule
from xData import link as linkModule

from . import base as baseModule

class Form( linkModule.Link, baseModule.Form ) :
    """
    This class is used to link to another distribution form.

    :param link:            The label for this form.
    :param root:            External file identifier. Required if the link points outside the current file.
    :param path:            A string representing the xPath to *link*.
    :param label:           The label for this form.
    :param relative:        If True use relative link when writing out xPath otherwise use absolute xPath.
    """

    moniker = 'reference'
    subformAttributes = []

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.Link.__init__( self, link = link, root = root, path = path, label = label, relative = relative )
        baseModule.Form.__init__( self, label, xDataEnumsModule.Frame.none, [] )

    @property
    def referenceInstance( self ):
        """
        Returns a reference to the instance referred to by *self*.
        """

        if self.link is None: raise Exception("Unresolved link!")
        return self.link

    @property
    def productFrame( self ):
        """
        Returns the product frame of the instance referred to by *self*.
        """

        return self.referenceInstance.productFrame

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        return( self.referenceInstance.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        pass

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in. 
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        
        :return:                    XYs1d instance for the energy spectrum.
        """

        return(self.referenceInstance.energySpectrumAtEnergy(energyIn, frame, **kwargs))

    def fixDomains(self, labels, energyMin, energyMax):
        """This method does nothing."""

        return 0

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an the results of calling **processMC_cdf** on the linked instance.

    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        # temporary solution:
        tempInfo['viaReference'] = True
        instance = self.referenceInstance.processMC_cdf(style, tempInfo, indent)
        tempInfo['viaReference'] = False

        return instance

        # better solution: add another link pointing to the processed version of what this points to:
        """
        newReference = Form( label=style.label, relative=True )
        self.ancestor.add( newReference )

        tempInfo['brokenLinks'].append( [self, newReference] )
        """

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        return( self.referenceInstance.processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of instance referenced by *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    Instance returned by *toPointwise_withLinearXYs* method the referenced instance.
        """

        return( self.referenceInstance.toPointwise_withLinearXYs( **kwargs ) )

class CoulombPlusNuclearElastic(Form):
    """
    An instance of this class links to a CoulombPlusNuclearElastic instance in the double differential cross section node.
    """

    moniker = 'CoulombPlusNuclearElastic'

class ThermalNeutronScatteringLaw( Form ) :
    """
    An instance of This class links to a ThermalNeutronScatteringLaw instance in the double differential cross section node.
    """

    moniker = 'thermalNeutronScatteringLaw'
