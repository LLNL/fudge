# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import enums as xDataEnumsModule
from xData import physicalQuantity as physicalQuantityModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as IDsPoPsModule

from . import base as baseModule
from . import incoherentElasticMisc as incoherentElasticMiscModule

"""
Thermal neutron scattering law incoherent elastic double difference cross section form and its supporting classes.
"""

class BoundAtomCrossSection(physicalQuantityModule.PhysicalQuantity):

    moniker = 'boundAtomCrossSection'

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        moniker = self.moniker
        if formatVersion == GNDS_formatVersionModule.version_1_10:
            moniker = 'characteristicCrossSection'
        kwargs['moniker'] = moniker

        return physicalQuantityModule.PhysicalQuantity.toXML_strList(self, indent, **kwargs)

class DebyeWallerIntegral( ancestryModule.AncestryIO ) :
    """
    For incoherent elastic sections, all we need is a characteristic cross section and a
    temperature-dependent list of the Debye-Waller integral.
    """

    moniker = "DebyeWallerIntegral"

    def __init__( self, function1d ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__function1d = function1d
        self.__function1d.setAncestor( self )

    @property
    def function1d( self ) :

        return( self.__function1d )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.__function1d.convertUnits( unitMap )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        moniker = self.moniker

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        if formatVersion == GNDS_formatVersionModule.version_1_10:
            moniker = 'DebyeWaller'

        xmlStringList = ['%s<%s>' % (indent, moniker)]
        xmlStringList += self.__function1d.toXML_strList(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % moniker

        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param node:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *node*.
        """

        xPath.append( node.tag )

        xys1d = XYs1dModule.XYs1d.parseNodeUsingClass(node[0], xPath, linkData, **kwargs)
        debyeWallerIntegral = cls(xys1d)

        xPath.pop( )
        return debyeWallerIntegral

class Form( baseModule.Form ) :

    moniker = 'thermalNeutronScatteringLaw_incoherentElastic'
    keyName = 'label'

    process = 'thermalNeutronScatteringLaw incoherent-elastic'
    subformAttributes = ( 'boundAtomCrossSection', 'DebyeWallerIntegral' )

    def __init__( self, label, boundAtomCrossSection, debyeWallerIntegral ) :

        baseModule.Form.__init__( self, IDsPoPsModule.neutron, label, xDataEnumsModule.Frame.lab, ( boundAtomCrossSection, debyeWallerIntegral, ) )

    @property
    def domainUnit( self ) :
        """
        This method rReturns the energy unit of the projectile.
        """

        invUnit = 1 / PQUModule.PQU( 1, self.DebyeWallerIntegral.function1d.axes[0].unit )
        return( str( invUnit.unit ) )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        temperature = style.temperature.getValueAs( self.DebyeWallerIntegral.function1d.axes[1].unit )

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )
        energyMax = PQUModule.PQU(  5e-6, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMax = kwargs.get( 'energyMax', energyMax )

        return( incoherentElasticMiscModule.process( self, style.label, energyMin, energyMax, temperature, kwargs ) )

    def temperatures( self, _temperatures ) :

        _temperatures['incoherent-elastic'] = [ self.DebyeWallerIntegral.function1d.axes[-1].unit, [ x for x, y in self.DebyeWallerIntegral.function1d ] ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param node:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        boundAtomCrossSection = BoundAtomCrossSection.parseNodeUsingClass(node[0], xPath, linkData, **kwargs)
        debyeWallerIntegral = DebyeWallerIntegral.parseNodeUsingClass(node[1], xPath, linkData, **kwargs)
        incoherentElastic = cls(node.get('label'), boundAtomCrossSection, debyeWallerIntegral)

        xPath.pop( )
        return( incoherentElastic )
