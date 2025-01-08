# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Thd module contain miscellaneous classes and function for charged particle elastive scattering.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | ChargedParticleElasticTerm        | This is a base class for the various terms (e.g., real or imaginary   |
    |                                   | interfernce term).                                                    |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | fixMuRange                        | This function ensures that the range of mu in *data* is [-1, 1].      |
    |                                   | bu added point if needed.
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import abc
import math

from LUPY import ancestry as ancestryModule

from xData import axes as axesModule


def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :
    """
    This function returns an :py:class:`axesModule.Axes` instance for a double difference cross section.

    :param energyUnit:          The unit for energy.
    :oaram crossSectionUnit:    The unit for cross section.
    """

    axes = axesModule.Axes(3)
    axes[0] = axesModule.Axis( 'dSigma_dMu', 0, crossSectionUnit )
    axes[1] = axesModule.Axis( 'mu', 1, '' )
    axes[2] = axesModule.Axis( 'energy_in', 2, energyUnit )

    return( axes )

class ChargedParticleElasticTerm(ancestryModule.AncestryIO, abc.ABC):
    """
    This is a base class for the various terms (e.g., real or imaginary interfernce term).
    """

    @property
    @abc.abstractmethod
    def allowedDataForms(self):
        """Abstract class that must be overwritten by the derived class."""

        pass

    def __init__(self, data):
        """
        :param data:    The data for the term.
        """

        ancestryModule.AncestryIO.__init__(self)
        self.data = data

    @property
    def data(self):
        """This method returns a reference to the data."""

        return self.__data

    @data.setter
    def data(self, value):
        """ 
        This method sets the data to *value*.
    
        :param value:       The nuclear term.
        """

        if not isinstance(value, self.allowedDataForms):
            raise TypeError("'%s' cannot contain data of type '%s'" % (self.moniker, type(value)))
        value.setAncestor( self )
        self.__data = value

    def convertUnits( self, unitMap ):
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXML_strList(indent+'  ', **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append(element.tag)

        child = element[0]
        parseClass = None
        for allowedForm in cls.allowedDataForms:
            if child.tag == allowedForm.moniker:
                parseClass = allowedForm
                break
        if parseClass is None:
            raise TypeError("Unexpected element '%s' encountered inside %s" % (child.tag, cls.moniker))
        data = parseClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        term_ = cls( data )

        xPath.pop()

        return term_

def fixMuRange(data, epsilon=1e-6):
    """
    This function ensures that the range of mu in *data* is [-1, 1] by adding 0.0's to the missing ranges. Note, *data* is modified.

    :param data:        A python list of :math:`P(\\mu)`.
    :param epsilon:     The relative amount to add a point below (above) the first (last) point in addition to the point at -1 (1).
    """

    m1 = data[0][0]
    if m1 > -1.0:
        if abs(m1) == 0.0:
            m2 = -epsilon
        else:
            eps = math.copysign(epsilon, -m1)
            m2 = (1 + eps) * m1
            if m2 <= -1.0:
                data[0][0] = (1 - eps) * m1
                m2 = m1
        data.insert(0, [m2, 0.0])
        data.insert(0, [-1.0, 0.0])

    m1 = data[-1][0]
    if m1 < 1.0:
        if abs(m1) == 0.0:
            m2 = epsilon
        else:
            eps = math.copysign(epsilon, m1)
            m2 = (1 + eps) * m1
            if m2 >= 1.0:
                data[-1][0] = (1 - eps) * m1
                m2 = m1
        data.append([m2, 0.0])
        data.append([1.0, 0.0])
