# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes for storing a multi-group distribution.

This module contains the following classes:

    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | Subform       | Abstract base class for multiGroup subforms.                          |
    +---------------+-----------------------------------------------------------------------+
    | Gridded3d     | Gridded 3d data representing the multi-group transfer matrices.       |
    +---------------+-----------------------------------------------------------------------+
    | Form          | Class representing GNDS multi-group distribution data.                |
    +---------------+-----------------------------------------------------------------------+
"""

import abc

from xData import gridded as griddedModule

from . import base as baseModule

class Subform( baseModule.Subform, abc.ABC ) :
    """
    Abstract base class for multiGroup subforms.
    """

    def __init__( self ) :

        baseModule.Subform.__init__( self )

class Gridded3d( Subform, griddedModule.Gridded3d ) :
    """
    The class for storing the multi-group distribution data (i.e., transfer matrices).

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | axes      | Unique label for the form within the distribution instance.   |
    +-----------+---------------------------------------------------------------+
    | array     | The frame the product data are in.                            |
    +-----------+---------------------------------------------------------------+
    | kwargs    | Instance storing the multi-group data.                        |
    +-----------+---------------------------------------------------------------+

    :param axes:        Axes for the data.
    :param array:       Multi-group data as an array.
    :param kwargs:      A dictionary that contains data to control the way this mmethod acts.
    """

    def __init__(self, axes, array, **kwargs):

        Subform.__init__(self)
        griddedModule.Gridded3d.__init__(self, axes, array, **kwargs)

    @staticmethod
    def dataToString( values, self, indent = '', **kwargs ) :

        valueFormatter = kwargs['valueFormatter']

        XMLStrList = [ ]
        start = 0
        for length in self.lengths :
            strList = [ valueFormatter( self.data[i1] ) for i1 in range( start, start + length ) ]
            XMLStrList.append( indent + ' '.join( strList ) )
            start += length
        return( XMLStrList )

class Form( baseModule.Form ) :
    """
    Class for storing the distribution that is the multi-group data needed by deterministic transport codes.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | label             | Unique label for the form within the distribution instance.   |
    +-------------------+---------------------------------------------------------------+
    | productFrame      | The frame the product data are in.                            |
    +-------------------+---------------------------------------------------------------+
    | multiGroupSubform | Instance storing the multi-group data.                        |
    +-------------------+---------------------------------------------------------------+

    :param label:           Unique label for the form within the distribution instance.
    :param productFrame:    The frame the product data are in.
    :param _energy:         The energy child node containing P(E'|E).
    :param _energyAngular:  The angular child node containing P(mu|E,E').
    """

    moniker = 'multiGroup3d'
    subformAttributes = ( 'multiGroupSubform', )

    def __init__( self, label, productFrame, multiGroupSubform ) :

        if( not( isinstance( multiGroupSubform, Subform ) ) ) :
            raise TypeError( 'instance is not a multiGroup subform: moniker = %s' % multiGroupSubform.moniker )
        baseModule.Form.__init__( self, label, productFrame, ( multiGroupSubform, ) )

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

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {
                Gridded3d.moniker : Gridded3d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Gridded3d subform: %s" % subformElement.tag )
        SubForm = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        instance = cls( element.get( 'label' ), element.get( 'productFrame' ), SubForm )

        xPath.pop( )

        return  instance
