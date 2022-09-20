# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the standard uncertainty class.
"""

from . import uncertainty as uncertaintyModule

class Standard( uncertaintyModule.Base ) :

    moniker = "standard"

    def __init__( self, value ) :

        uncertaintyModule.Base.__init__( self )

        if not isinstance(value, uncertaintyModule.Quantity): raise TypeError( 'Invalid quantity for %s uncertainty' % self.moniker )
        self.__value = value
        self.__value.setAncestor( self )

    @property
    def value( self ) :

        return( self.__value )

    def copy( self ) :

        return( self.__class__( self.__value.copy( ) ) )

    def parentConvertingUnits( self, factors ) :

        self.__value.parentConvertingUnits( factors )

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XML_strList += self.value.toXML_strList(indent = indent2, **kwargs)
        XML_strList[-1] += "</%s>" % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        child = node[0]
        if child.tag == uncertaintyModule.Double.moniker:
            uncertainty = uncertaintyModule.Double.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else :
            raise Exception('Invalid child node = "%s".' % child.tag)

        instance = Standard(uncertainty)

        xPath.pop()

        return instance
