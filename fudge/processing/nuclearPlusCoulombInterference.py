# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
Defines the NuclearPlusCoulombInterference class which is used to store the elastic scattering reaction for a protare
with a charged particle as the projectile where only the 'nuclear + interference' data are included.  Ergo, the 
Rutherford scattering term is ignored. This reaction is equivalent to the ENDL C=9 reaction.
"""

from LUPY import ancestry as ancestryModule

from .. import enums as enumsModule
from .. import outputChannel as outputChannelModule
from ..reactions import reaction as reactionModule


class NuclearPlusCoulombInterference( ancestryModule.AncestryIO ) :
    """
    This class is designed to store an LLNL ENDL C=9 reaction which the elatic scattering between two charged particle
    (e.g., 'p + H2') but without the Rutherford scattering term. This reaction is designed to support a legacy LLNL
    type reaction and should not be used in general since its cross section can be negative.
    The reaction in this class is not supported by ENDF of GNDS and is only put in the application node by FUDGE.
    In ENDL, this reaction is called a 'nuclear + interference' reaction.

    The following table list the primary members of this class:
    
    +---------------+-----------------------------------------------------------+
    | Member        | Description                                               |
    +===============+===========================================================+
    | reaction      | The ENDL C=9 reaction.                                    |
    +---------------+-----------------------------------------------------------+
    """

    moniker = 'nuclearPlusCoulombInterference'
    ancestryMembers = ( 'reaction', )

    def __init__(self, label):
        """
        :param label:       The label for *self* within the LLNL institution.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__reaction = reactionModule.Reaction(label, enumsModule.Genre.twoBody, 2)
        self.__reaction.setAncestor( self )

    @property
    def reaction( self ) :
        """This function returns a reference to the reaction."""

        return( self.__reaction )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a python list of str instances representing the XML lines of *self*.
        
        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                Python list of str instances.
        """ 

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLList += self.__reaction.toXML_strList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """ 
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Dictionary that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        instance = cls(node[0].get('label'))
        instance.reaction.parseNode(node[0], xPath, linkData, **kwargs)

        return instance
