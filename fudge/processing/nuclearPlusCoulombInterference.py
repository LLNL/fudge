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
    This class has only one member which is the 'nuclear + interference' reaction.
    """

    moniker = 'nuclearPlusCoulombInterference'
    ancestryMembers = ( 'reaction', )

    def __init__(self, label):

        ancestryModule.AncestryIO.__init__(self)

        self.__reaction = reactionModule.Reaction(label, enumsModule.Genre.twoBody, 2)
        self.__reaction.setAncestor( self )

    @property
    def reaction( self ) :

        return( self.__reaction )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLList += self.__reaction.toXML_strList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        instance = cls(node[0].get('label'))
        instance.reaction.parseNode(node[0], xPath, linkData, **kwargs)

        return instance
