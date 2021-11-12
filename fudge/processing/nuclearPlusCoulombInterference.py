# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
Defines the NuclearPlusCoulombInterference class which is used to store the elastic scattering reaction for a protare
with a charged particle as the projectile where only the 'nuclear + interference' data are included.  Ergo, the 
Rutherford scattering term is ignored. This reaction is equivalent to the ENDL C=9 reaction.
"""

from xData import ancestry as ancestryModule

from .. import outputChannel as outputChannelModule
from ..reactions import reaction as reactionModule

__metaclass__ = type

class NuclearPlusCoulombInterference( ancestryModule.Ancestry2 ) :
    """
    This class has only one member which is the 'nuclear + interference' reaction.
    """

    moniker = 'nuclearPlusCoulombInterference'
    ancestryMembers = [ 'reaction' ]

    def __init__( self ) :

        ancestryModule.Ancestry2.__init__( self )

        self.__reaction = reactionModule.reaction( outputChannelModule.Genre.twoBody, 2 )
        self.__reaction.setAncestor( self )

    @property
    def reaction( self ) :

        return( self.__reaction )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLList += self.__reaction.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        return( NuclearPlusCoulombInterference( ) )
