# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from . import multiplicity as multiplicityModule
from . import decayConstant as decayConstantModule

# FIXME CMM: class appears to be unused. Delete?

class DelayedNeutronData(ancestryModule.AncestryIO):

    # FIXME does not agree with GNDS-2.0 specs. Should just be 'delayedNeutron'?
    moniker = 'delayedNeutronData'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__multiplicity = multiplicityModule.Suite( )
        self.__multiplicity.setAncestor( self )

        self.__decayConstants = decayConstantModule.Suite( )
        self.__decayConstants.setAncestor( self )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    @property
    def decayConstants( self ) :

        return( self.__decayConstants )

    def convertUnits( self, unitMap ) :

        self.__multiplicity.convertUnits( unitMap )
        self.__decayConstants.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__multiplicity.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__decayConstants.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : return( [] )
        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        for child in element:
            if child.tag == multiplicityModule.Suite.moniker:
                self.multiplicity.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == decayConstantModule.Suite.moniker:
                self.decayConstants.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return self

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return self
