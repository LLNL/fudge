# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""This module contains the referred class."""

from fudge.core.ancestry import ancestry

__metaclass__ = type

monikerReferredData = 'referredData'

def parseXMLNode( element, linkData={} ):
    """ translate <referredData> element from xml. Currently only crossSection elements allowed """
    import fudge
    referred = referredData()
    for dat in element:
        xsc = fudge.gnd.reactionData.crossSection.parseXMLNode( dat[0], linkData )
        referred.appendDatum( xsc )
    return referred

class referredData( ancestry ) :

    def __init__( self, parent = None ) :

        ancestry.__init__( self, monikerReferredData, parent )
        self.data = []

    def __len__( self ) :

        return( len( self.data ) )

    def __getitem__( self, index ) :

        return( self.data[index] )

    def appendDatum( self, data ) :

        label = str( len(self.data) )
        self.data.append( referredDatum( self, label, data ) )
        data.setParent( self.data[-1] )

    def toXMLList( self, indent = "" ) :

        if( len( self ) == 0 ) : return( [] )
        indent2, xmlString = indent + '  ', []
        xmlString.append( '%s<%s>' % ( indent, self.moniker ) )
        for datum in self : xmlString += datum.toXMLList( indent2 )
        xmlString[-1] += '</%s>' % ( self.moniker )
        return( xmlString )

class referredDatum( ancestry ) :

    def __init__( self, parent, label, data ) :

        ancestry.__init__( self, 'key', parent, attribute = "label" )
        self.label = label
        self.data = data
        data.setParent( self )

    def toXMLList( self, indent = "" ) :

        xmlString = [ '%s<key label="%s">' % ( indent, self.label ) ]
        xmlString += self.data.toXMLList( indent + '  ' )
        xmlString[-1] += "</%s>" % self.moniker
        return( xmlString )
