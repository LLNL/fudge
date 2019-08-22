# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the decay classes.
"""

import abc

from xData import ancestry as ancestryModule

from .. import suite as suiteModule

from . import probability as probabilityModule
from . import product as productModule
from . import Q as QModule

class decay( ancestryModule.ancestry ) :

    moniker = 'decay'

    def __init__( self, label ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label not str' )
        self.__label = label

        self.__probability = probabilityModule.suite( )
        self.__probability.setAncestor( self )

        self.__Q = QModule.suite( )
        self.__Q.setAncestor( self )

        self.__products = productModule.suite( )
        self.__products.setAncestor( self )

    @property
    def key( self ) :

        return( self.__label )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.__label = value

    @property
    def label( self ) :

        return( self.__label )

    @property
    def probability( self ) :

        return( self.__probability )

    @property
    def products( self ) :

        return( self.__products )

    @property
    def Q( self ) :

        return( self.__Q )

    def convertUnits( self, unitMap ) :

        self.__probability.convertUnits( unitMap )
        self.__Q.convertUnits( unitMap )
        self.__products.convertUnits( unitMap )

    def copy( self ) :

        _decay = self.__class__( self.label )
        self.__copyStandardQuantities( _decay )
        return( _decay )

    def __copyStandardQuantities( self, other ) :

        if( other.moniker != self.moniker ) : TypeError( '''other = "%s" not same as self = "%s".''' %
            ( other.moniker, self.moniker ) )
        for item in self.__probability : other.probability.add( item.copy( ) )
        for item in self.__Q : other.Q.add( item.copy( ) )
        for item in self.__products : other.products.add( item.copy( ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        XMLStringList += self.probability.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.Q.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.products.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if( child.tag == probabilityModule.suite.moniker ) :
                self.probability.parseXMLNode( child, xPath, linkData )
            elif( child.tag == QModule.suite.moniker ) :
                self.Q.parseXMLNode( child, xPath, linkData )
            elif( child.tag == productModule.suite.moniker ) :
                self.products.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

class suite( suiteModule.suite ) :

    moniker = 'decays'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( decay, ), key = 'label' )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            _decay = decay( child.attrib['label'] )
            self.add( _decay.parseXMLNode( child, xPath, linkData ) )

        xPath.pop( )
        return( self )
