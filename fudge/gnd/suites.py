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

__metaclass__ = type

import abc
import xData.ancestry as ancestryModule

from fudge.core.utilities import brb

class suite( ancestryModule.ancestry ) :
    """
    Base class for a class member that is list like. For example, the lists inside the class 
    reactionSuite ('reactions', 'sums', 'productions' and 'fissionComponents').
    """
    __metaclass__ = abc.ABCMeta

    def __init__( self, allowedClasses, replace = False ) :

        ancestryModule.ancestry.__init__( self )
        self.__allowedClasses = [ cls for cls in allowedClasses ]
        self.__replace = replace
        self.__items = []

    def __contains__( self, label ) :

        for item in self :
            if( item.label == label ) : return( True )
        return( False )

    def __getitem__( self, label ) :

        if( isinstance( label, int ) ) : return( self.__items[label] )   # BRB - FIXME, Temp fix until Caleb gets neutrons/n-017_Cl_035.endf working without it.
        if( not( isinstance( label, str ) ) ) : raise TypeError( "label must be a string" )
        for item in self :
            if( item.label == label ) : return( item )
        # requested style not found, but what about styles it derives from?
        requestedStyle = self.getRootAncestor().styles[ label ]
        for dstyle in requestedStyle.derivedStyles:
            for item in self :
                if( item.label == dstyle.label ) : return( item )
        raise KeyError( "item with label '%s' not found in suite '%s'" % ( label, self.moniker ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__items[i1]

    def __len__( self ) :

        return( len( self.__items ) )

    @property
    def allowedClasses( self ) :

        return( self.__allowedClasses )

    @property
    def replace( self ) :

        return( self.__replace )

    def add( self, newItem ) :

        found = False
        for cls in self.__allowedClasses :
            if( isinstance( newItem, cls ) ) :
                found = True
                break
        if( not( found ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( newItem.__class__, self.moniker ) )

        index = len( self )
        for i1, item in enumerate( self.__items ) :
            if( item.label == newItem.label ) :
                if( self.replace ) :
                    index = i1
                    break
                else :
                    raise KeyError( 'item with label = "%s" already present in suite' % item.label )

        newItem.setAncestor( self, attribute = 'label' )
        self.__items.insert( index, newItem )

    def remove( self, label ) :

        for i1, item in enumerate( self.__items ) :
            if( item.label == label ) :
                del self.__items[i1]
                return( True )
        return( False )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( len( self ) == 0 ) : return( [] )
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for item in self : xmlString += item.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            parseClass = None
            for class_ in self.__allowedClasses:
                if (child.tag == class_.moniker):
                    parseClass = class_
                    break
            if parseClass is None:
                raise TypeError( "Invalid element '%s' encountered in suite '%s'" % (child.tag, self.moniker) )

            self.add( parseClass.parseXMLNode( child, xPath, linkData ) )

        xPath.pop()

class reactions( suite ) :

    moniker = 'reactions'

    def __init__( self ) :

        from fudge.gnd.reactions import reaction as reactionModule
        from fudge.gnd.reactions import production as productionModule

        suite.__init__( self, [reactionModule.reaction, productionModule.production] )

class sums( suite ) :

    moniker = 'sums'

    def __init__( self ) :

        from fudge.gnd import sums as sumsModule

        suite.__init__( self, [sumsModule.crossSectionSum, sumsModule.multiplicitySum] )

class productions( suite ) :

    moniker = 'productions'

    def __init__( self ) :

        from fudge.gnd.reactions import production as productionModule

        suite.__init__( self, [productionModule.production] )

class fissionComponents( suite ) :

    moniker = 'fissionComponents'

    def __init__( self ) :

        from fudge.gnd.reactions import fissionComponent as fissionComponentModule

        suite.__init__( self, [fissionComponentModule.fissionComponent] )
