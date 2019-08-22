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

import string
import abc

import xData.ancestry as ancestryModule

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

    def __delitem__( self, key ) :

        if( isinstance( key, int ) ) :
            del self.__items[key]
            return
        status = self.remove( key )
        if( not( status ) ) : KeyError( 'key not in suite' )

    def __getitem__( self, label ) :

        if( isinstance( label, int ) ) : return( self.__items[label] )
        if( not( isinstance( label, str ) ) ) : raise TypeError( "label must be a string" )
        for item in self :
            if( item.label == label ) : return( item )

        # requested style not found, but what about styles it derives from?
        requestedStyle = self.getRootAncestor().styles[ label ]

        for item in self :
            if( item.label == requestedStyle.label ) : return( item )
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
        """
        Add newItem to the suite. If another item in the suite has the same label as newItem, raise KeyError
        :param newItem:
        :return:
        """

        if( not( isinstance( newItem.label, str ) ) ) : IndexError( '''Item's label must be a string''' )
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

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        for index, item in enumerate( self.__items ) : item.convertUnits( unitMap )

    def checkAncestry( self, verbose = 0, level = 0 ) :

        for item in self : item.checkAncestry( verbose = verbose, level = level )

    def pop( self, label, *args ) :
        """
        Remove item by label, and return it. If the label is not found, raise KeyError
        :param label:
        :param args:
        :return:
        """

        if( len( args ) > 2 ) : raise Exception( 'Only one default value is allowed: got %s' % len( args ) )
        for i1, item in enumerate( self.__items ) :
            if( item.label == label ) : return( self.__items.pop( i1 ) )
        if( len( args ) == 1 ) : return( args[0] )
        raise KeyError( label )

    def remove( self, label ) :
        """
        Remove item by label. Returns True if label was present, otherwise returns False
        :param label: str
        :return: bool
        """

        for i1, item in enumerate( self.__items ) :
            if( item.label == label ) :
                del self.__items[i1]
                return( True )
        return( False )

    def clear( self ):
        """
        Remove all members of the suite
        :return:
        """
        for idx in range( len( self.__items ) ):
            self.__items.pop()

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( len( self ) == 0 ) : return( [] )
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for item in self : xmlString += item.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def uniqueLabel( self, item ) :
        """
        If item's label is the same as another item's label in self, construct a new unique label
        based on item's label appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        def incrementSuffix( suffix ) :

            if( len( suffix ) == 0 ) : return( 'a' )
            index = string.ascii_lowercase.find( suffix[-1] ) + 1
            if( index != 26 ) : return( suffix[:-1] + string.ascii_lowercase[index] )
            return( incrementSuffix( suffix[:-1] ) + 'a' )

        if( item.label in self ) :
            label = item.label
            label__ = label + '__'
            n1 = len( label__ )
            l1 = 0
            suffixes = []
            for _item in self :          # Find list of longest labels that start with label__.
                if( _item.label[:n1] == label__ ) :
                    suffix = _item.label[n1:]
                    if( not( suffix.islower( ) ) ) : continue       # Ignore non-standard labels.
                    l2 = len( suffix )
                    if( l2 < l1 ) : continue
                    if( l2 > l1 ) :
                        l1 = l2
                        suffixes = []
                    suffixes.append( suffix )
            if( len( suffixes ) == 0 ) :
                suffix = 'a'
            else :
                suffix = incrementSuffix( sorted( suffixes )[-1] )
            item.label = label__ + suffix
        return( item )

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

        suite.__init__( self, [ reactionModule.reaction ] )

class orphanProducts( suite ) :

    moniker = 'orphanProducts'

    def __init__( self ) :

        from fudge.gnd.reactions import reaction as reactionModule

        suite.__init__( self, [ reactionModule.reaction ] )

class crossSections( suite ) :

    moniker = 'crossSections'

    def __init__( self ) :

        from fudge.gnd import sums as sumsModule

        suite.__init__( self, [sumsModule.crossSectionSum] )

class multiplicities( suite) :

    moniker = 'multiplicities'

    def __init__( self ) :

        from fudge.gnd import sums as sumsModule

        suite.__init__( self, [sumsModule.multiplicitySum] )

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
