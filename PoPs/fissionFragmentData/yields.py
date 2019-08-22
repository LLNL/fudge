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

from xData import ancestry as ancestryModule
from xData import values as valuesModule
from xData import array as arrayModule

from .. import misc as miscModule
from .. import suite as suiteModule

class values( miscModule.classWithLabelKey ) :

    moniker = 'values'

    def __init__( self, _values ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( _values, ( tuple, list ) ) ) ) : raise TypeError( 'values must be a values instance.' )
        self.__values = tuple( value for value in _values )

    def __getitem__( self, index ) :

        return( self.__values[index] )

    def __iter__( self ) :

        n1 = len( self.__values )
        for i1 in range( n1 ) : yield self.__values[i1]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList[-1] += ' '.join( [ '%s' % value for value in self.__values ] )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass(cls, element, xPath, linkData):

        xPath.append( element.tag )
        values_ = cls( map(float, element.text.split()) )
        xPath.pop()
        return values_

class covariance( ancestryModule.ancestry ) :

    moniker = 'covariance'

    def __init__( self, _matrix ) :

        ancestryModule.ancestry.__init__( self )

        self.__matrix = _matrix

    def matrix( self ) :

        return( self.__matrix )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( 'valuesPerLine' not in kwargs ) : kwargs['valuesPerLine'] = 1000
        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__matrix.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == arrayModule.arrayBase.moniker:
            _matrix = arrayModule.arrayBase.parseXMLNode(child, xPath, linkData)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        covar_ = cls(_matrix)
        xPath.pop()
        return covar_

class uncertainty( ancestryModule.ancestry ) :

    moniker = 'uncertainty'

    def __init__( self, form ) :

        ancestryModule.ancestry.__init__( self )

        self.__form = form

    def form( self ) :

        return( self.__form )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__form.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == covariance.moniker:
            form = covariance.parseXMLNodeAsClass(child, xPath, linkData)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        uncert_ = cls(form)
        xPath.pop()
        return uncert_

class yields( ancestryModule.ancestry ) :

    moniker = 'yields'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__values = None
        self.__uncertainty = None

    @property
    def values( self ) :

        return( self.__values )

    @values.setter
    def values( self, _values ) :

        if( not( isinstance( _values, values ) ) ) : raise TypeError( 'Invalid values instance.' )
        self.__values = _values
        self.__values.setAncestor( self )

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( not( isinstance( _uncertainty, uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )
        self.__uncertainty = _uncertainty
        self.__uncertainty.setAncestor( self )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        if( self.__values is not None ) : XMLStringList += self.__values.toXMLList( indent2, **kwargs )
        if( self.__uncertainty is not None ) : XMLStringList += self.__uncertainty.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append(element.tag)
        for child in element:
            if child.tag == values.moniker:
                self.values = values.parseXMLNodeAsClass(child, xPath, linkData)
            elif child.tag == uncertainty.moniker:
                self.uncertainty = uncertainty.parseXMLNodeAsClass(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return (self)

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseXMLNode(element, xPath, linkData)

        return self
