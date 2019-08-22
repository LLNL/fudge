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

from pqu import PQU as PQUModule

import ancestry as ancestryModule

class physicalQuantity( ancestryModule.ancestry ) :

    moniker = 'physicalQuantity'

    def __init__( self, value, unit, label = None ) :

        ancestryModule.ancestry.__init__( self )
        self.__PQ = PQUModule.PQU( value, unit )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

    def __str__( self ) :

        return( self.toString( ) )

    def __float__( self ) :

        return( float( self.__PQ ) )

    def __eq__( self, other ) :
        """
        If other is not a physicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, physicalQuantity ) ) : return( False )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( False )
        return( self.compare( other, 5 ) == 0 )

    def __ne__( self, other ) :
        """
        If other is not a physicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, physicalQuantity ) ) : return( True )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( True )
        return( self.compare( other, 5 ) != 0 )

    def __lt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) < 0 )

    def __le__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) <= 0 )

    def __gt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) > 0 )

    def __ge__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) >= 0 )

    def compare( self, other, epsilonFactor = 0 ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        if( not isinstance( other, physicalQuantity ) ) : raise TypeError( 'other not instance of physicalQuantity' )
        return( self.__PQ.compare( other.__PQ, epsilonFactor ) )

    @property
    def key( self ) :

        return( self.__label )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def value( self ) :

        return( float( self ) )

    @property
    def unit( self ) :

        return( self.__PQ.unit )

    def convertToUnit( self, unit ) :

        self.__PQ.convertToUnit( unit )

    def convertUnits( self, unitMap ) :

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        self.__PQ.convertToUnit( unit )

    def copy( self ) :

        return( self.__class__( self.value, self.unit, self.label ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    def copyToUnit( self, unit ) :

        _copy = self.copy( )
        _copy.convertToUnit( unit )
        return( _copy )

    def getValueAs( self, unit ) :

        return( self.__PQ.getValueAs( unit ) )

    def toString( self, significantDigits = None, keepPeriod = True ) :

        return( self.__PQ.toString( significantDigits = significantDigits, keepPeriod = keepPeriod ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )

        label = ''
        if( self.label is not None ) : label = ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s" unit="%s"/>' % ( indent, self.moniker, label, self.value, self.unit ) ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        value = element.get( 'value' )
        unit = element.get( 'unit' )
        label = element.get( 'label', None )
        _cls = cls( value, unit, label )

        xPath.pop( ) 
        return( _cls )
