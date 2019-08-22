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

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from . import base as baseModule
from . import axes as axesModule

class constant( baseModule.xDataFunctional )  :

    def __init__( self, _constant, domainMin, domainMax, axes = None, label = None ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, label = label, axes = axes )

        self.constant = _constant

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : TypeError( 'domainMin not a float instance' )
        self.__domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : TypeError( 'domainMax not a float instance' )
        self.__domainMax = domainMax

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = self.axes.copy( )
        return( self.__class__( self.constant, self.domainMin, self.domainMax, axes = axes, label = self.label ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    @property
    def constant( self ) :

        return( self.__constant )

    @constant.setter
    def constant( self, _constant ) :

        if( isinstance( _constant, int ) ) : _constant = float( _constant )
        if( not( isinstance( _constant, float ) ) ) : TypeError( 'constant not a float instance' )
        self.__constant = _constant

    @property
    def domainMin( self ) :

        return( self.__domainMin )

    @property
    def domainMax( self ) :

        return( self.__domainMax )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :

        return( self.constant )

    rangeMax = rangeMin

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def fixDomainPerUnitChange( self, factors ) :

        self.__domainMin *= factors[self.dimension]
        self.__domainMax *= factors[self.dimension]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        valueFormatter = kwargs.get( 'valueFormatter', floatToShortestString )
        significantDigits = kwargs.get( 'significantDigits', 15 )

        attributeStr = baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s%s constant="%s" domainMin="%s" domainMax="%s">' % ( indent, self.moniker, attributeStr, 
                valueFormatter( self.constant, significantDigits = significantDigits ),
                valueFormatter( self.domainMin, significantDigits = significantDigits ),
                valueFormatter( self.domainMax, significantDigits = significantDigits ) ) ]

        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent = indent2, **kwargs )

        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None, **kwargs ) :

        attrs =      { 'constant' :  None, 'domainMin' :  None, 'domainMax' :  None, 'label' : None }
        attributes = { 'constant' : float, 'domainMin' : float, 'domainMax' : float, 'label' : str }
        for key, item in xDataElement.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        uncertainty = None
        for subElement in xDataElement :
            if( subElement.tag == 'axes' ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'uncertainty' ) :
                from . import uncertainties as uncertaintiesModule
                uncertainty = uncertaintiesModule.uncertainty.parseXMLNode( subElement, xPath, linkData )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        _constant = attrs.pop( 'constant' )
        newConstant = cls( _constant, axes = axes, **attrs )
        newConstant.uncertainty = uncertainty
        return newConstant

    @staticmethod
    def parseXMLString( XMLString ) :

        from xml.etree import cElementTree

        return( constant.parseXMLNode( cElementTree.fromstring( XMLString ), xPath = [], linkData = {} ) )

class constant1d( constant ) :

    moniker = 'constant1d'
    dimension = 1

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        factors = self.axes.convertUnits( unitMap )
        self.constant *= factors[0]
        self.fixDomainPerUnitChange( factors )
        self.fixValuePerUnitChange( factors )

    def evaluate( self, x ) :

        return( self.constant )
