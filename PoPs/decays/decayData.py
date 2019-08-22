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
This module contains the particle decay classes.
"""

import abc

from xData import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import misc as miscModule

from ..decays import probability as probabilityModule
from ..decays import Q as QModule
from ..decays import product as productModule
from ..decays import spectrum as spectrumModule
from ..decays import averageEnergy as averageEnergyModule

decayModesIT = "isomeric transition"
decayModesSF = "spontaneous fission"
decayModesParticle = ""

decayModeTypes = [decayModesIT, decayModesSF, decayModesParticle]

class decay( miscModule.classWithIndexKey ) :

    moniker = 'decay'

    def __init__( self, index, type, complete = True ) :

        miscModule.classWithIndexKey.__init__( self, index )

        if( not( isinstance( type, str ) ) ) : raise TypeError( 'decay mode type not an instance of str' )
        if( type not in decayModeTypes ) : raise ValueError( 'decay mode type invalid' )
        self.__type = type

        if( not( isinstance( complete, bool ) ) ) : raise TypeError( 'complete not an instance of bool' )
        self.__complete = complete

        self.__products = productModule.suite( )
        self.__products.setAncestor( self )

    @property
    def complete( self ) :

        return( self.__complete )

    @property
    def products( self ) :

        return( self.__products )

    @property
    def type( self ) :

        return( self.__type )

    def convertUnits( self, unitMap ) :

        self.__products.convertUnits( unitMap )

    def copy( self ) :

        _decay = self.__class__( self.index, self.type, self.complete )
        for item in self.__products : _decay.products.add( item.copy( ) )
        return( _decay )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        if( self.type != '' ) : attributes = ' type="%s"' % self.type
        if( not( self.complete ) ) : attributes = ' complete="false"'
        XMLStringList = [ '%s<%s index="%s"%s>' % ( indent, self.moniker, self.index, attributes ) ]
        XMLStringList += self.products.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        self.products.parseXMLNode( element.find( self.products.moniker ), xPath, linkData )

        xPath.pop()
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        complete = True
        if element.attrib.get('complete') == 'false': complete = False
        self = cls( element.attrib['index'], element.attrib.get( 'type', '' ), complete )
        xPath.pop( )
        self.parseXMLNode( element, xPath, linkData )

        return( self )

class decayPath( suiteModule.suite ) :

    moniker = 'decayPath'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( decay, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( decay.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )

class decayMode( miscModule.classWithLabelKey ) :

    moniker = 'decayMode'

    def __init__( self, label, mode ) :

        miscModule.classWithLabelKey.__init__( self, label )

        if( not( isinstance( mode, str ) ) ) : raise TypeError( 'mode not an instance of str' )
        self.__mode = mode

        self.__probability = probabilityModule.suite( )
        self.__probability.setAncestor( self )

        self.__internalConversionCoefficients = spectrumModule.internalConversionCoefficients()
        self.__internalConversionCoefficients.setAncestor( self )

        self.__photonEmissionProbabilities = spectrumModule.photonEmissionProbabilities()
        self.__photonEmissionProbabilities.setAncestor( self )

        self.__Q = QModule.suite( )
        self.__Q.setAncestor( self )

        self.__decayPath = decayPath( )
        self.__decayPath.setAncestor( self )

        self.__spectra = spectrumModule.spectra( )
        self.__spectra.setAncestor( self )

    @property
    def mode( self ) :

        return( self.__mode )

    @property
    def probability( self ) :

        return( self.__probability )

    @property
    def internalConversionCoefficients( self ):

        return( self.__internalConversionCoefficients )

    @property
    def photonEmissionProbabilities( self ):

        return( self.__photonEmissionProbabilities )

    @property
    def decayPath( self ) :

        return( self.__decayPath )

    @property
    def Q( self ) :

        return( self.__Q )

    @property
    def spectra( self ) :

        return( self.__spectra )

    def convertUnits( self, unitMap ) :

        self.__probability.convertUnits( unitMap )
        self.__Q.convertUnits( unitMap )
        self.__decayPath.convertUnits( unitMap )
        self.__spectra.convertUnits( unitMap )

    def copy( self ) :

        other = self.__class__( self.label, self.mode )
        for item in self.__probability : other.probability.add( item.copy( ) )
        for item in self.__Q : other.Q.add( item.copy( ) )
        for item in self.__decayPath : other.decayPath.add( item.copy( ) )
        for item in self.__spectra : other.spectra.add( item.copy( ) )
        return( other )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" mode="%s">' % ( indent, self.moniker, self.label, self.mode ) ]
        XMLStringList += self.probability.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.internalConversionCoefficients.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.photonEmissionProbabilities.toXMLList(indent=indent2, **kwargs)
        XMLStringList += self.Q.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.decayPath.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.spectra.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )

        for child in element :
            if( child.tag == probabilityModule.suite.moniker ) :
                self.probability.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.internalConversionCoefficients.moniker ) :
                self.internalConversionCoefficients.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.photonEmissionProbabilities.moniker ) :
                self.photonEmissionProbabilities.parseXMLNode( child, xPath, linkData )
            elif( child.tag == QModule.suite.moniker ) :
                self.Q.parseXMLNode( child, xPath, linkData )
            elif( child.tag == decayPath.moniker ) :
                self.decayPath.parseXMLNode( child, xPath, linkData )
            elif( child.tag == spectrumModule.spectra.moniker ) :
                self.spectra.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( element.attrib['label'], element.attrib['mode'] )
        self.parseXMLNode( element, xPath, linkData )

        xPath.pop( )
        return( self )

class decayModes( suiteModule.suite ) :

    moniker = 'decayModes'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( decayMode, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( decayMode.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )

class decayData( ancestryModule.ancestry ) :

    moniker = 'decayData'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__decayModes = decayModes( )
        self.__decayModes.setAncestor( self )

        self.__averageEnergies = averageEnergyModule.averageEnergies( )
        self.__averageEnergies.setAncestor( self )

    @property
    def decayModes( self ) :

        return( self.__decayModes )

    @property
    def averageEnergies( self ) :

        return( self.__averageEnergies )

    def convertUnits( self, unitMap ) :

        self.__decayModes.convertUnits( unitMap )
        self.__averageEnergies.convertUnits( unitMap )

    def copy( self ) :

        _decayData = decayData( )
        self.copyItems( _decayData )
        return( _decayData )

    def copyItems( self, other ) :

        for item in self.__decayModes : other.decayModes.add( item.copy( ) )
        for item in self.__averageEnergies : other.averageEnergies.add( item.copy( ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList_subs  = self.__decayModes.toXMLList( indent = indent2, **kwargs )
        XMLStringList_subs += self.averageEnergies.toXMLList( indent = indent2, **kwargs )
        if( len( XMLStringList_subs ) == 0 ) : return( [] )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ] + XMLStringList_subs
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        children = { 'decayModes' : self.decayModes, 'averageEnergies' : self.averageEnergies }

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseXMLNode( child, xPath, linkData )
            else:
                if( not( self.parseExtraXMLElement( child, xPath, linkData ) ) ) :
                    raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )
