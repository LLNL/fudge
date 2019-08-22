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
This module contains the spectrum classes.
"""

import abc

from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import physicalQuantity as physicalQuantityModule

from .. import suite as suiteModule

gammaProduct = "gamma"
betaMinusProduct = "beta-"
betaPlusProduct = "beta+"
betaPlusOrElectronCaptureProduct = "beta+ or electronCapture"
alphaProduct = "alpha"
neutronProduct = "neutron"
SFProduct = "spontaneous fission fragments"
protronProduct = "proton"
discreteElectronProduct = "discrete electron"
xRayProduct = "x-ray"

class unitless( physicalQuantityModule.physicalQuantity ) :

    def __init__( self, value, label = None ) :

        physicalQuantityModule.physicalQuantity.__init__( self, value, '', label = label )

class transitionType :

    allowed = 'allowed'
    firstForbidden = 'first-forbidden'
    secondForbidden = 'second-forbidden'

    types = [ allowed, firstForbidden, secondForbidden ]

class intensity( unitless ) :

    moniker = 'intensity'

class energy( physicalQuantityModule.physicalQuantity ) :

    moniker = 'energy'

class shell( unitless ) :

    moniker = 'shell'
    total = 'total'
    KShell = 'K'
    LShell = 'L'

class internalConversionCoefficients( suiteModule.suite ) :

    moniker = 'internalConversionCoefficients'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( shell, ) )

class internalPairFormationCoefficient( unitless ) :

    moniker = 'internalPairFormationCoefficient'

class discrete( ancestryModule.ancestry ) :

    moniker = 'discrete'

    def __init__( self, _intensity, _energy, type = None, _internalPairFormationCoefficient = None ) :

        if( not( isinstance( _intensity, intensity ) ) ) : raise TypeError( '_intensity must be an instance of intensity' )
        self.__intensity = _intensity

        if( not( isinstance( _energy, energy ) ) ) : raise TypeError( '_energy must be an instance of energy' )
        self.__energy = _energy

        if( type is not None ) :
            if( type not in transitionType.types ) : raise ValueError( 'invalid type' )
        self.__type = type

        self.__internalConversionCoefficients = internalConversionCoefficients( )
        self.__internalConversionCoefficients.setAncestor( self )

        if( _internalPairFormationCoefficient is not None ) :
            if( not( isinstance( _internalPairFormationCoefficient, internalPairFormationCoefficient ) ) ) :
                raise TypeError( '_internalPairFormationCoefficient must be an instance of internalPairFormationCoefficient' )
        self.__internalPairFormationCoefficient = _internalPairFormationCoefficient

    @property
    def energy( self ) :

        return( self.__energy )

    @property
    def internalConversionCoefficients( self ) :

        return( self.__internalConversionCoefficients )

    @property
    def internalPairFormationCoefficient( self ) :

        return( self.__internalPairFormationCoefficient )

    @property
    def intensity( self ) :

        return( self.__intensity )

    @property
    def type( self ) :

        return( self.__type )

    def convertUnits( self, unitMap ) :

        self.__energy.convertUnits( unitMap )
        self.__intensity.convertUnits( unitMap )

    def copy( self ) :

        _discrete = self.__class__( self.__energy.copy( ), self.__energy.copy( ), type = self.__type, _internalPairFormationCoefficient = self.__internalPairFormationCoefficient )
        for icc in self.internalConversionCoefficients : _discrete.internalConversionCoefficients.add( icc.copy( ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        type = ''
        if( self.type is not None ) : type = ' type="%s"' % self.type
        XMLStringList = [ '%s<%s%s>' % ( indent, self.moniker, type ) ]
        XMLStringList += self.__intensity.toXMLList( indent2, **kwargs )
        XMLStringList += self.__energy.toXMLList( indent2, **kwargs )
        XMLStringList += self.__internalConversionCoefficients.toXMLList( indent2, **kwargs )
        if( self.__internalPairFormationCoefficient is not None ) : XMLStringList += self.__internalPairFormationCoefficient.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class continuum( ancestryModule.ancestry ) :

    moniker = 'continuum'

    def __init__( self, spectrum ) :

        self.__spectrum = spectrum

    @property
    def spectrum( self ) :

        return( self.__spectrum )

    def convertUnits( self, unitMap ) :

        self.__spectrum.convertUnits( unitMap )

    def copy( self ) :

        return( self.__class__( self.__spectrum.copy( ) ) )

    def toXML( self, indent = "", **kwargs ) :
    
        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__spectrum.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    @staticmethod
    def defaultAxes( energyUnit ) :

        axes = axesModule.axes( rank = 2 )
        axes[0] = axesModule.axis( "P(energy_out)", 0, "1/%s" % energyUnit )
        axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
        return( axes )

class spectrum( ancestryModule.ancestry ) :

    moniker = 'spectrum'

    def __init__( self, label, pid ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label not str' )
        self.__label = label

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid not str' )
        self.__pid = pid

        self.__spectra = []

    def __len__( self ) :

        return( len( self.__spectra ) )

    def __getitem__( self, index ) :

        return( self.__spectra[index] )

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
    def pid( self ) :

        return( self.__pid )

    def append( self, spectrum ) :

        if( not( isinstance( spectrum, ( discrete, continuum ) ) ) ) : raise TypeError( 'spectrum must be instance of discrete or continuum' )
        self.__spectra.append( spectrum )

    def convertUnits( self, unitMap ) :

        for spectrum in self.__spectra : spectrum.convertUnits( unitMap )

    def copy( self ) :
        """Returns a deep copy of self."""

        _spectrum = self.__class__( self.label, self.pid )
        for __spectrum in self.__spectra : _spectrum.addSpectrum( __spectrum.copy( ) )
        return( _spectrum )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" pid="%s">' % ( indent, self.moniker, self.label, self.pid ) ]
        for _spectrum in self.__spectra : XMLStringList += _spectrum.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        xPath.pop( )
        return( self )

class spectra( suiteModule.suite ) :

    moniker = 'spectra'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( spectrum, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            _product = product( child.attrib['label'], child.attrib['pid'] )
            self.add( _product.parseXMLNode( child, xPath, linkData ) )

        xPath.pop( )
        return( self )
