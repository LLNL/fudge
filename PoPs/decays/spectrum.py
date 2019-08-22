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
from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule
from xData import XYs as XYsModule

from .. import suite as suiteModule
from .. import misc as miscModule

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

#
# FIXME Need a physicalQuantity class that has keyName.
#
class unitless( physicalQuantityModule.physicalQuantity ) :

    keyName = 'label'

    def __init__( self, value, label = None ) :

        physicalQuantityModule.physicalQuantity.__init__( self, value, '', label = label )

    def copy( self ) :              # overrides required since this __init__ takes different arguments:

        cls = self.__class__( self.value, self.label )
        if( self.uncertainty is not None ) : cls.uncertainty = self.uncertainty.copy( )
        return( cls )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        value = element.get( 'value' )
        label = element.get( 'label', None )
        _cls = cls( value, label )
        for child in element :
            if( child.tag == uncertaintyModule.uncertainty.moniker ) :
                _cls.uncertainty = uncertaintyModule.uncertainty.parseXMLNodeAsClass( child, xPath, linkData )

        xPath.pop( )
        return( _cls )

    @classmethod    # Need parseXMLNodeAsClass to handle internalConversionCoefficients
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        _cls = cls.parseXMLNode( element, xPath, linkData )
        return( _cls )

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

class photonEmissionProbabilities( suiteModule.suite ) :

    moniker = 'photonEmissionProbabilities'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( shell, ) )

class internalPairFormationCoefficient( unitless ) :

    moniker = 'internalPairFormationCoefficient'

class positronEmissionIntensity( unitless ) :

    moniker = 'positronEmissionIntensity'

class discrete( ancestryModule.ancestry ) :

    moniker = 'discrete'

    def __init__( self, _intensity, _energy, type = None, _internalPairFormationCoefficient = None,
                  _positronEmissionIntensity = None ) :

        ancestryModule.ancestry.__init__(self)
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

        if( _positronEmissionIntensity is not None ) :
            if( not( isinstance( _positronEmissionIntensity, positronEmissionIntensity ) ) ) :
                raise TypeError( '_positronEmissionIntensity must be an instance of positronEmissionIntensity' )
        self.__positronEmissionIntensity = _positronEmissionIntensity


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
    def positronEmissionIntensity( self ) :

        return( self.__positronEmissionIntensity )

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

        _discrete = self.__class__( self.__intensity.copy( ), self.__energy.copy( ), type = self.__type, _internalPairFormationCoefficient = self.__internalPairFormationCoefficient )
        for icc in self.internalConversionCoefficients : _discrete.internalConversionCoefficients.add( icc.copy( ) )
        return _discrete

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
        if( self.__positronEmissionIntensity is not None): XMLStringList += self.__positronEmissionIntensity.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        kwargs = {'type': element.get('type')}
        ICCelement = None
        for child in element :
            found = False
            for subclass in ( intensity, energy, internalPairFormationCoefficient, positronEmissionIntensity ) :
                if( child.tag == subclass.moniker ) :
                    kwargs['_'+subclass.moniker] = subclass.parseXMLNode( element.find( subclass.moniker ), xPath, linkData )
                    found = True
            if not found:
                if child.tag == internalConversionCoefficients.moniker:
                    ICCelement = child
                else:
                    raise ValueError("Encountered unexpected child element '%s'" % child.tag)
        _discrete = cls( **kwargs )
        if ICCelement is not None:
            _discrete.internalConversionCoefficients.parseXMLNode( ICCelement, xPath, linkData )
        xPath.pop()
        return _discrete

class continuum( ancestryModule.ancestry ) :

    moniker = 'continuum'

    def __init__( self, spectrum ) :

        ancestryModule.ancestry.__init__( self )
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

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )
        _spectrum = XYsModule.XYs1d.parseXMLNode( element.find(XYsModule.XYs1d.moniker), xPath, linkData )
        instance = continuum( _spectrum )
        xPath.pop()
        return instance

    @staticmethod
    def defaultAxes( energyUnit ) :

        axes = axesModule.axes( rank = 2 )
        axes[0] = axesModule.axis( "P(energy_out)", 0, "1/%s" % energyUnit )
        axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
        return( axes )

class spectrum( miscModule.classWithLabelKey ) :

    moniker = 'spectrum'

    def __init__( self, label, pid ) :

        miscModule.classWithLabelKey.__init__( self, label )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid not str' )
        self.__pid = pid

        self.__emissions = []

    def __len__( self ) :

        return( len( self.__emissions ) )

    def __getitem__( self, index ) :

        return( self.__emissions[index] )

    @property
    def pid( self ) :

        return( self.__pid )

    def append( self, emission ) :

        if( not( isinstance( emission, ( discrete, continuum ) ) ) ) : raise TypeError( 'emission must be instance of discrete or continuum' )
        self.__emissions.append( emission )

    def convertUnits( self, unitMap ) :

        for emission in self.__emissions : emission.convertUnits( unitMap )

    def copy( self ) :
        """Returns a deep copy of self."""

        _spectrum = self.__class__( self.label, self.pid )
        for __emission in self.__emissions : _spectrum.append( __emission.copy( ) )
        return( _spectrum )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" pid="%s">' % ( indent, self.moniker, self.label, self.pid ) ]
        for _emission in self.__emissions : XMLStringList += _emission.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )
        for child in element:
            if( child.tag == discrete.moniker ) :
                self.append( discrete.parseXMLNode( child, xPath, linkData ) )
            elif( child.tag == continuum.moniker ) :
                self.append( continuum.parseXMLNode( child, xPath, linkData ) )
            else :
                raise ValueError("Encountered unexpected child element '%s'" % child.tag)
        xPath.pop( )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )
        _spectrum = cls( element.get('label'), element.get('pid') )
        xPath.pop()

        _spectrum.parseXMLNode( element, xPath, linkData )
        return _spectrum

class spectra( suiteModule.suite ) :

    moniker = 'spectra'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( spectrum, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

#
# FIXME, this method is broken. What is product?
#
        for child in element :
            _spectrum = spectrum.parseXMLNodeAsClass( child, xPath, linkData )
            self.add( _spectrum )

        xPath.pop( )
        return( self )
