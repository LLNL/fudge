# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Angular/Energy double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import angular as angularModule

__metaclass__ = type

class xs_pdf_cdf1d( xs_pdf_cdfModule.xs_pdf_cdf1d ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( xs_pdf_cdf1d, ) )

class XYs3d( multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class subform( baseModule.subform ) :
    """Abstract base class for angularEnergyMC subforms."""

    def __init__( self, data ) :

        baseModule.subform.__init__( self )
        if( not( isinstance( data, self.allowedSubElements ) ) ) : raise TypeError( 'Invalid instance: %s' % type( data ) )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.data.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class angular( subform ) :

    moniker = 'angular'
    allowedSubElements = ( angularModule.XYs2d, )
    ancestryMembers = ( 'angular', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    angularModule.XYs2d.moniker       : angularModule.XYs2d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _angular= angular( angularSubform )

        xPath.pop( )
        return( _angular )

class angularEnergy( subform ) :

    moniker = 'angularEnergy'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'angularEnergy', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    XYs3d.moniker       : XYs3d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angularEnergy subform "%s"' % subformElement.tag )
        angularEnergySubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _angularEnergy = angularEnergy( angularEnergySubform )

        xPath.pop( )
        return( _angularEnergy )

class form( baseModule.form ) :

    moniker = 'angularEnergyMC'
    subformAttributes = ( 'angular', 'angularEnergy' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angular, _angularEnergy ) :

        if( not( isinstance( _angular, angular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angular ) )
        if( not( isinstance( _angularEnergy, angularEnergy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angularEnergy ) )
        baseModule.form.__init__( self, label, productFrame, ( _angular, _angularEnergy ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angular.convertUnits( unitMap )
        self.angularEnergy.convertUnits( unitMap )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <angularEnergyMC> element from xml."""

        xPath.append( element.tag )

        _angular = None
        _angularEnergy = None
        for child in element :
            if( child.tag == angular.moniker ) :
                _angular = angular.parseXMLNode( child, xPath, linkData )
            elif( child.tag == angularEnergy.moniker ) :
                _angularEnergy = angularEnergy.parseXMLNode( child, xPath, linkData )
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        angularEnergyMC = form( element.get( "label" ), element.get('productFrame'), 
                _angular, _angularEnergy )

        xPath.pop( )
        return( angularEnergyMC ) 
