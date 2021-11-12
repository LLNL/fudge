# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Energy/angular double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import energy as energyModule

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
    """Abstract base class for energyAngular subforms."""

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

class energy( subform ) :

    moniker = 'energy'
    allowedSubElements = ( energyModule.XYs2d, )
    ancestryMembers = ( 'energy', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    energyModule.XYs2d.moniker       : energyModule.XYs2d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown energy subform "%s"' % subformElement.tag )
        energySubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _energy = energy ( energySubform )

        xPath.pop( )
        return( _energy )

class energyAngular( subform ) :

    moniker = 'energyAngular'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'energyAngular', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    XYs3d.moniker       : XYs3d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown energyAngular subform "%s"' % subformElement.tag )
        energyAngularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _energyAngular = energyAngular( energyAngularSubform )

        xPath.pop( )
        return( _energyAngular )

class form( baseModule.form ) :

    moniker = 'energyAngularMC'
    subformAttributes = ( 'energy', 'energyAngular' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _energy, _energyAngular ) :

        if( not( isinstance( _energy, energy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energy ) )
        if( not( isinstance( _energyAngular, energyAngular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energyAngular ) )
        baseModule.form.__init__( self, label, productFrame, ( _energy, _energyAngular ) )

    @property
    def domainMin( self ) :

        return( self.energy.data.domainMin )

    @property
    def domainMax( self ) :

        return( self.energy.data.domainMax )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <energyAngularMC> element from xml."""

        xPath.append( element.tag )

        _energy = None
        _energyAngular = None
        for child in element :
            if( child.tag == energy.moniker ) :
                _energy = energy.parseXMLNode( child, xPath, linkData )
            elif( child.tag == energyAngular.moniker ) :
                _energyAngular = energyAngular.parseXMLNode( child, xPath, linkData )
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        energyAngularMC = form( element.get( "label" ), element.get( "productFrame" ), _energy, _energyAngular )

        xPath.pop( )
        return( energyAngularMC )
