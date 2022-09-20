# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Energy/angular double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import energy as energyModule


class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( Xs_pdf_cdf1d, ) )

class XYs3d( multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class Subform( baseModule.Subform ) :
    """Abstract base class for energyAngular subforms."""

    moniker = 'dummy'                   # This is not used but added to stop lylint from reporting.

    def __init__( self, data ) :

        baseModule.Subform.__init__( self )
        if( not( isinstance( data, self.allowedSubElements ) ) ) : raise TypeError( 'Invalid instance: %s' % type( data ) )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.data.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class Energy( Subform ) :

    moniker = 'energy'
    allowedSubElements = ( energyModule.XYs2d, )
    ancestryMembers = ( 'energy', )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {energyModule.XYs2d.moniker: energyModule.XYs2d}.get(subformElement.tag)
        if( subformClass is None ) : raise Exception( 'unknown energy subform "%s"' % subformElement.tag )
        energySubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _energy = cls(energySubform)

        xPath.pop( )
        return( _energy )

class EnergyAngular( Subform ) :

    moniker = 'energyAngular'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'energyAngular', )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {XYs3d.moniker: XYs3d}.get(subformElement.tag)
        if( subformClass is None ) : raise Exception( 'unknown energyAngular subform "%s"' % subformElement.tag )
        energyAngularSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _energyAngular = cls(energyAngularSubform)

        xPath.pop( )
        return( _energyAngular )

class Form( baseModule.Form ) :

    moniker = 'energyAngularMC'
    subformAttributes = ( 'energy', 'energyAngular' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _energy, _energyAngular ) :

        if( not( isinstance( _energy, Energy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energy ) )
        if( not( isinstance( _energyAngular, EnergyAngular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energyAngular ) )
        baseModule.Form.__init__( self, label, productFrame, ( _energy, _energyAngular ) )

    @property
    def domainMin( self ) :

        return( self.energy.data.domainMin )

    @property
    def domainMax( self ) :

        return( self.energy.data.domainMax )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Translate a GNDS energyAngularMC node."""

        xPath.append(node.tag)

        _energy = None
        _energyAngular = None
        for child in node:
            if( child.tag == Energy.moniker ) :
                _energy = Energy.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == EnergyAngular.moniker ) :
                _energyAngular = EnergyAngular.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        energyAngularMC = cls(node.get("label"), node.get("productFrame"), _energy, _energyAngular)

        xPath.pop( )
        return( energyAngularMC )
