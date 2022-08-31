# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Angular/Energy double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import angular as angularModule


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
    """Abstract base class for angularEnergyMC subforms."""

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

class Angular( Subform ) :

    moniker = 'angular'
    allowedSubElements = ( angularModule.XYs2d, )
    ancestryMembers = ( 'angular', )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {angularModule.XYs2d.moniker: angularModule.XYs2d}.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _angular= cls(angularSubform)

        xPath.pop( )
        return( _angular )

class AngularEnergy( Subform ) :

    moniker = 'angularEnergy'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'angularEnergy', )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {XYs3d.moniker: XYs3d}.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angularEnergy subform "%s"' % subformElement.tag )
        angularEnergySubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _angularEnergy = cls(angularEnergySubform)

        xPath.pop( )
        return( _angularEnergy )

class Form( baseModule.Form ) :

    moniker = 'angularEnergyMC'
    subformAttributes = ( 'angular', 'angularEnergy' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angular, _angularEnergy ) :

        if( not( isinstance( _angular, Angular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angular ) )
        if( not( isinstance( _angularEnergy, AngularEnergy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angularEnergy ) )
        baseModule.Form.__init__( self, label, productFrame, ( _angular, _angularEnergy ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angular.convertUnits( unitMap )
        self.angularEnergy.convertUnits( unitMap )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Translate a GNDS angularEnergyMC node."""

        xPath.append(node.tag)

        _angular = None
        _angularEnergy = None
        for child in node:
            if( child.tag == Angular.moniker ) :
                _angular = Angular.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == AngularEnergy.moniker ) :
                _angularEnergy = AngularEnergy.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        angularEnergyMC = cls(node.get("label"), node.get('productFrame'), _angular, _angularEnergy)

        xPath.pop( )
        return angularEnergyMC
