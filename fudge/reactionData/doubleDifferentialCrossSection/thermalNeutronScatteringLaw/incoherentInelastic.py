# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Thermal neutron scattering law coherent elastic double difference cross section form and its supporting classes.
"""

__metaclass__ = type

import math

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import physicalQuantity as physicalQuantityModule
from xData import XYs as XYs1dModule
from xData import gridded as griddedModule

from PoPs import IDs as IDsPoPsModule

from fudge import suites as suitesModule

from . import base as baseModule
from . import incoherentInelasticMisc as incoherentInelasticMiscModule

class options( ancestryModule.ancestry ) :
    """
    Stores options calculatedAtThermal and asymmetric for inelastic incoherent scattering.
    """

    moniker = 'options'

    def __init__( self, calculatedAtThermal = False, asymmetric = False ) :

        ancestryModule.ancestry.__init__( self )
        self.calculatedAtThermal = calculatedAtThermal
        self.asymmetric = asymmetric

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        pass

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        calculatedAtThermal = 'false'
        if( self.calculatedAtThermal ) : calculatedAtThermal = 'true'

        asymmetric = 'false'
        if( self.asymmetric ) : asymmetric = 'true'

        xmlStringList = [ '%s<%s calculatedAtThermal="%s" asymmetric="%s"/>' % ( indent, self.moniker, calculatedAtThermal, asymmetric ) ]

        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )

        kwargs = {  'calculatedAtThermal'   : element.get( 'calculatedAtThermal' ) == 'true',
                    'asymmetric'            : element.get( 'asymmetric' ) == 'true' }
        _options = options( **kwargs )

        xPath.pop( )
        return( _options )

class mass( physicalQuantityModule.physicalQuantity ) :

    moniker = 'mass'

class freeAtomCrossSection( physicalQuantityModule.physicalQuantity ) :

    moniker = 'freeAtomCrossSection'

class e_critical( physicalQuantityModule.physicalQuantity ) :

    moniker = 'e_critical'

class e_max( physicalQuantityModule.physicalQuantity ) :

    moniker = 'e_max'

class T_effective( ancestryModule.ancestry ) :
    """
    In incoherent inelastic sections, each scattering atom using the short collision time (SCT)
    approximation needs an effective temperature.
    """

    moniker = 'T_effective'

    def __init__( self, function1d ) :

        ancestryModule.ancestry.__init__( self )

        self.__function1d = function1d
        self.__function1d.setAncestor( self )

    @property
    def function1d( self ) :

        return( self.__function1d )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__function1d.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.__function1d.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
    
        xPath.append( element.tag )

        xys1d = XYs1dModule.XYs1d.parseXMLNode( element[0], xPath, linkData )
        T_effective1 = T_effective( xys1d )

        xPath.pop( )
        return( T_effective1 )

class scatteringAtom( ancestryModule.ancestry ):
    """
    Inelastic incoherent scattering requires a description of each atom that is scattered off.
    """

    moniker = 'scatteringAtom'

    def __init__(self, label, numberPerMolecule, mass, freeAtomCrossSection, e_critical = None, e_max = None,
            functionalForm = None, T_effective = None ):

        ancestryModule.ancestry.__init__( self )
        self.label = label
        self.numberPerMolecule = numberPerMolecule
        self.functionalForm = functionalForm

        self.mass = mass
        self.mass.setAncestor( self )

        self.freeAtomCrossSection = freeAtomCrossSection
        self.freeAtomCrossSection.setAncestor( self )

        self.__boundCrossSection = None

        self.e_critical = e_critical
        if( self.e_critical is not None ) : self.e_critical.setAncestor( self )

        self.e_max = e_max
        if( self.e_max is not None ) : self.e_max.setAncestor( self )

        self.T_effective = T_effective
        if( self.T_effective is not None ) : self.T_effective.setAncestor( self )

    def boundCrossSection( self ) :

        if self.__boundCrossSection is None:
            neutronMassAMU = self.findAttributeInAncestry('PoPs')['n'].mass.float('amu')
            self.__boundCrossSection = self.freeAtomCrossSection.value * ((self.mass.value + neutronMassAMU) / self.mass.value)**2
        return self.__boundCrossSection

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.mass.convertUnits( unitMap )
        self.freeAtomCrossSection.convertUnits( unitMap )
        if( self.e_critical is not None ) : self.e_critical.convertUnits( unitMap )
        if( self.e_max is not None ) : self.e_max.convertUnits( unitMap )
        if( self.T_effective is not None ) : self.T_effective.convertUnits( unitMap )

    def shortCollisionTime( self, alpha, beta, T ):
        """
        Equation 7.8 in ENDF manual
        :param alpha:
        :param beta:
        :param T:
        :return:
        """

        Teff = self.T_effective.function1d.evaluate( T )
            # num = math.exp( -( alpha - abs( beta ) )**2 * T / ( 4 * alpha * Teff ) - ( beta + abs( beta ) ) / 2 )
            # reformulate using T/Teff = 1-f, should be more numerically stable
        f = 1 - T / Teff
        num = math.exp( ( f * ( alpha - abs( beta ) )**2 - ( alpha + beta )**2 ) / ( 4 * alpha ) )
        den = math.sqrt( 4 * math.pi * alpha * Teff / T )
        return( num / den )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ''
        if( self.functionalForm is not None ) : attrs = ' functionalForm="%s"' % self.functionalForm
        xmlStringList = [ '%s<%s label="%s" numberPerMolecule="%d"%s>' % ( indent, self.moniker, self.label, self.numberPerMolecule, attrs ) ]
        for attribute in ( 'mass', 'freeAtomCrossSection', 'e_critical', 'e_max', 'T_effective' ) :
            if( getattr( self, attribute ) is not None ) : xmlStringList += getattr( self, attribute ).toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        kwargs = {  'label' : element.get( 'label' ),
                    'numberPerMolecule' : int( element.get( 'numberPerMolecule' ) ),
                    'functionalForm' : element.get( 'functionalForm' ) }
        for child in element:
            _class = {  mass.moniker : mass,
                        freeAtomCrossSection.moniker : freeAtomCrossSection,
                        e_critical.moniker : e_critical,
                        e_max.moniker : e_max,
                        T_effective.moniker : T_effective }.get( child.tag, None )
            if( _class is None ) : print( "Warning: encountered unexpected element '%s' in scatteringAtom!" % child.tag )
            kwargs[child.tag] = _class.parseXMLNode( child, xPath, linkData )
        SA = scatteringAtom( **kwargs )
        xPath.pop( )
        return( SA )

class scatteringAtoms( suitesModule.suite ) :

    moniker = 'scatteringAtoms'

    def __init__( self ) :

        suitesModule.suite.__init__( self, ( scatteringAtom, ) )

class S_alpha_beta( ancestryModule.ancestry ):
    """
    Inelastic incoherent section contains S as a function of (alpha,beta,T). Currently supports gridded3d.
    """

    moniker = 'S_alpha_beta'

    def __init__( self, gridded3d ) :

        ancestryModule.ancestry.__init__( self )
        self.gridded3d = gridded3d

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.gridded3d.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.gridded3d.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        gridded3d = griddedModule.gridded3d.parseXMLNode( element[0], xPath, linkData )
        Sab = S_alpha_beta( gridded3d )

        xPath.pop( )
        return( Sab )

class form( baseModule.form ) :

    moniker = 'thermalNeutronScatteringLaw_incoherentInelastic'
    process = 'thermalNeutronScatteringLaw incoherent-inelastic'
    subformAttributes = ( 'options', 'scatteringAtoms', 'S_alpha_beta' )

    def __init__( self, label, _options, _S_alpha_beta, atoms = () ) :

        baseModule.form.__init__( self, IDsPoPsModule.neutron, label, standardsModule.frames.labToken, ( _options, scatteringAtoms( ), _S_alpha_beta ) )

        if( not( isinstance( _options, options ) ) ) : raise TypeError( "Invalid options for %s." % self.moniker )
        self.options = _options

        if( not( isinstance( _S_alpha_beta, S_alpha_beta ) ) ) : raise TypeError( "Invalid S_table for %s." % self.moniker )
        self.S_alpha_beta = _S_alpha_beta

        for atom in atoms : self.scatteringAtoms.add( atom )

    @property
    def domainUnit( self ) :

        return( str( self.scatteringAtoms[0].e_max.unit ) )

    def energyMaximum( self ) :

        return( float( self.scatteringAtoms[0].e_max ) )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        from fudge import styles as stylesModule

        temperature = style.temperature.getValueAs( kwargs['temperatureUnit'] )

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )

        energyDomain = style.findDerivedFromStyle(stylesModule.evaluated).projectileEnergyDomain
        energyMax = PQUModule.PQU( energyDomain.max, energyDomain.unit ).getValueAs( kwargs['incidentEnergyUnit'] )

        return( incoherentInelasticMiscModule.process( self, style, energyMin, energyMax, temperature, kwargs ) )

    def temperatures( self, _temperatures ) :

        _temperatures['incoherent-inelastic'] = [ self.S_alpha_beta.gridded3d.axes[3].unit, self.S_alpha_beta.gridded3d.axes[3].values.values ]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )

        label = element.get( 'label' )
        _options = options.parseXMLNode( element.find( options.moniker ), xPath, linkData )
        Sab = S_alpha_beta.parseXMLNode( element.find( S_alpha_beta.moniker ), xPath, linkData )
        atoms = [ scatteringAtom.parseXMLNode( child, xPath, linkData ) for child in element.find( 'scatteringAtoms' ) ]
        incoherentInelastic = form( label, _options, Sab, atoms = atoms )

        xPath.pop( )
        return( incoherentInelastic )
