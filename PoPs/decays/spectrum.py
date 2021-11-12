# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes for storing an outgoing energy spectrum for decay products.
The spectrum may include a combination of discrete and continuum contributions.

Spectral data are organized as follows::

   - decayMode
      - ...
      - spectra
         - spectrum (1 or more)
            - discrete (0 or more)
               - intensity
               - energy
               - ...
            - continuum (0 or more)

Each spectrum corresponds to one type of decay product (e.g. photon, positron or neutron).
Sometimes a product is broken out into more than one spectrum sections, however. For example, photons
may be divided into 'gamma' and 'x-ray' parts.

Discrete radiations each have an associated intensity, indicating how often that discrete radiation is emitted
per decay. The intensity may be > 1, for example for 511 keV photons emitted by beta+ decay.
Discrete radiation may be accompanied by internalConversionCoefficients, internalPairFormationCoefficients, etc.

Continuum radiation is described by a probability distribution of outgoing energy.
    FIXME is it really a probability distribution? Not always normalized...
"""

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
    """Stores a decay particle emitted with discrete energy"""

    moniker = 'discrete'

    def __init__( self, _intensity, _energy, type = None, _internalPairFormationCoefficient = None,
                  _positronEmissionIntensity = None ) :
        """
        :param _intensity: relative intensity for emission of this discrete particle
        :param _energy: outgoing product energy
        :param type: optional type of transition (i.e. 'allowed', 'first-forbidden', etc.
        :param _internalPairFormationCoefficient: optional coefficient for forming an electron/positron pair
            Should only be used with electromagnetic decays
        :param _positronEmissionIntensity: optional coefficient for emitting a positron
            Should only be used with electron capture / beta+ decays
        """

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
    """Describes continuum particle emission"""

    moniker = 'continuum'

    def __init__( self, spectrum ) :
        """
        :param spectrum: outgoing energy probability distribution
        """

        ancestryModule.ancestry.__init__( self )
        self.__spectrum = spectrum

    @property
    def spectrum( self ) :

        return( self.__spectrum )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        self.__spectrum.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

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
    """
    Contains the outgoing spectrum for one type of decay product.
    """

    moniker = 'spectrum'

    def __init__( self, label, pid ) :
        """
        :param label: unique label for the outgoing product (i.e. 'x-ray' or 'gamma')
        :param pid: PoPs particle id (i.e. 'photon')
        """

        miscModule.classWithLabelKey.__init__( self, label )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid not str' )
        self.__pid = pid

        self.__emissions = []

    def __len__( self ) :

        return( len( self.__emissions ) )

    def __getitem__( self, index ) :

        return( self.__emissions[index] )

    @property
    def emissions( self ) :
        """Returns the list of all discrete and continuum emissions"""

        return( self.__emissions )      # FIXME accesses mutable member

    @property
    def pid( self ) :

        return( self.__pid )

    def append( self, emission ) :
        """
        Add new emitted particle to self
        :param emission: 'discrete' or 'continuum' instance
        """

        if( not( isinstance( emission, ( discrete, continuum ) ) ) ) : raise TypeError( 'emission must be instance of discrete or continuum' )
        self.__emissions.append( emission )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        for emission in self.__emissions : emission.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

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
    """Contains the outgoing spectrum list for a decayMode"""

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
