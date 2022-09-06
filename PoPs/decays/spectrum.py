# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

from LUPY import ancestry as ancestryModule
from LUPY import enums as enumsModule

from xData import axes as axesModule
from xData import physicalQuantity as physicalQuantityModule
from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule
from xData import XYs1d as XYs1dModule

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

class TransitionType(enumsModule.Enum):

    none = ''
    allowed = 'allowed'
    firstForbidden = 'first-forbidden'
    secondForbidden = 'second-forbidden'

    @staticmethod
    def ENDF_index(transitionType):

        if not isinstance(transitionType, TransitionType):
            raise TypeError('transitionType must be a TransitionType: got "%s".' % type(transitionType))

        if transitionType == TransitionType.none:
            return 0
        if transitionType == TransitionType.allowed:
            return 1
        if transitionType == TransitionType.firstForbidden:
            return 2
        if transitionType == TransitionType.secondForbidden:
            return 3

        raise TypeError('This should never happen.')

#
# FIXME Need a PhysicalQuantity class that has keyName.
#
class Unitless( physicalQuantityModule.PhysicalQuantity ) :

    keyName = 'label'

    def __init__( self, value, label = None ) :

        physicalQuantityModule.PhysicalQuantity.__init__( self, value, '', label = label )

    def copy( self ) :              # overrides required since this __init__ takes different arguments:

        cls = self.__class__( self.value, self.label )
        if( self.uncertainty is not None ) : cls.uncertainty = self.uncertainty.copy( )
        return( cls )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        value = element.get( 'value' )
        label = element.get( 'label', None )
        _cls = cls( value, label )
        for child in element :
            if( child.tag == uncertaintyModule.Uncertainty.moniker ) :
                _cls.uncertainty = uncertaintyModule.Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        xPath.pop( )
        return( _cls )

class Intensity( Unitless ) :

    moniker = 'intensity'

class Energy( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'energy'

class Shell(physicalQuantityModule.PhysicalQuantity):

    moniker = 'shell'
    total = 'total'
    KShell = 'K'
    LShell = 'L'

class InternalConversionCoefficients( suiteModule.Suite ) :

    moniker = 'internalConversionCoefficients'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( Shell, ) )

class PhotonEmissionProbabilities( suiteModule.Suite ) :

    moniker = 'photonEmissionProbabilities'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( Shell, ) )

class InternalPairFormationCoefficient( Unitless ) :

    moniker = 'internalPairFormationCoefficient'

class PositronEmissionIntensity( Unitless ) :

    moniker = 'positronEmissionIntensity'

class Discrete(ancestryModule.AncestryIO):
    """Stores a decay particle emitted with discrete energy"""

    moniker = 'discrete'

    def __init__( self, _intensity, _energy, type = TransitionType.none, _internalPairFormationCoefficient = None,
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

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( _intensity, Intensity ) ) ) : raise TypeError( '_intensity must be an instance of intensity' )
        self.__intensity = _intensity

        if( not( isinstance( _energy, Energy ) ) ) : raise TypeError( '_energy must be an instance of Energy' )
        self.__energy = _energy

        self.__type = TransitionType.checkEnumOrString(type)

        self.__internalConversionCoefficients = InternalConversionCoefficients( )
        self.__internalConversionCoefficients.setAncestor( self )

        if( _internalPairFormationCoefficient is not None ) :
            if( not( isinstance( _internalPairFormationCoefficient, InternalPairFormationCoefficient ) ) ) :
                raise TypeError( '_internalPairFormationCoefficient must be an instance of InternalPairFormationCoefficient' )
        self.__internalPairFormationCoefficient = _internalPairFormationCoefficient

        if( _positronEmissionIntensity is not None ) :
            if( not( isinstance( _positronEmissionIntensity, PositronEmissionIntensity ) ) ) :
                raise TypeError( '_positronEmissionIntensity must be an instance of PositronEmissionIntensity' )
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        type = ''
        if( self.type is not TransitionType.none ) : type = ' type="%s"' % self.type
        XMLStringList = [ '%s<%s%s>' % ( indent, self.moniker, type ) ]
        XMLStringList += self.__intensity.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__energy.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__internalConversionCoefficients.toXML_strList( indent2, **kwargs )
        if( self.__internalPairFormationCoefficient is not None ) : XMLStringList += self.__internalPairFormationCoefficient.toXML_strList( indent2, **kwargs )
        if( self.__positronEmissionIntensity is not None): XMLStringList += self.__positronEmissionIntensity.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        kwargs = {'type': element.get('type', TransitionType.none)}
        ICCelement = None
        for child in element :
            found = False
            for subclass in ( Intensity, Energy, InternalPairFormationCoefficient, PositronEmissionIntensity ) :
                if( child.tag == subclass.moniker ) :
                    kwargs['_'+subclass.moniker] = subclass.parseNodeUsingClass(element.find(subclass.moniker), xPath, linkData, **kwargs)
                    found = True
            if not found:
                if child.tag == InternalConversionCoefficients.moniker:
                    ICCelement = child
                else:
                    raise ValueError("Encountered unexpected child element '%s'" % child.tag)

        _discrete = cls( **kwargs )

        if ICCelement is not None:
            _discrete.internalConversionCoefficients.parseNode(ICCelement, xPath, linkData, **kwargs)

        xPath.pop()

        return _discrete

class Continuum(ancestryModule.AncestryIO):
    """Describes continuum particle emission"""

    moniker = 'continuum'

    def __init__( self, spectrum ) :
        """
        :param spectrum: outgoing energy probability distribution
        """

        ancestryModule.AncestryIO.__init__(self)
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__spectrum.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        _spectrum = XYs1dModule.XYs1d.parseNodeUsingClass(element.find(XYs1dModule.XYs1d.moniker), xPath, linkData, **kwargs)
        instance = Continuum( _spectrum )
        xPath.pop()
        return instance

    @staticmethod
    def defaultAxes( energyUnit ) :

        axes = axesModule.Axes(2)
        axes[0] = axesModule.Axis( "P(energy_out)", 0, "1/%s" % energyUnit )
        axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
        return( axes )

class Spectrum( miscModule.ClassWithLabelKey ) :
    """
    Contains the outgoing spectrum for one type of decay product.
    """

    moniker = 'spectrum'

    def __init__( self, label, pid ) :
        """
        :param label: unique label for the outgoing product (i.e. 'x-ray' or 'gamma')
        :param pid: PoPs particle id (i.e. 'photon')
        """

        miscModule.ClassWithLabelKey.__init__( self, label )

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

        if( not( isinstance( emission, ( Discrete, Continuum ) ) ) ) : raise TypeError( 'emission must be instance of Discrete or Continuum' )
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" pid="%s">' % ( indent, self.moniker, self.label, self.pid ) ]
        for _emission in self.__emissions : XMLStringList += _emission.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )

        for child in element:
            if( child.tag == Discrete.moniker ) :
                self.append(Discrete.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            elif( child.tag == Continuum.moniker ) :
                self.append(Continuum.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            else :
                raise ValueError("Encountered unexpected child element '%s'" % child.tag)

        xPath.pop( )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )
        _spectrum = cls( element.get('label'), element.get('pid') )
        xPath.pop()

        _spectrum.parseNode(element, xPath, linkData, **kwargs)

        return _spectrum

class Spectra( suiteModule.Suite ) :
    """Contains the outgoing spectrum list for a decayMode"""

    moniker = 'spectra'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( Spectrum, ) )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

#
# FIXME, this method is broken. What is product?
#
        for child in element :
            _spectrum = Spectrum.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            self.add( _spectrum )

        xPath.pop( )
        return( self )
