# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
Angular distribution classes with forms :py:class:`form`, :py:class:`TwoBody` and :py:class:`CoulombElasticForm`.

This details most of the stuff in :py:mod:`angular.py`, and it also applies to :py:mod:`energy.py`,
:py:mod:`energyAngular.py` and :py:mod:`angularEnergy.py`. Here we only describe some of the old classes.
Namely,

    * ``pdfOfMu.pointwise``

        - Stores a :math:`P(\mu)` as a list of :math:`\mu_i`, :math:`P_i`

    * ``pdfOfMu.Legendre``

        - Stores a :math:`P(\mu)`  as a Legendre series

    * ``pointwise``

        - Stores the :math:`P(\mu|E)` as a list of :math:`P(\mu)` which can be either
          pdfOfMu.pointwise of pdfOfMu.Legendre

    * ``piecewise``

        - Stores :math:`P(\mu|E)` as a contigous list of ``pointwise`` instances.

and their new equivalent classes

    * :py:class:`XYs1d`

        - Same as ``pdfOfMu.pointwise``

    * :py:class:`Legendre`

        - Same as ``pdfOfMu.Legendre``

    * :py:class:`XYs2d`

        - Same as ``pointwise``

    * :py:class:`Regions2d`

        - Same as ``piecewise``


Note, there was no ``pdfOfMu.piecewise`` and its new equivalent :py:class:`Regions1d`, as
they are not currently needed. There are other classes in :py:mod:`angular.py`
(e.g., :py:class:`Isotropic2d`, :py:class:`Forward`, :py:class:`Recoil`) but we can ignore them for now.

So what are the classes in :py:mod:`angular.py` for? They are to represent the
physics function :math:`P(\mu|E)` where :math:`\mu` is the outgoing :math:`\cos( \theta )` and :math:`E` is
the projectile's energy. Currently, we only store these as list of :math:`P(\mu)`
or :math:`C_l`'s (or regions of these).

For the simple case with only :math:`P(\mu`'s (no :math:`C_l`'s), and no change of
interpolation or step in function, a :math:`P(\mu|E)` looks like::

    E_1, P_1(mu)
    E_2, P_2(mu)
    ...
    E_n-1, P_n-1(mu)
    E_1, P_1(mu)

In the old code each :math:`P_i(\mu)` was stored as an instance of a
``pdfOfMu.pointwise`` and the P(mu|E) was an instance of ``pointwise``. In the new
code each :math:'P_i(\mu)` is stored as an instance of :py:class:`XYs1d`, and :math:`P(\mu|E)` is an
instance of :py:class:`XYs2d`. There are several nice things about the new structure.
For example, the names here match the names in ``xData``. Also, the internal
functions do not need to be nested as, for example, ``pdfOfMu.pointwise`` was.
This latter point is a result of adding the dimension to the name of each
class (e.g., :py:class:`XYs1d` vs. :py:class:`XYs2d`) so there is no name clashing.


As a further note, in MF=4, the LTT = 1, 2 and 3 are stored using

LTT=1
    XYs1d for :math:`P(\mu)` and :py:class:`XYs2d` for the list of :math:`{E_i,P(\mu)}`.

LTT=2
    Legendre for :math:`P(\mu)` and :py:class:`XYs2d` for the list of :math:`{E_i,P(\mu)}`.

LTT=3
    :py:class:`Regions2d` that stores the LTT=1 in the first region and LTT=2 in the second region.

"""

import math
from fudge.core.math import fudgemath

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from numericalFunctions import Legendre as LegendreModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import series1d as series1dModule
from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule
from xData import link as linkModule

from .. import averageProductEnergy as averageProductEnergyModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule

from . import base as baseModule
from . import reference as referenceModule
from . import miscellaneous as miscellaneousModule
from . import probabilities as probabilitiesModule
from . import energy as energyModule


def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(3)
    axes[0] = axesModule.Axis( 'P(mu|energy_in)', 0, '' )
    axes[1] = axesModule.Axis( 'mu', 1, '' )
    axes[2] = axesModule.Axis( 'energy_in', 2, energyUnit )
    return( axes )

class XYs1d( XYs1dModule.XYs1d ) :

    def averageMu( self ) :

        allowedInterpolations = [xDataEnumsModule.Interpolation.linlin, xDataEnumsModule.Interpolation.flat]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYs1dModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

    def invert( self ) :

        reflected = self.copy( )
        pdf = [ [ -x, y ] for x, y in self ]
        pdf.reverse( )
        reflected.setData( pdf )
        return( reflected )

    def isIsotropic( self ) :

        return( self.rangeMin == self.rangeMax )

    def LegendreCoefficient( self, LegendreOrder ) :

        return( LegendreModule.from_pointwiseXY_C( self, LegendreOrder )[LegendreOrder] )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    def to_xs_pdf_cdf1d(self, style, tempInfo, indent):

        return Xs_pdf_cdf1d.fromXYs(self, outerDomainValue=self.outerDomainValue, thinEpsilon=1e-14)

class Legendre( series1dModule.LegendreSeries ) :

    def averageMu( self ) :

        return( self.getCoefficientSafely( 1 ) )

    def LegendreCoefficient( self, LegendreOrder ) :

        if( LegendreOrder < len( self ) ) : return( self[LegendreOrder] )
        return( 0.0 )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( accuracy = XYs1dModule.defaultAccuracy, lowerEps = 0, upperEps = 1e-8 )
        linear.outerDomainValue = self.outerDomainValue
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :

    def isIsotropic( self ) :

        return( min( self.pdf.values ) == max( self.pdf.values ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class Isotropic1d(ancestryModule.AncestryIO):

    moniker = 'isotropic1d'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

    @property
    def domainMin( self ) :

        return( -1.0 )

    @property
    def domainMax( self ) :

        return( 1.0 )

    def evaluate( self, mu ) :

        return( 0.5 )

    def isIsotropic( self ) :

        return( True )

    def integrate(self, domainMin, domainMax):

        return 0.5 * (domainMax - domainMin)

    def LegendreCoefficient( self, LegendreOrder ) :

        if( LegendreOrder == 0 ) : return( 1.0 )
        return( 0.0 )

    def toXML_strList(self, **kwargs):

        raise Exception('Not supported')

    @classmethod
    def parseNodeUsingClass(self, **kwargs):

        raise Exception('Not supported')

class Subform( baseModule.Subform ) :
    """Abstract base class for angular subforms."""

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( None )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        The method does nothing.
        """

        return 0

class Isotropic2d( Subform ) :

    moniker = 'isotropic2d'

    def __init__( self ) :

        Subform.__init__( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( Isotropic2d( ) )

    __copy__ = copy

    @property
    def domainGrid( self ) :

        return( [ self.domainMin, self.domainMax ] )

    @property
    def domainMin( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMin )

    @property
    def domainMax( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMax )

    @property
    def domainUnit( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainUnit )

    def evaluate( self, energy ) :

        _isotropic1d = Isotropic1d( )
        _isotropic1d.setAncestor( self )
        return( _isotropic1d )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def invert( self ) :

        return( self.copy( ) )

    def isIsotropic( self ) :

        return( True )

    def averageMu( self, E, accuracy = None ) :

        return( 0. )

    def check( self, info ) :

        return []

    def integrate( self, energyIn, muOut ) :

        if( self.domainMin <= energyIn <= self.domainMax ) :
            domainMin, domainMax = miscellaneousModule.domainLimits( muOut, -1.0, 1.0 )
            if( domainMax is None ) : return( 1.0 )
            if( domainMax == 0.0 ) : return( 0.5 )
            return( 0.5 * ( domainMax - domainMin ) )

        return( 0.0 )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        ptw = XYs2d( axes = defaultAxes( self.domainUnit ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], outerDomainValue = self.domainMin ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], outerDomainValue = self.domainMax ) )
        return( ptw )

    def toXML_strList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        return cls()

class Forward( Subform ) :

    moniker = 'forward'

    def __init__( self ) :

        Subform.__init__( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( Forward( ) )

    __copy__ = copy

    @property
    def domainMin( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMin )

    @property
    def domainMax( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMax )

    @property
    def domainUnit( self ) :

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainUnit )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def invert( self ) :

        return( self.toPointwise_withLinearXYs( ).invert( ) )

    def isIsotropic( self ) :

        return( False )

    def averageMu( self, E, accuracy = None ) :

        return( 1. )

    def check( self, info ) :

        return []

    def toPointwise_withLinearXYs( self, **kwargs ) :

        accuracy = xDataBaseModule.getArguments( kwargs, { 'accuracy' : 1e-6 } )['accuracy']

        ptw = XYs2d(axes=defaultAxes(self.domainUnit))
        ptw.append( XYs1dModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMin ) )
        ptw.append( XYs1dModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMax ) )
        return( ptw )

    def toXML_strList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        return cls()

class Recoil( linkModule.Link, Subform ) :

    moniker = 'recoil'

    def __init__( self, link, root = None, path = None, label = None, relative = False ) :
        """
        Distributions for this particle can be computed from its recoil partner using kinematics.
        Only meant for angular 2-body reactions.
        The 'root' and 'path' usually do not need to be specified: they will be computed from the link.

        :param link: distribution form for recoil partner
        :type link: distribution.base.Form
        :param root: filename where recoil partner is defined (usually None)
        :param path: xpath to recoil partner (usually will be computed from the link)
        :param label: label for this subform (usually None)
        :param relative: whether to write xpath as relative or absolute
        :return:
        """

        linkModule.Link.__init__( self, link=link, root=root, path=path, label=label, relative=relative )
        Subform.__init__( self )

    @property
    def domainUnit( self ) :

        return( self.partner.angularSubform.domainUnit )

    @property
    def domainMin(self):
        """Returns the domainMin."""

        return self.partner.angularSubform.domainMin

    @property
    def domainMax(self):
        """Returns the domainMax."""

        return self.partner.angularSubform.domainMax

    @property
    def partner( self ) :

        if( self.link is None ) : raise Exception( "Encountered unresolved link!" )
        return( self.link )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( Recoil( self.partner ) )

    __copy__ = copy

    def evaluate( self, energy ) :

        return( self.partner.angularSubform.evaluate( energy ).invert( ) )

    def getNumericalDistribution( self ) :

        partnerForm = self.partner.toPointwise_withLinearXYs( upperEps = 1e-8 ).invert( )
        return( partnerForm )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.partner.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def isIsotropic( self ) :

        return( self.partner.isIsotropic( ) )

    def averageMu( self, E, accuracy = None ) :

        return( -self.partner.averageMu( E, accuracy = accuracy ) )

    def check( self, info ) :

        from fudge import warning

        warnings = []
        if( not( isinstance( self.partner, ( TwoBody, referenceModule.CoulombPlusNuclearElastic ) ) ) ) :
            warnings.append( warning.MissingRecoilDistribution( self ) )

        return warnings

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.partner.invert( ).toPointwise_withLinearXYs( **kwargs ) )

class XYs2d( Subform, probabilitiesModule.PofX1GivenX2 ) :

    def __init__( self, **kwargs ) :

        probabilitiesModule.PofX1GivenX2.__init__( self, **kwargs )
        Subform.__init__( self )

    def getAtEnergy( self, energy ) :

        return self.interpolateAtValue(energy, extrapolation=xDataEnumsModule.Extrapolation.flat)

    def fixDomains(self, domainMin, domainMax, fixToDomain, tweakLower = True):
        """
        This method call **fixDomains* on the *self* via the **probabilitiesModule.PofX1GivenX2** class.
        """

        return probabilitiesModule.PofX1GivenX2.fixDomains(self, domainMin, domainMax, fixToDomain, tweakLower = tweakLower)

    def invert( self ) :

        ptw = XYs2d( axes = self.axes.copy( ) )
        for POfMu in self : ptw.append( POfMu.invert( ) )
        return( ptw )

    def isIsotropic( self ) :

        for energy_in in self :
            if( not( energy_in.isIsotropic( ) ) ) : return( False )
        return( True )

    def averageMu( self, E, accuracy = None ) :

        mode, function1, function2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( E )
        if( mode is None ) : raise ValueError( 'angular has no data' )
        mu1 = function1.normalize( ).averageMu( )
        if( mode == '' ) :
            mu2 = function2.normalize( ).averageMu( )
            mu1 = ( 1 - frac ) * mu1 + frac * mu2
            if   interpolationQualifier == xDataEnumsModule.InterpolationQualifier.unitBase:
                raise TypeError( 'Unsupported interpolation qualifier "%s"' % interpolationQualifier )
            elif interpolationQualifier != xDataEnumsModule.InterpolationQualifier.none:
                raise TypeError( 'Unsupported interpolation qualifier "%s"' % interpolationQualifier )
        return( mu1 )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.outerDomainValue for data in self ]
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def check( self, info ) :

        from fudge import warning

        warnings = []
        for idx,function in enumerate(self):
            xys = function.toPointwise_withLinearXYs( accuracy = info['normTolerance'], lowerEps = 1e-6 )

            if isinstance(function, Legendre):
                integral = function.coefficients[0]
            elif isinstance(function, XYs1d):
                integral = xys.integrate(-1,1)
            else:
                raise NotImplementedError("checking function of type %s" % type(function))  # FIXME support xs_pdf_cdf1d

            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.UnnormalizedDistribution( PQUModule.PQU( function.outerDomainValue, self.axes[-1].unit ), idx, integral, function ) )

            if xys.rangeMin < 0.0:
                warnings.append(warning.NegativeProbability(xys.rangeMin, PQUModule.PQU(function.outerDomainValue, self.axes[-1].unit), obj=function))

        return warnings

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        subform = XYs2d( axes = self.axes.copy( ), interpolation = self.interpolation )
        for xys in self : subform.append( xys.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( subform )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, Legendre, Xs_pdf_cdf1d ) )

class Regions2d( Subform, regionsModule.Regions2d ) :

    def __init__( self, **kwargs ) :

        regionsModule.Regions2d.__init__( self, **kwargs )
        Subform.__init__( self )

    def check( self, info ):

        from fudge import warning

        warnings = []
        for idx,region in enumerate(self):
            regionWarnings = region.check( info )
            if regionWarnings: warnings.append( warning.Context('region index %i: %s'
                % (idx, region.moniker), regionWarnings ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        EMax = kwargs['EMax']
        aveEnergies = []
        aveMomenta = []
        for i1, region in enumerate( self ) :
            kwargs['EMax'] = region.domainMax
            if( i1 == ( len( self ) - 1 ) ) : kwargs['EMax'] = EMax
            aveEnergy, aveMomentum = calculateAverageProductData( region, style, indent, **kwargs ) 
            kwargs['EMin'] = kwargs['EMax']
            aveEnergies.append( aveEnergy[0] )
            aveMomenta.append( aveMomentum[0] )
        return( aveEnergies, aveMomenta )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        This method call **fixDomains* on the *self* via the **regionsModule.Regions2d** class.
        """

        return regionsModule.Regions2d.fixDomains(self, domainMin, domainMax, fixToDomain)

    def isIsotropic( self ) :

        for region in self :
            if( not( region.isIsotropic( ) ) ) : return( False )
        return( True )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ self[0][0].outerDomainValue ]
        for region in self: Es.extend( [ data.outerDomainValue for data in region[1:] ] )
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def averageMu( self, E, accuracy = None ) :

        # FIXME what if E lies on boundary between regions?
        for region in self:
            if( region.domainMax > E ) :
                return region.averageMu( E, accuracy )
        return self[-1].averageMu( E, accuracy )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        _regions2d = Regions2d( axes = self.axes )
        for region in self : _regions2d.append( region.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( _regions2d )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        arguments = self.getArguments( kwargs, { 'lowerEps' : 0, 'upperEps' : 0 } )
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']

        for region in self :
            if region.interpolation != xDataEnumsModule.Interpolation.linlin:
                print ( '    WARNING: ignoring interpolation "%s" for toPointwise_withLinearXYs' % region.interpolation )
        if( ( lowerEps < 0 ) or ( upperEps < 0 ) ) :
            raise ValueError( 'lowerEps = %s and upperEps = %s must be >= 0.' % ( lowerEps, upperEps ) )
        if( lowerEps == upperEps == 0 ) : raise ValueError( 'lowerEps and upperEps cannot both be 0.' )

        minEps = 5e-16
        if( lowerEps != 0 ) : lowerEps = max( minEps, lowerEps )
        if( upperEps != 0 ) : upperEps = max( minEps, upperEps )

        linear = self.toLinearXYsClass( )( axes = self.axes )
        n1 = len( self )
        for i1, region in enumerate( self ) :
            _region = region.toPointwise_withLinearXYs( **kwargs )
            if( i1 < n1 - 1 ) :
                if( lowerEps > 0  ) : _region[-1].outerDomainValue = _region[-1].outerDomainValue * ( 1 - lowerEps )
            if( i1 > 0 ) :
                if( upperEps > 0  ) : _region[0].outerDomainValue = _region[0].outerDomainValue * ( 1 + upperEps )
            for xys in _region : linear.append( xys )
        return( linear )

    def toLinearXYsClass( self ) :

        return( XYs2d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

"""
# BRB need to implement
    def getAtEnergy( self, energy ) :
    def invert( self ) :
"""

class Form( baseModule.Form ) :

    moniker = 'angular'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        if( not( isinstance( angularSubform, Subform ) ) ) : raise TypeError( 'instance is not an angular subform' )
        baseModule.Form.__init__( self, label, productFrame, ( angularSubform, ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angularSubform.convertUnits( unitMap )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise NotImplementedError('FIXME, when am I called?')

        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        product = kwargs['product']
        reactionSuite = kwargs['reactionSuite']

        if( product.pid != IDsPoPsModule.photon ) :
            raise ValueError( 'For form %s, calculateAverageProductData is only for gammas, not %s' % ( self.moniker, product.pid ) )

        depData = []
        angularSubform = self.angularSubform

        Es = angularSubform.getEnergyArray( kwargs['EMin'], kwargs['EMax'] )
        if( ( 'discrete' in product.attributes ) or ( 'primary' in product.attributes ) ) :
            massRatio = 0.
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
            else :
                Eg = product.attributes['primary'].getValueAs( energyUnit )
                mass1 = kwargs['projectileMass']
                mass2 = kwargs['targetMass']
                massRatio = mass2 / ( mass1 + mass2 )
            depEnergy = [ [ E, Eg + massRatio * E ] for E in Es ]
            depMomentum = [ [ E, ( Eg + massRatio * E ) * angularSubform.averageMu( E, accuracy = 0.1 * momentumAccuracy ) ] for E in Es ]
        else :
            raise ValueError( 'Unsupported gamma; gamma must be "discrete" or "primary"' )

        axes = averageProductEnergyModule.defaultAxes( energyUnit = energyUnit )
        depData.append( averageProductEnergyModule.XYs1d( data = depEnergy, axes = axes, label = style.label ) )
        axes = averageProductMomentumModule.defaultAxes( energyUnit = energyUnit, momentumUnit = momentumUnit )
        depData.append( averageProductMomentumModule.XYs1d( data = depMomentum, axes = axes, label = style.label ) )

        return( depData )

    def processMultiGroup( self, style, tempInfo, indent ) :

        raise NotImplementedError( 'need to implement' )

    @classmethod
    def parseNodeUsingClass(cls, angularElement, xPath, linkData, **kwargs):
        """Translate <angular> element from xml."""

        xPath.append( angularElement.tag )
        subformClass = {
                Isotropic2d.moniker   : Isotropic2d,
                Forward.moniker     : Forward,
                XYs2d.moniker       : XYs2d,
                Regions2d.moniker   : Regions2d,
            }.get( angularElement.tag )
        if( subformClass is None ) : raise ValueError( "encountered unknown angular subform: %s" % angularElement.tag )
        angularSubform = subformClass.parseNodeUsingClass(angularElement, xPath, linkData, **kwargs)
        xPath.pop( )
        return( angularSubform )

class TwoBody( baseModule.Form ) :

    moniker = 'angularTwoBody'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        angularSubform.label = None                 # FIXME: Need to check the type of angularSubform.
        baseModule.Form.__init__( self, label, productFrame, ( angularSubform, ) )

    @property
    def data(self):
        '''Temporary method until I redo all the classes that have a data-like member.'''

        return(self.angularSubform)

    @data.setter
    def data(self, data):
        '''Temporary method until I redo all the classes that have a data-like member. Need to check data type.'''

        self.angularSubform = data

    @property
    def domainMin(self):

        return self.angularSubform.domainMin

    def averageMu( self, E, accuracy = None ) :
 
        return( self.angularSubform.averageMu( E, accuracy = accuracy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        kwargs['productFrame'] = self.productFrame
        return( calculateAverageProductData( self.angularSubform, style, indent, **kwargs ) )

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        self.angularSubform.convertUnits( unitMap )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """Returns the energy spectrum in the lab frame for the specified incident energy."""

        from fudge import product as productModule
        from ... import outputChannel as outputChannelModule

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)
        partialMuDomain = muMin != -1.0 or muMax != 1.0

        reactionSuite = self.rootAncestor
        energyUnit = self.angularSubform.domainUnit
        massUnit = energyUnit + '/c**2'
        projectileMass = reactionSuite.PoPs[reactionSuite.projectile].getMass(massUnit)
        targetMass = reactionSuite.PoPs[reactionSuite.target].getMass(massUnit)
        if projectileMass == 0.0 or targetMass == 0.0:
            print('    WARNING: skipping unsupported relativistic TwoBody.energySpectrumAtEnergy: projectile mass = %s, target mass = %s' % 
                    (projectileMass, targetMass))
            return energyModule.XYs1d(axes=energyModule.defaultAxes(energyUnit))
        product = self.findClassInAncestry(productModule.Product)
        outputChannel = self.findClassInAncestry(outputChannelModule.OutputChannel)

        productMass = reactionSuite.PoPs[outputChannel.products[0].pid].getMass(massUnit)
        residualMass = reactionSuite.PoPs[outputChannel.products[1].pid].getMass(massUnit)
        if productMass == 0.0:                      # This is capture with product being primary photon.
            Q = outputChannel.Q[0].value
            inputMass = projectileMass + targetMass
            gammaEnergy = Q + energyIn * targetMass / inputMass * (1 - Q / inputMass)
            halfWidth = kwargs['twoBodyCOM_energyResolution']
            data = [[gammaEnergy - halfWidth, 0.0], [gammaEnergy, 1 / halfWidth], [gammaEnergy + halfWidth, 0.0]]
            return energyModule.XYs1d(data=data, axes=energyModule.defaultAxes(energyUnit))
        elif residualMass == 0.0:
            print('    WARNING: skipping unsupported relativistic TwoBody.energySpectrumAtEnergy: product mass = %s, target mass = %s' % 
                (productMass, residualMass))
            return energyModule.XYs1d(axes=energyModule.defaultAxes(energyUnit))

        if outputChannel.products[1] == product:
            productMass, residualMass = residualMass, productMass

        comKineticEnergy = energyIn * targetMass / (projectileMass + targetMass) + projectileMass + targetMass - productMass - residualMass
        comKineticEnergy *= residualMass / (productMass + residualMass)
        if comKineticEnergy < 0.0:
            return energyModule.XYs1d(axes=energyModule.defaultAxes(energyUnit))

        PofMu = self.angularSubform.evaluate(energyIn)
        if frame == xDataEnumsModule.Frame.centerOfMass:
            twoBodyCOMResolution = min(0.1, max(1e-6, kwargs.get('twoBodyCOMResolution', 1e-2)))
            data = [[comKineticEnergy * (1 - twoBodyCOMResolution), 0.0],
                    [comKineticEnergy, 1.0], 
                    [comKineticEnergy * (1 + twoBodyCOMResolution), 0.0]]
            return energyModule.XYs1d(data, axes=energyModule.defaultAxes(energyUnit)).normalize() * PofMu.integrate(domainMin=muMin, domainMax=muMax)

        projectileSpeed = math.sqrt(2.0 * energyIn / projectileMass)
        comSpeed = projectileMass / (projectileMass + targetMass) * projectileSpeed

        comSpeedProductEnergy = 0.5 * productMass * comSpeed * comSpeed
        productComEnergy = comKineticEnergy * residualMass / (productMass + residualMass)
        productCOM_speed = math.sqrt(2.0 * productComEnergy / productMass)
        sqrtEnergiesTimes2 = 2.0 * math.sqrt(comSpeedProductEnergy * productComEnergy)

        mus = self.productMuInCOM_fromProductMusInLab(muMin, muMax, comSpeed, productCOM_speed)

        PofMu = self.angularSubform.evaluate(energyIn)
        PofEnergyPrime = []
        numberOfPoints = 51
        numberOfPointsMinusOne = numberOfPoints - 1
        norm = 0.0
        for index, (mu1, mu2) in enumerate(mus):
            for muIndex in range(numberOfPoints):
                mu = (mu1 * (numberOfPointsMinusOne - muIndex) + mu2 * muIndex) / numberOfPointsMinusOne
                if abs(mu) < 1e-12:
                    mu = 0.0
                energyPrime = comSpeedProductEnergy + productComEnergy + mu * sqrtEnergiesTimes2
                PatMu = PofMu.evaluate(mu)
                PatE = PatMu / sqrtEnergiesTimes2
                PofEnergyPrime.append([energyPrime, PatE])
            if partialMuDomain:
                norm += PofMu.integrate(domainMin=mu1, domainMax=mu2)
                if len(mus) > 1:
                    if index == 0:
                        firstPofEnergyPrime = energyModule.XYs1d(data=PofEnergyPrime, axes=energyModule.defaultAxes(energyUnit))
                        PofEnergyPrime = []
                    else:
                        PofEnergyPrime = energyModule.XYs1d(data=PofEnergyPrime, axes=energyModule.defaultAxes(energyUnit))
                        firstPofEnergyPrime, PofEnergyPrime = firstPofEnergyPrime.mutualify(0.0, 1e-6, True, PofEnergyPrime, 1e-6, 0, True)
                        PofEnergyPrime = firstPofEnergyPrime + PofEnergyPrime
        PofEnergyPrime = energyModule.XYs1d(data=PofEnergyPrime, axes=energyModule.defaultAxes(energyUnit))
        if len(PofEnergyPrime) > 0:
            PofEnergyPrime.normalize(insitu=True)
            if partialMuDomain:
                PofEnergyPrime *= norm

        return PofEnergyPrime

    def evaluate( self, energy, mu, phi = 0, frame = None ) :

        if( frame is None ) : frame = self.productFrame
        if( frame == self.productFrame ) :
            if( isinstance( self.angularSubform, Recoil ) ) :
                if( isinstance( mu, float ) ) : 
                    mu = -mu
                else :
                    raise TypeError( 'Currently, mu can only be a float. Got "%s"' % mu )
                return( self.angularSubform.partner.evaluate( energy, mu, phi, frame ) )
            else :
                return( self.angularSubform.evaluate( energy ).evaluate( mu ) )
        else :
            raise ValueError( 'Only frame of product is currently supported' )

    def findLinks( self, links ) :

        item = self.angularSubform
        if( isinstance( item, ( linkModule.Link, linkModule.Link2 ) ) ) : links.append( [ item, item.link, item.path ] )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        This method call \*\*fixDomains\* on the \*\*angularSubform\* member.
        """

        return self.angularSubform.fixDomains(energyMin, energyMax, domainToFix)

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.angularSubform.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :

        from fudge import product as productModule
        from ... import outputChannel as outputChannelModule

        def product_mu_in_com_from_final_product_energy_lab( product_energy_prime_lab, product_mass, com_speed, product_com_speed ) :

            product_speed_prime_lab = math.sqrt( 2.0 * product_energy_prime_lab / product_mass )
            mu = 0.5 * ( product_speed_prime_lab * product_speed_prime_lab - com_speed * com_speed - product_com_speed * product_com_speed ) / ( com_speed * product_com_speed )
            return( min( 1.0, max( -1.0, mu ) ) )

        energy_unit = self.angularSubform.domainUnit
        mass_unit = energy_unit + '/c**2'
        projectile_mass = reaction_suite.PoPs[reaction_suite.projectile].getMass( mass_unit )
        target_mass = reaction_suite.PoPs[reaction_suite.target].getMass( mass_unit )
        product = self.findClassInAncestry(productModule.Product)
        output_channel = self.findClassInAncestry(outputChannelModule.OutputChannel)

        product_mass = reaction_suite.PoPs[output_channel.products[0].pid].getMass( mass_unit )
        residual_mass = reaction_suite.PoPs[output_channel.products[1].pid].getMass( mass_unit )
        if( output_channel.products[1] == product ) : product_mass, residual_mass = residual_mass, product_mass

        com_kinetic_energy = energyIn * target_mass / ( projectile_mass + target_mass ) + projectile_mass + target_mass - product_mass - residual_mass
        if( com_kinetic_energy < 0.0 ) : return( 0.0 )

        projectile_speed = math.sqrt( 2.0 * energyIn / projectile_mass )
        com_speed = projectile_mass / ( projectile_mass + target_mass ) * projectile_speed

        product_com_energy = com_kinetic_energy * residual_mass / ( product_mass + residual_mass )
        product_com_speed = math.sqrt( 2.0 * product_com_energy / product_mass )

        product_energy_prime_min = 0.5 * product_mass * math.pow( com_speed - product_com_speed, 2.0 )      # Minimum product energy in the lab frame occurs for mu_com = -1.
        product_energy_prime_max = 0.5 * product_mass * math.pow( com_speed + product_com_speed, 2.0 )      # Maximum product energy in the lab frame occurs for mu_com = 1.
        energy_out_min, energy_out_max = miscellaneousModule.domainLimits( energyOut, product_energy_prime_min, product_energy_prime_max )
        if( product_energy_prime_max <= energy_out_min ) : return( 0.0 )
        mu_energy_prime_1 = product_mu_in_com_from_final_product_energy_lab( energy_out_min, product_mass, com_speed, product_com_speed )

        mu_energy_prime_2 = None
        if( energy_out_max is not None ) :
            if( product_energy_prime_min >= energy_out_max ) : return( 0.0 )
            mu_energy_prime_2 = product_mu_in_com_from_final_product_energy_lab( energy_out_max, product_mass, com_speed, product_com_speed )
        else :
            mu_energy_prime_2 = min(  1.0, mu_energy_prime_1 + 1e-12 * math.fabs( mu_energy_prime_1 ) )
            mu_energy_prime_1 = max( -1.0, mu_energy_prime_1 - 1e-12 * math.fabs( mu_energy_prime_1 ) )

        angular1d = self.angularSubform.evaluate( energyIn )

        if( muOut is None ) :
            raise ValueError( 'Legendre treatment not currently supported.' )
        else :
            mu_prime_min, mu_prime_max = miscellaneousModule.domainLimits( muOut, angular1d.domainMin, angular1d.domainMax )
            if( product_com_speed <= com_speed ) :
                mu_lab_max = math.sqrt( com_speed * com_speed - product_com_speed * product_com_speed ) / com_speed     # Maximum mu in lab frame that product can have.
                if( mu_prime_min >= mu_lab_max ) : return( 0.0 )
                mu_prime_1, mu_prime_4 = self.productMuInCOM_fromProductMuInLab( mu_prime_min, com_speed, product_com_speed )
                if( mu_prime_max < mu_lab_max ) :
                    mu_prime_2, mu_prime_3 = self.productMuInCOM_fromProductMuInLab( mu_prime_max, com_speed, product_com_speed )
                    muLists = [ [ mu_prime_1, mu_prime_2 ], [ mu_prime_3, mu_prime_4 ] ]
                else :
                    muLists = [ [ mu_prime_1, mu_prime_4 ] ]
            else :
                mu_prime_1, dummy = self.productMuInCOM_fromProductMuInLab( mu_prime_min, com_speed, product_com_speed )
                mu_prime_2 = None
                if mu_prime_max is not None:
                    mu_prime_2, dummy = self.productMuInCOM_fromProductMuInLab( mu_prime_max, com_speed, product_com_speed )
                muLists = [ [ mu_prime_1, mu_prime_2 ] ]

        angular_evaluate = 0.0

        for mu_prime_1, mu_prime_2 in muLists :
            if( mu_prime_2 is None ) :
                if( mu_energy_prime_1 <= mu_prime_1 <= mu_energy_prime_2 ) : angular_evaluate += angular1d.evaluate( mu_prime_1 )
            else :
                mu_min = max( mu_prime_1, mu_energy_prime_1 )
                mu_max = min( mu_prime_2, mu_energy_prime_2 )
                if( mu_max > mu_min ) : angular_evaluate += angular1d.integrate( mu_min, mu_max )

        phi_evaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( angular_evaluate * phi_evaluate )

    def invert( self ) :

        return( self.angularSubform.invert( ) )

    def isIsotropic( self ) :

        return( self.angularSubform.isIsotropic( ) )

    def isTwoBody( self ) :

        return( True )

    def processMC_cdf( self, style, tempInfo, indent ) :

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        if( verbosity > 2 ) : print('%sGrouping %s' % ( indent, self.moniker ))

        subform = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( subform is None ) : return( None )

        form = TwoBody( style.label, self.productFrame, subform )
        return( form )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']
        outputChannel = tempInfo['reaction'].outputChannel

        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

        angularSubform = self.angularSubform
        if( isinstance( angularSubform, Recoil ) ) :
            angularSubform = angularSubform.getNumericalDistribution( )

        Q = outputChannel.Q.evaluated
        Q = Q.evaluate( Q.domainMin )

        residual = outputChannel.products[1]
        if( tempInfo['productIndex'] == '1' ) : residual = outputChannel.products[0]
        targetMass = tempInfo['masses']['Target']
        residualMass = tempInfo['masses']['Residual']
        if( tempInfo['isInfiniteTargetMass'] ) :                      # Happends for electro-atomic reactions.
            tempInfo['masses']['Target'] = 1e9 * tempInfo['masses']['Projectile']
            tempInfo['masses']['Residual'] = tempInfo['masses']['Target']
        else :
            tempInfo['masses']['Residual'] = residual.getMass( tempInfo['massUnit'] )
        TM_1, TM_E = transferMatricesModule.twoBodyTransferMatrix( style, tempInfo, self.productFrame, tempInfo['crossSection'], angularSubform, 
                Q, comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )

        tempInfo['masses']['Target'] = targetMass
        tempInfo['masses']['Residual'] = residualMass
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.angularSubform.toPointwise_withLinearXYs( **kwargs ) )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Translate <angularTwoBody> element from xml. Returns instance of angular.TwoBody.
        """

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {    Isotropic2d.moniker   : Isotropic2d,
                            Recoil.moniker      : Recoil,
                            XYs2d.moniker       : XYs2d,
                            Regions2d.moniker   : Regions2d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise ValueError( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        angularForm = cls(node.get('label'), node.get('productFrame'), angularSubform)

        xPath.pop( )
        return( angularForm )

    @staticmethod
    def productMuInCOM_fromProductMuInLab(muLab, COM_speed, productCOM_speed):
        '''
        This function returns the tuple (mu1, mu2, isTangent) where mu1 and mu2 are the center-of-
        mass (COM) frame mus for an outgoing particle (product) that correspond to the lab frame mu *muLab*.
        The interaction is two-body, *COM_speed* is the COM mass speed and *productCOM_speed*
        is the product's speed in the COM frame.  If *COM_speed* <= productCOM_speed, 
        there is only one COM mu (mu1) for *muLab*. In this case, the returned tuple is (mu1, None, False).

        Otherwise, there are 0, 1 or 2 solutions as follows. If *muLab* is too large, there is no solution
        and the tuple (None, None, False) is returned. If *muLab* is tangent to the circle of the product's
        velolicty in the COM frame then (mu1, None, True) are returned. Otherwise (mu1, mu2, False) is
        returned with mu1 < mu2.
        '''

        speedRatio = COM_speed / productCOM_speed
        factor = speedRatio * (1.0 - muLab * muLab)
        sqrtArgument = 1.0 - speedRatio * factor

        if speedRatio >= 1 and muLab < 0.0:
            return None, None, False

        if sqrtArgument <= 0.0:
            if sqrtArgument < 0.0:                      # No solution.
                return None, None, False
            else:                                       # Tangent to productCOM_speed circle in COM frame.
                return math.sqrt( 1 - 1 / (speedRatio * speedRatio)), None, True
        term_pm = muLab * math.sqrt(sqrtArgument)
        mu2 = term_pm - factor                          # Solution for productCOM_speed > COM_speed or one solution for productCOM_speed < COM_speed.

        if productCOM_speed >= COM_speed:
            return mu2, None, False

        mu1 = -term_pm - factor                         # Other solution for productCOM_speed < COM_speed (mu1 < mu2).

        return mu1, mu2, False

    @staticmethod
    def productMuInCOM_fromProductMusInLab(muMinLab, muMaxLab, COM_speed, productCOM_speed):
        '''
        This function determines all center-of-mass (COM) frame mu domains (i.e., [mu1, mu2] pairs) that lay 
        between the lab frame mus *muMinLab* and *muMaxLab*. The returned objects is a list that can contain 0, 
        1 or 2 pairs of mu values. Note if *muMinLab* >= *muMaxLab* an empty list is returned. If 
        *COM_speed* < *productCOM_speed* the list [[mu1, mu2]] is returned. Otherwise, there are 3 possible solutions:
        1) no COM mus corresonding to *muMaxLab* (and hence *muMinLab*) in which case an empty list is also returned, 
        2) only *muMaxLab* has corresonding COM mus and [[mu1, mu2]] is returned where mu1 < mu2 or 3) both *muMinLab*
        and *muMaxLab* have corresonding COM mus and [[mu1, mu2], [mu3, mu4]] where mu1 < mu2 < mu3 < mu4 is returned.
        '''

        if muMinLab >= muMaxLab:
            return []

        mu1, mu2, isTangentMin = TwoBody.productMuInCOM_fromProductMuInLab(muMinLab, COM_speed, productCOM_speed)
        mu3, mu4, isTangentMax = TwoBody.productMuInCOM_fromProductMuInLab(muMaxLab, COM_speed, productCOM_speed)

        if mu3 is None or isTangentMax:
            return []

        if mu1 is None or isTangentMin:
            return [[mu3, mu4]]

        if mu2 is None:
            return [[mu1, mu3]]

        return [[mu3, mu1], [mu2, mu4]]

def calculateAverageProductData( self, style, indent, **kwargs ) :
    """
    This function calculates average product data for two-body reactions.
    """

    def calculateDepositionEnergyAtE( angularData, E, parameters ) :

        a1x, a2y, Qp = parameters['a1x'], parameters['a2y'], parameters['Qp']
        dE = max( 0., E - Qp )
        return( a1x * E + a2y * dE + 2. * math.sqrt( a1x * a2y * E * dE ) * angularData.averageMu( E, accuracy = energyAccuracy ) )

    def calculateDepositionEnergyAtEForPhotoProjectile( angularData, E, parameters ) :

        mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
        E__E_m2 = E / ( E + mass2 )
        dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
        return( 0.5 * E__E_m2**2 * massx + dE + E__E_m2 * math.sqrt( 2. * dE * massx ) * angularData.averageMu( E, accuracy = energyAccuracy ) )

    def calculateDepositionEnergyAtE_forInfiniteTargetMass( angularData, E, parameters ) :

        return( E )

    class CalculateDepositionEnergyThicken :

        def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

            self.data = data
            self.angular = angular
            self.func = func
            self.parameters = parameters
            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

        def evaluateAtX( self, E ) :

            return( self.func( self.angular, E, self.parameters ) )

    def calculateDepositionMomentumAtE( angularData, E, parameters ) :

        mass1, massx, b1x, a2y, Qp = parameters['m1'], parameters['mx'], parameters['b1x'], parameters['a2y'], parameters['Qp']
        dE = max( 0., E - Qp )
        return( b1x * math.sqrt( 2. * E / mass1 ) + math.sqrt( 2. * massx * a2y * dE ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    def calculateDepositionMomentumAtEForPhotoProjectile( angularData, E, parameters ) :

        mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
        E__E_m2 = E / ( E + mass2 )
        dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
        return( E__E_m2 * massx + math.sqrt( 2. * dE * massx ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    def calculateDepositionMomentumAtE_forInfiniteTargetMass( angularData, E, parameters ) :

        return( math.sqrt( 2.0 * mass1 * E ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    class CalculateDepositionMomentumThicken :

        def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

            self.data = data
            self.angular = angular
            self.func = func
            self.parameters = parameters
            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

        def evaluateAtX( self, E ) :

            return( self.func( self.angular, E, self.parameters ) )

    if( hasattr( self, 'calculateAverageProductData' ) ) :      # This happends when, for example, the angular is a Regions2d form.
        return( self.calculateAverageProductData( style, indent = indent, **kwargs ) )

    energyUnit = kwargs['incidentEnergyUnit']
    massUnit = kwargs['massUnit']
    energyAccuracy = kwargs['energyAccuracy']
    momentumAccuracy = kwargs['momentumAccuracy']
    reactionSuite = kwargs['reactionSuite']
    reaction = kwargs['reaction']
    outputChannel = kwargs['outputChannel']
    productIndex = kwargs['productIndex']

    mass1 = kwargs['projectileMass']
    mass2 = kwargs['targetMass']

    productx = reactionSuite.PoPs[outputChannel.products[0].pid]
    if productx.id in reactionSuite.PoPs.aliases: productx = reactionSuite.PoPs[productx.id]
    massx = productx.getMass( massUnit )

    producty = reactionSuite.PoPs[outputChannel.products[1].pid]
    if( not kwargs['isInfiniteTargetMass'] ) :
        if producty.id in reactionSuite.PoPs.aliases: producty = reactionSuite.PoPs[producty.id]
        massy = producty.getMass( massUnit )

        if( productIndex == '1' ) : massx, massy = massy, massx
        m12 = mass1 + mass2
        mxy = massx + massy
        b1x = mass1 * massx / m12
        a1x = mass1 * massx / ( m12 * m12 )
        a2y = mass2 * massy / ( m12 * mxy )
        Qm = m12 - mxy
        Q = reaction.getQ( energyUnit, final = False )
        Qp = -float( Q ) * m12 / mass2              # This is the threshold in the COM frame.

    if( kwargs['isInfiniteTargetMass'] ) :
        energyFunc = calculateDepositionEnergyAtE_forInfiniteTargetMass
        momentumFunc = calculateDepositionMomentumAtE_forInfiniteTargetMass
        parameters = {}
    elif( mass1 == 0. ) :                         # Photo as projectile
        energyFunc = calculateDepositionEnergyAtEForPhotoProjectile
        momentumFunc = calculateDepositionMomentumAtEForPhotoProjectile
        parameters = { 'm2' : mass2, 'mx' : massx, 'my' : massy, 'Q' : Q }
    else :
        energyFunc = calculateDepositionEnergyAtE
        momentumFunc = calculateDepositionMomentumAtE
        parameters = { 'm1' : mass1, 'mx' : massx, 'a1x' : a1x, 'a2y' : a2y, 'b1x' : b1x, 'Qp' : Qp }

    Es = self.getEnergyArray( kwargs['EMin'], kwargs['EMax'] )
    aveEnergy = [ [ E, energyFunc( self, E, parameters ) ] for E in Es ]
    aveMomentum = [ [ E, momentumFunc( self, E, parameters ) ] for E in Es ]

    if( not kwargs['isInfiniteTargetMass'] ) :
        aveEnergy = fudgemath.thickenXYList( aveEnergy, CalculateDepositionEnergyThicken( aveEnergy, self, energyFunc, 
                parameters, energyAccuracy, 1e-10 ) )
        aveMomentum = fudgemath.thickenXYList( aveMomentum, CalculateDepositionMomentumThicken( aveMomentum, self, momentumFunc, 
                parameters, momentumAccuracy, 1e-10 ) )

    return( [ aveEnergy ], [ aveMomentum ] )
