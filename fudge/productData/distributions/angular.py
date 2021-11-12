# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
Angular distribution classes with forms :py:class:`form`, :py:class:`twoBodyForm` and :py:class:`CoulombElasticForm`.

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

    * :py:class:`regions2d`

        - Same as ``piecewise``


Note, there was no ``pdfOfMu.piecewise`` and its new equivalent :py:class:`regions1d`, as
they are not currently needed. There are other classes in :py:mod:`angular.py`
(e.g., :py:class:`isotropic2d`, :py:class:`forward`, :py:class:`recoil`) but we can ignore them for now.

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
    :py:class:`regions2d` that stores the LTT=1 in the first region and LTT=2 in the second region.

"""

import math
from fudge.core.math import fudgemath

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from numericalFunctions import Legendre as LegendreModule

from xData import standards as standardsModule
from xData import ancestry as ancestryModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import series1d as series1dModule
from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule
from xData import link as linkModule

from .. import energyDeposition as energyDepositionModule
from fudge.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule
from . import reference as referenceModule
from . import miscellaneous as miscellaneousModule
from . import probabilities as probabilitiesModule
from . import energy as energyModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 3 )
    axes[0] = axesModule.axis( 'P(mu|energy_in)', 0, '' )
    axes[1] = axesModule.axis( 'mu', 1, '' )
    axes[2] = axesModule.axis( 'energy_in', 2, energyUnit )
    return( axes )

class XYs1d( XYsModule.XYs1d ) :

    def averageMu( self ) :

        allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                  standardsModule.interpolation.flatToken ]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYsModule.defaultAccuracy )
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

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( xs_pdf_cdf1d.fromXYs( self, outerDomainValue = self.outerDomainValue ) )


class Legendre( series1dModule.LegendreSeries ) :

    def averageMu( self ) :

        return( self.getCoefficientSafely( 1 ) )

    def LegendreCoefficient( self, LegendreOrder ) :

        if( LegendreOrder < len( self ) ) : return( self[LegendreOrder] )
        return( 0.0 )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, lowerEps = 0, upperEps = 1e-8 )
        linear.outerDomainValue = self.outerDomainValue
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

class xs_pdf_cdf1d( xs_pdf_cdfModule.xs_pdf_cdf1d ) :

    def isIsotropic( self ) :

        return( min( self.pdf.values ) == max( self.pdf.values ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class isotropic1d( ancestryModule.ancestry ) :

    moniker = 'isotropic1d'
    ancestryMembers = ( '', )

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

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

    def integrate( self, muMin, muMax ) :

        return( 0.5 * ( muMax - muMin ) )

    def LegendreCoefficient( self, LegendreOrder ) :

        if( LegendreOrder == 0 ) : return( 1.0 )
        return( 0.0 )

class subform( baseModule.subform ) :
    """Abstract base class for angular subforms."""

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( None )

class isotropic2d( subform ) :

    moniker = 'isotropic2d'
    ancestryMembers = ( '', )

    def __init__( self ) :

        subform.__init__( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( isotropic2d( ) )

    __copy__ = copy

    @property
    def domainGrid( self ) :

        return( [ self.domainMin, self.domainMax ] )

    @property
    def domainMin( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainMin )

    @property
    def domainMax( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainMax )

    @property
    def domainUnit( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainUnit )

    def evaluate( self, energy ) :

        _isotropic1d = isotropic1d( )
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

    def toXMLList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        return( isotropic2d( ) )

class forward( subform ) :

    moniker = 'forward'
    ancestryMembers = ( '', )

    def __init__( self ) :

        subform.__init__( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( forward( ) )

    __copy__ = copy

    @property
    def domainMin( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainMin )

    @property
    def domainMax( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainMax )

    @property
    def domainUnit( self ) :

        from fudge import product
        return( self.findClassInAncestry( product.product ).domainUnit )

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

        ptw = XYs2d( defaultAxes( self.domainUnit ) )
        ptw.append( XYsModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMin ) )
        ptw.append( XYsModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMax ) )
        return( ptw )

    def toXMLList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        return forward()

class recoil( linkModule.link, subform ) :

    moniker = 'recoil'
    ancestryMembers = ( '', )
    def __init__( self, link, root = None, path = None, label = None, relative = False ) :
        """
        Distributions for this particle can be computed from its recoil partner using kinematics.
        Only meant for angular 2-body reactions.
        The 'root' and 'path' usually do not need to be specified: they will be computed from the link.

        :param link: distribution form for recoil partner
        :type link: distribution.base.form
        :param root: filename where recoil partner is defined (usually None)
        :param path: xpath to recoil partner (usually will be computed from the link)
        :param label: label for this subform (usually None)
        :param relative: whether to write xpath as relative or absolute
        :return:
        """

        linkModule.link.__init__( self, link=link, root=root, path=path, label=label, relative=relative )
        subform.__init__( self )

    @property
    def domainUnit( self ) :

        return( self.partner.angularSubform.domainUnit )

    @property
    def partner( self ) :

        if( self.link is None ) : raise Exception( "Encountered unresolved link!" )
        return( self.link )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def copy( self ):

        return( recoil( self.partner ) )

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
        if( not( isinstance( self.partner, ( twoBodyForm, referenceModule.CoulombPlusNuclearElastic ) ) ) ) :
            warnings.append( warning.missingRecoilDistribution( self ) )

        return warnings

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.partner.invert( ).toPointwise_withLinearXYs( **kwargs ) )

class XYs2d( subform, probabilitiesModule.PofX1GivenX2 ) :

    def __init__( self, **kwargs ) :

        probabilitiesModule.PofX1GivenX2.__init__( self, **kwargs )
        subform.__init__( self )

    def getAtEnergy( self, energy ) :

        return( self.interpolateAtValue( energy, extrapolation = standardsModule.flatExtrapolationToken ) )

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
            if( interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
                raise TypeError( 'Unsupported interpolation qualifier "%s"' % interpolationQualifier )
            elif( interpolationQualifier != standardsModule.interpolation.noneQualifierToken ) :
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
                warnings.append( warning.unnormalizedDistribution( PQUModule.PQU( function.outerDomainValue, self.axes[-1].unit ), idx, integral, function ) )

            if( xys.rangeMin < 0.0 ) :
                warnings.append( warning.negativeProbability( PQUModule.PQU( function.outerDomainValue, self.axes[-1].unit), value = xys.rangeMin, obj=function ) )

        return warnings

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        subform = XYs2d( axes = self.axes.copy( ), interpolation = self.interpolation )
        for xys in self : subform.append( xys.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( subform )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, Legendre, xs_pdf_cdf1d ) )

class regions2d( subform, regionsModule.regions2d ) :

    def __init__( self, **kwargs ) :

        regionsModule.regions2d.__init__( self, **kwargs )
        subform.__init__( self )

    def check( self, info ):

        from fudge import warning

        warnings = []
        for idx,region in enumerate(self):
            regionWarnings = region.check( info )
            if regionWarnings: warnings.append( warning.context('region index %i: %s'
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

        _regions2d = regions2d( axes = self.axes )
        for region in self : _regions2d.append( region.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( _regions2d )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        arguments = self.getArguments( kwargs, { 'lowerEps' : 0, 'upperEps' : 0 } )
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']

        for region in self :
            if( region.interpolation != standardsModule.interpolation.linlinToken ) :
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

class form( baseModule.form ) :

    moniker = 'angular'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        if( not( isinstance( angularSubform, subform ) ) ) : raise TypeError( 'instance is not an angular subform' )
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, ) )

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

        if( product.id != IDsPoPsModule.photon ) :
            raise ValueError( 'For form %s, calculateAverageProductData is only for gammas, not %s' % ( self.moniker, product.id ) )

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

        axes = energyDepositionModule.defaultAxes( energyUnit = energyUnit )
        depData.append( energyDepositionModule.XYs1d( data = depEnergy, axes = axes, label = style.label ) )
        axes = momentumDepositionModule.defaultAxes( energyUnit = energyUnit, momentumUnit = momentumUnit )
        depData.append( momentumDepositionModule.XYs1d( data = depMomentum, axes = axes, label = style.label ) )

        return( depData )

    def processMultiGroup( self, style, tempInfo, indent ) :

        raise NotImplementedError( 'need to implement' )

#    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge import miscellaneous as miscellaneousModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        product = tempInfo['product']
        if( product.id != IDsPoPsModule.photon ) : raise ValueError( 'For form %s, process is only for gammas, not %s' % ( self.moniker, product.id ) )

        newForms = []
        angularSubform = self.angularSubform
        crossSection = tempInfo['crossSection']
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print ('%sGrouping %s' % ( verbosityIndent, self.moniker ))
            projectileName, productName = processInfo.getProjectileName( ), product.id
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
                # FIXME passing too many arguments!
                TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( processInfo, projectileName, productName, Eg, crossSection,
                   angularSubform, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            elif( 'primary' in product.attributes ) :
                projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
                mass1, mass2 = projectile.getMass( massUnit ), target.getMass( massUnit )
                massRatio = mass2 / ( mass1 + mass2 )
                ELevel = product.attributes['primary'].getValueAs( energyUnit )
                # FIXME passing too many arguments!
                TM_1, TM_E = transferMatricesModule.primaryGammaAngularData( processInfo, projectileName, productName, ELevel, massRatio, crossSection,
                   angularSubform, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            else :
                raise ValueError( 'Unsupported gamma; gamma must be "discrete" or "primary"' )
            miscellaneousModule.TMs2Form( processInfo, tempInfo, newForms, TM_1, TM_E, crossSection.axes )

        return( newForms )

    @staticmethod
    def parseXMLNode( angularElement, xPath, linkData ) :
        """Translate <angular> element from xml."""

        xPath.append( angularElement.tag )
        subformClass = {
                isotropic2d.moniker   : isotropic2d,
                forward.moniker     : forward,
                XYs2d.moniker       : XYs2d,
                regions2d.moniker   : regions2d,
            }.get( angularElement.tag )
        if( subformClass is None ) : raise ValueError( "encountered unknown angular subform: %s" % angularElement.tag )
        angularSubform = subformClass.parseXMLNode( angularElement, xPath, linkData )
        xPath.pop( )
        return( angularSubform )

class twoBodyForm( baseModule.form ) :

    moniker = 'angularTwoBody'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        angularSubform.label = None
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, ) )

    def averageMu( self, E, accuracy = None ) :
 
        return( self.angularSubform.averageMu( E, accuracy = accuracy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        kwargs['productFrame'] = self.productFrame
        return( calculateAverageProductData( self.angularSubform, style, indent, **kwargs ) )

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        self.angularSubform.convertUnits( unitMap )

    def energySpectrumAtEnergy( self, energyIn, frame, **kwargs ) :
        """Returns the energy spectrum in the lab frame for the specified incident energy."""

        from fudge import product as productModule
        from ... import outputChannel as outputChannelModule

        reactionSuite = self.getRootAncestor( )
        energyUnit = self.angularSubform.domainUnit
        massUnit = energyUnit + '/c**2'
        projectileMass = reactionSuite.PoPs[reactionSuite.projectile].getMass( massUnit )
        targetMass = reactionSuite.PoPs[reactionSuite.target].getMass( massUnit )
        if( ( projectileMass == 0.0 ) or ( targetMass == 0.0 ) ) :
            print( '    WARNING: skipping unsupported relativistic twoBodyForm.energySpectrumAtEnergy: projectile mass = %s, target mass = %s' % ( projectileMass, targetMass ) )
            return( energyModule.XYs1d( axes = energyModule.defaultAxes( energyUnit ) ) )
        product = self.findClassInAncestry( productModule.product )
        outputChannel = self.findClassInAncestry( outputChannelModule.outputChannel )

        productMass = reactionSuite.PoPs[outputChannel.products[0].pid].getMass( massUnit )
        residualMass = reactionSuite.PoPs[outputChannel.products[1].pid].getMass( massUnit )
        if( ( productMass == 0.0 ) or ( residualMass == 0.0 ) ) :
            print( '    WARNING: skipping unsupported relativistic twoBodyForm.energySpectrumAtEnergy: product mass = %s, target mass = %s' % ( productMass, residualMass ) )
            return( energyModule.XYs1d( axes = energyModule.defaultAxes( energyUnit ) ) )
        if( outputChannel.products[1] == product ) : productMass, residualMass = residualMass, productMass

        comKineticEnergy = energyIn * targetMass / ( projectileMass + targetMass ) + projectileMass + targetMass - productMass - residualMass
        comKineticEnergy *= residualMass / ( productMass + residualMass )
        if( comKineticEnergy < 0.0 ) : return( energyModule.XYs1d( axes = energyModule.defaultAxes( energyUnit ) ) )

        if( frame == standardsModule.frames.centerOfMassToken ) :
            twoBodyCOMResolution = min( 0.1, max( 1e-6, kwargs.get( 'twoBodyCOMResolution', 1e-2 ) ) )
            data = [ [ comKineticEnergy * ( 1 - twoBodyCOMResolution ), 0.0 ],
                     [ comKineticEnergy, 1.0 ], 
                     [ comKineticEnergy * ( 1 + twoBodyCOMResolution ), 0.0 ] ]
            return( energyModule.XYs1d( data, axes = energyModule.defaultAxes( energyUnit ) ).normalize( ) )

        projectileSpeed = math.sqrt( 2.0 * energyIn / projectileMass )
        comSpeed = projectileMass / ( projectileMass + targetMass ) * projectileSpeed

        comSpeedProductEnergy = 0.5 * productMass * comSpeed * comSpeed
        productComEnergy = comKineticEnergy * residualMass / ( productMass + residualMass )
        sqrtEnergiesTimes2 = 2.0 * math.sqrt( comSpeedProductEnergy * productComEnergy )

        PofMu = self.angularSubform.evaluate( energyIn )
        PofEnergyPrime = []
        numberOfPoints = 51
        for muIndex in range( numberOfPoints ) :
            mu = -1.0 + 2.0 * muIndex / ( numberOfPoints - 1 )
            if( mu > 0.99999 ) : mu = 1.0
            if( abs( mu ) < 1e-8 ) : mu = 0.0
            energyPrime = comSpeedProductEnergy + productComEnergy + mu * sqrtEnergiesTimes2
            PatMu = PofMu.evaluate( mu )
            PatE = PatMu / sqrtEnergiesTimes2
            PofEnergyPrime.append( [ energyPrime, PatE ] )
        PofEnergyPrime = energyModule.XYs1d( data = PofEnergyPrime, axes = energyModule.defaultAxes( energyUnit ) )
        PofEnergyPrime.normalize( insitu = True )

        return( PofEnergyPrime )

    def evaluate( self, energy, mu, phi = 0, frame = None ) :

        if( frame is None ) : frame = self.productFrame
        if( frame == self.productFrame ) :
            if( isinstance( self.angularSubform, recoil ) ) :
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
        if( isinstance( item, ( linkModule.link, linkModule.link2 ) ) ) : links.append( [ item, item.link, item.path ] )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.angularSubform.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = standardsModule.frames.productToken, LegendreOrder = 0 ) :

        from fudge import product as productModule
        from ... import outputChannel as outputChannelModule

        def product_mu_in_com_from_product_mu_in_lab( mu_prime_lab, com_speed, product_com_speed ) :

            speedRatio = com_speed / product_com_speed
            factor = speedRatio * ( 1.0 - mu_prime_lab * mu_prime_lab )
            term_pm = mu_prime_lab * math.sqrt( 1.0 - speedRatio * factor )
            mu1 = term_pm - factor                      # Solution for mu_prime_lab > com_speed and one solution for mu_prime_lab < com_speed.
            if( product_com_speed >= com_speed ) : return( mu1 )
            mu2 = -term_pm - factor                     # Other solution for mu_prime_lab < com_speed (mu2 < mu1).

            return( mu2, mu1 )

        def product_mu_in_com_from_final_product_energy_lab( product_energy_prime_lab, product_mass, com_speed, product_com_speed ) :

            product_speed_prime_lab = math.sqrt( 2.0 * product_energy_prime_lab / product_mass )
            mu = 0.5 * ( product_speed_prime_lab * product_speed_prime_lab - com_speed * com_speed - product_com_speed * product_com_speed ) / ( com_speed * product_com_speed )
            return( min( 1.0, max( -1.0, mu ) ) )

        energy_unit = self.angularSubform.domainUnit
        mass_unit = energy_unit + '/c**2'
        projectile_mass = reaction_suite.PoPs[reaction_suite.projectile].getMass( mass_unit )
        target_mass = reaction_suite.PoPs[reaction_suite.target].getMass( mass_unit )
        product = self.findClassInAncestry( productModule.product )
        output_channel = self.findClassInAncestry( outputChannelModule.outputChannel )

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
                mu_prime_1, mu_prime_4 = product_mu_in_com_from_product_mu_in_lab( mu_prime_min, com_speed, product_com_speed )
                if( mu_prime_max < mu_lab_max ) :
                    mu_prime_2, mu_prime_3 = product_mu_in_com_from_product_mu_in_lab( mu_prime_max, com_speed, product_com_speed )
                    muLists = [ [ mu_prime_1, mu_prime_2 ], [ mu_prime_3, mu_prime_4 ] ]
                else :
                    muLists = [ [ mu_prime_1, mu_prime_4 ] ]
            else :
                mu_prime_1 = product_mu_in_com_from_product_mu_in_lab( mu_prime_min, com_speed, product_com_speed )
                mu_prime_2 = None
                if( mu_prime_max is not None ) : mu_prime_2 = product_mu_in_com_from_product_mu_in_lab( mu_prime_max, com_speed, product_com_speed )
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

        form = twoBodyForm( style.label, self.productFrame, subform )
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
        if( isinstance( angularSubform, recoil ) ) :
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

    @staticmethod
    def parseXMLNode( angularElement, xPath, linkData ) :
        """
        Translate <angularTwoBody> element from xml. Returns instance of angular.twoBodyForm.
        """

        xPath.append( angularElement.tag )

        subformElement = angularElement[0]
        subformClass = {    isotropic2d.moniker   : isotropic2d,
                            recoil.moniker      : recoil,
                            XYs2d.moniker       : XYs2d,
                            regions2d.moniker   : regions2d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise ValueError( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )
        angularForm = twoBodyForm( angularElement.get( 'label' ), 
                angularElement.get( 'productFrame' ), angularSubform )
        xPath.pop( )
        return( angularForm )

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

    class calculateDepositionEnergyThicken :

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

    class calculateDepositionMomentumThicken :

        def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

            self.data = data
            self.angular = angular
            self.func = func
            self.parameters = parameters
            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

        def evaluateAtX( self, E ) :

            return( self.func( self.angular, E, self.parameters ) )

    if( hasattr( self, 'calculateAverageProductData' ) ) :      # This happends when, for example, the angular is a regions2d form.
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

    productx = reactionSuite.PoPs[outputChannel.products[0].id]
    if productx.id in reactionSuite.PoPs.aliases: productx = reactionSuite.PoPs[productx.pid]
    massx = productx.getMass( massUnit )

    producty = reactionSuite.PoPs[outputChannel.products[1].id]
    if( not kwargs['isInfiniteTargetMass'] ) :
        if producty.id in reactionSuite.PoPs.aliases: producty = reactionSuite.PoPs[producty.pid]
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
        aveEnergy = fudgemath.thickenXYList( aveEnergy, calculateDepositionEnergyThicken( aveEnergy, self, energyFunc, 
                parameters, energyAccuracy, 1e-10 ) )
        aveMomentum = fudgemath.thickenXYList( aveMomentum, calculateDepositionMomentumThicken( aveMomentum, self, momentumFunc, 
                parameters, momentumAccuracy, 1e-10 ) )

    return( [ aveEnergy ], [ aveMomentum ] )
