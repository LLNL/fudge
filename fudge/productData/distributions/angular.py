# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
The module contains the Angular distribution classes for store :math:`P(\mu)` and :math:`P(\mu|E)` distributions
for two-body and uncorrelated distributions, as well as the :math:`P(\mu|E)` when the
correlated distribution is represented as :math:`P(\mu|E)` * :math:`P(E'|E,\mu)`.

This module contains the following classes:

    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | XYs1d         | Stores :math:`P(\mu)` as a list of :math:`\mu_i`, :math:`P_i`.        |
    +---------------+-----------------------------------------------------------------------+
    | Legendre      | Stores :math:`P(\mu)` as a Legendre series.                           |
    +---------------+-----------------------------------------------------------------------+
    | Xs_pdf_cdf1d  | Stores :math:`P(\mu)` (i.e., :math:`pdf(\mu)`) and :math:`cdf(\mu)`   |
    |               | as a list of :math:`\mu_i`, :math:`pdf_i` and :math:`cdf_i`.          |
    +---------------+-----------------------------------------------------------------------+
    | Isotropic1d   | Specifies that the :math:`P(\mu)` is isotopic.                        |
    +---------------+-----------------------------------------------------------------------+
    | Subform       | Base class for all 2d angular distributions.                          |
    |               | and :py:class:`FSubform` classes.                                     |
    +---------------+-----------------------------------------------------------------------+
    | Isotropic2d   | Specifies that the angular distribution :math:`P(\mu|E)` is isotopic  |
    |               | for all E.                                                            |
    +---------------+-----------------------------------------------------------------------+
    | Forward       | Specifies that the angular distribution :math:`P(\mu|E)` is           |
    |               | :math:`\delta( mu - 1 )` for all E.                                   |
    +---------------+-----------------------------------------------------------------------+
    | Recoil        | Specifies that the angular distribution :math:`P(\mu|E)` is           |
    |               | the recoil of the other two-body outgoing particle.                   |
    +---------------+-----------------------------------------------------------------------+
    | XYs2d         | Stores :math:`P(\mu|E)` as a list of :math:`E`, :math:`P(\mu)`.       |
    +---------------+-----------------------------------------------------------------------+
    | Regions2d     | Stores :math:`P(\mu|E)` as a list of other 2d angular distributions,  |
    |               | for contiguous regions of E.                                          |
    +---------------+-----------------------------------------------------------------------+
    | Form          | This class does not seem to be used anywhere.                         |
    +---------------+-----------------------------------------------------------------------+
    | TwoBody       | Specifies that the :math:`P(\mu|E)` is for a two-body reaction and    |
    | Form          | that the data are in the center-of-mass frame.                        |
    +---------------+-----------------------------------------------------------------------+
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
    """
    Returns an :py:class:`Axes` instance for angular data (i.e., P(mu|E).

    :param energyUnit:  Unit for the energy data.

    :return:            An :py:class:`Axes` instance.
    """

    axes = axesModule.Axes(3)
    axes[0] = axesModule.Axis( 'P(mu|energy_in)', 0, '' )
    axes[1] = axesModule.Axis( 'mu', 1, '' )
    axes[2] = axesModule.Axis( 'energy_in', 2, energyUnit )
    return( axes )


class XYs1d( XYs1dModule.XYs1d ) :
    r"""
    Class to store :math:`P(\mu)` as a list of :math:`\mu_i` and :math:`P_i` points.
    """

    def averageMu( self ) :
        r"""
        This method returns the average value of :math:`\mu` for *self*.
        """

        allowedInterpolations = [xDataEnumsModule.Interpolation.linlin, xDataEnumsModule.Interpolation.flat]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYs1dModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

    def invert( self ) :
        r"""
        This methed returns an :py:class:`XYs1d` instance that is the reflection of *self* (i.e., :math:`P(-\mu)`).
        """

        reflected = self.copy( )
        pdf = [ [ -x, y ] for x, y in self ]
        pdf.reverse( )
        reflected.setData( pdf )
        return( reflected )

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic and False otherwise.
        """

        return( self.rangeMin == self.rangeMax )

    def LegendreCoefficient( self, LegendreOrder ) :
        """
        This method returns the Legenre coefficient representing *self* at Legenre order *LegendreOrder*.

        :param LegendreOrder:   The Legenre order of the desired Legenre coefficient.
        """

        return( LegendreModule.from_pointwiseXY_C( self, LegendreOrder )[LegendreOrder] )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

    def to_xs_pdf_cdf1d(self, style, tempInfo, indent):
        """
        This method returns a copy of *self* as an :py:class:`Xs_pdf_cdf1d` representation.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`Xs_pdf_cdf1d` instance.
        """

        return Xs_pdf_cdf1d.fromXYs(self, outerDomainValue=self.outerDomainValue, thinEpsilon=1e-14)

class Legendre( series1dModule.LegendreSeries ) :
    r"""
    Class to store :math:`P(\mu)` as Legendre series coefficients.
    """

    def averageMu( self ) :
        """
        This function returns the average value of mu for *self*.
        """

        return( self.getCoefficientSafely( 1 ) )

    def LegendreCoefficient( self, LegendreOrder ) :
        """
        This method returns the Legenre coefficient of *self* at Legenre order *LegendreOrder*.

        :param LegendreOrder:   The Legenre order of the desired Legenre coefficient.
        """

        if( LegendreOrder < len( self ) ) : return( self[LegendreOrder] )
        return( 0.0 )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as an :py:class:`Xs_pdf_cdf1d` representation.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`Xs_pdf_cdf1d` instance.
        """

        linear = self.toPointwise_withLinearXYs( accuracy = XYs1dModule.defaultAccuracy, lowerEps = 0, upperEps = 1e-8 )
        linear.outerDomainValue = self.outerDomainValue
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :
    r"""
    This class stores both the :math:`P(\mu|E)` (i.e., pdf) and its cumlative distribution function cdf as three lists: one list
    for the :math:`\mu` values, and one list for the associated pdf value and one list for the associated cdf value for each :math:`\mu` value,
    """

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic and False otherwise.
        """

        return( min( self.pdf.values ) == max( self.pdf.values ) )

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )

class Isotropic1d(ancestryModule.AncestryIO):
    r"""
    This class specifies that the angular distribution :math:`P(\mu)` is isotopic. 
    """

    moniker = 'isotropic1d'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

    @property
    def domainMin( self ) :
        """
        Returns the minimum mu for the product.
        """

        return( -1.0 )

    @property
    def domainMax( self ) :
        """
        Returns the maximum mu for the product.
        """

        return( 1.0 )

    def evaluate( self, mu ) :
        r"""
        This method returns the value of :math:`P(\mu)` which is always 0.5 for an isotropic distribution.

        :param mu:      Value of :math:`\mu`  to evaluated :math:`P(\mu)` at.
        """

        return( 0.5 )

    def isIsotropic( self ) :
        """
        This method always returns True since the distribution is isotropic.
        """

        return( True )

    def integrate(self, domainMin, domainMax):
        """
        This meethod integrates the angular distribution over the specified outgoing mu range.
        
        :param domainMin:               The minimum mu for the range.
        :param domainMax:               The maximum mu for the range.

        :return:                    A float representing the value of the integration.
        """

        return 0.5 * (domainMax - domainMin)

    def LegendreCoefficient( self, LegendreOrder ) :

        if( LegendreOrder == 0 ) : return( 1.0 )
        return( 0.0 )

    def toXML_strList(self, **kwargs):
        """
        This method currently does a **raise** as is not in the GNDS specifications.

        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        raise Exception('Not supported')

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        This should not be called and it executes a **raise** if called as **GNDS** does not define this class.

        :param cls:         Form class to return. 
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        raise Exception('Not supported')

class Subform( baseModule.Subform ) :
    """Abstract base class for 2d angular subforms."""

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns None.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                None
        """

        return( None )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        The method does nothing.
        """

        return 0

class Isotropic2d( Subform ) :
    r"""
    Specifies that the angular distribution :math:`P(\mu|E)` is isotopic for all E. 
    """

    moniker = 'isotropic2d'

    def __init__( self ) :

        Subform.__init__( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        pass

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return( Isotropic2d( ) )

    __copy__ = copy

    @property
    def domainGrid( self ) :
        """
        This method returns the minimum and maximum energy domain for *self* as a Python list.
        """

        return( [ self.domainMin, self.domainMax ] )

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMax )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainUnit )

    def evaluate( self, energy ) :
        """
        This method returns an :py:class:`Isotropic1d` since the distribtion is isotropic at all projectile energies.

        :param energy:  This parameter is not used.
        """

        _isotropic1d = Isotropic1d( )
        _isotropic1d.setAncestor( self )
        return( _isotropic1d )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method just returns [EMin, EMax].

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return( [ EMin, EMax ] )

    def invert( self ) :
        """
        This methed returns another  :py:class:`Isotropic2d` instance since the reflection of an isotropic angular
        distribution is the same as itself.
        """

        return( self.copy( ) )

    def isIsotropic( self ) :
        """
        This method always returns True since the distribution is isotropic.
        """

        return( True )

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E. This
        method always returned 0 since the distribution is isotropic.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """
# FIXME. what frame is mu wanted and what frame are the data in.

        return( 0. )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        return []

    def integrate( self, energyIn, muOut ) :
        """
        This meethod integrates the angular distribution at projectile energy over the specified outgoing mu range.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *muOut* range.
        
        :param energyIn:            The energy of the projectile.
        :param muOut:               The outgoing mu range to integrate over.

        :return:                    A float representing the value of the integration.
        """

        if( self.domainMin <= energyIn <= self.domainMax ) :
            domainMin, domainMax = miscellaneousModule.domainLimits( muOut, -1.0, 1.0 )
            if( domainMax is None ) : return( 1.0 )
            if( domainMax == 0.0 ) : return( 0.5 )
            return( 0.5 * ( domainMax - domainMin ) )

        return( 0.0 )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        ptw = XYs2d( axes = defaultAxes( self.domainUnit ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], outerDomainValue = self.domainMin ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], outerDomainValue = self.domainMax ) )
        return( ptw )

    def toXML_strList( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance *cls*.

        :param cls:         Form class to return. 
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        return cls()

class Forward( Subform ) :
    r"""
    Specifies that the angular distribution :math:`P(\mu|E)` is :math:`\delta( mu - 1 )` for all E. 
    """

    moniker = 'forward'

    def __init__( self ) :

        Subform.__init__( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        pass

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return( Forward( ) )

    __copy__ = copy

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy for *self*.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy for *self*.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainMax )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        from fudge import product as productModule
        return( self.findClassInAncestry( productModule.Product ).domainUnit )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method just returns [EMin, EMax].

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return( [ EMin, EMax ] )

    def invert( self ) :
        r"""
        This methed returns an :py:class:`XYs1d` instance that is the reflection of *self* (i.e., :math:`P(-\mu)`).
        """

        return( self.toPointwise_withLinearXYs( ).invert( ) )

    def isIsotropic( self ) :
        """
        This method always returns False since the distribution is not isotropic.
        """

        return( False )

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E. This
        method always returned 1.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """

        return( 1. )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        return []

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        accuracy = xDataBaseModule.getArguments( kwargs, { 'accuracy' : 1e-6 } )['accuracy']

        ptw = XYs2d(axes=defaultAxes(self.domainUnit))
        ptw.append( XYs1dModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMin ) )
        ptw.append( XYs1dModule.XYs1d( [ [ 1 - accuracy, 0. ], [ 1, 2 / accuracy ] ], outerDomainValue = self.domainMax ) )
        return( ptw )

    def toXML_strList( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance *cls*.

        :param cls:         Form class to return. 
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        return cls()

class Recoil( linkModule.Link, Subform ) :
    r"""
    Specifies that the angular distribution :math:`P(\mu|E)` is the recoil of the other two-body outgoing particle. 
    """

    moniker = 'recoil'

    def __init__( self, link, root = None, path = None, label = None, relative = False ) :
        """
        Distributions for this particle can be computed from its recoil partner using kinematics.
        Only meant for angular 2-body reactions.
        The 'root' and 'path' usually do not need to be specified: they will be computed from the link.

        :param link:        distribution form for recoil partner.
        :type link:         distribution.base.Form.
        :param root:        File path where recoil partner is defined (usually None).
        :param path:        xpath to recoil partner (usually will be computed from the link).
        :param label:       label for this subform (usually None).
        :param relative:    Whether to write xpath as a relative or absolute path.
        :return:
        """

        linkModule.Link.__init__( self, link=link, root=root, path=path, label=label, relative=relative )
        Subform.__init__( self )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        return( self.partner.angularSubform.domainUnit )

    @property
    def domainMin(self):
        """
        Returns the minimum projectile energy for *self*.
        """

        return self.partner.angularSubform.domainMin

    @property
    def domainMax(self):
        """
        Returns the maximum projectile energy for *self*.
        """

        return self.partner.angularSubform.domainMax

    @property
    def partner( self ) :
        """
        This method returns a reference to the angular form for the other outgoing particle of the two-body channel.
        """

        if( self.link is None ) : raise Exception( "Encountered unresolved link!" )
        return( self.link )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        pass

    def copy( self ):
        """
        Returns a copy of *self*.
        """

        return( Recoil( self.partner ) )

    __copy__ = copy

    def evaluate( self, energy ) :
        r"""
        This method returns a :math:`P(\mu)` of *self* evaluated at projectile energy *energy*.

        :param energy:      The projectile energy to evaluate *self* at.
        """

        return( self.partner.angularSubform.evaluate( energy ).invert( ) )

    def getNumericalDistribution( self ) :
        """
        This method returns a pointwise representation of *self*.
        """

        partnerForm = self.partner.toPointwise_withLinearXYs( upperEps = 1e-8 ).invert( )
        return( partnerForm )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* (i.e., its partner) is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return( self.partner.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic and False otherwise.
        """

        return( self.partner.isIsotropic( ) )

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """

        return( -self.partner.averageMu( E, accuracy = accuracy ) )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        if( not( isinstance( self.partner, ( TwoBody, referenceModule.CoulombPlusNuclearElastic ) ) ) ) :
            warnings.append( warning.MissingRecoilDistribution( self ) )

        return warnings

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        return( self.partner.invert( ).toPointwise_withLinearXYs( **kwargs ) )

class XYs2d( Subform, probabilitiesModule.PofX1GivenX2 ) :
    r"""
    Stores :math:`P(\mu|E)` as a list of :math:`E`, :math:`P(\mu)`. The :math:`P(\mu)` functions can be any value 1d angular function.

    :param kwargs:              A dictionary that contains data to control the way this method acts.
    """

    def __init__( self, **kwargs ) :

        probabilitiesModule.PofX1GivenX2.__init__( self, **kwargs )
        Subform.__init__( self )

    def getAtEnergy( self, energy ) :
        r"""
        This method returns :math:`P(\mu)` evaluated at projectile energy *energy*.

        :param energy:      The energy of the projectile.

        :return:            A :math:`P(\mu)`.
        """

        return self.interpolateAtValue(energy, extrapolation=xDataEnumsModule.Extrapolation.flat)

    def fixDomains(self, domainMin, domainMax, fixToDomain, tweakLower = True):
        """
        This method calls *fixDomains* on *self* via the **probabilitiesModule.PofX1GivenX2** class.
        """

        return probabilitiesModule.PofX1GivenX2.fixDomains(self, domainMin, domainMax, fixToDomain, tweakLower = tweakLower)

    def invert( self ) :
        r"""
        This method returns an :py:class:`XYs2d` instance that is the reflection of *self* (i.e., :math:`P(-\mu|E)`).
        """

        ptw = XYs2d( axes = self.axes.copy( ) )
        for POfMu in self : ptw.append( POfMu.invert( ) )
        return( ptw )

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic for all projectile energies and False otherwise.
        """

        for energy_in in self :
            if( not( energy_in.isIsotropic( ) ) ) : return( False )
        return( True )

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """

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
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        Es = [ data.outerDomainValue for data in self ]
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

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
        """
        This method returns a copy of *self* as a :py:class:`XYs2d` instance with the P(mu) data
        returned as :py:class:`Xs_pdf_cdf1d` instances.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`XYs2d` instance.
        """

        subform = XYs2d( axes = self.axes.copy( ), interpolation = self.interpolation )
        for xys in self : subform.append( xys.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( subform )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs1d, Legendre, Xs_pdf_cdf1d ) )

class Regions2d( Subform, regionsModule.Regions2d ) :
    r"""
    Stores :math:`P(\mu|E)` as a list of other 2d angular distributions that abut at specific E values and form
    a contiguous function of E.

    :param kwargs:              A dictionary that contains data to control the way this method acts.
    """

    def __init__( self, **kwargs ) :

        regionsModule.Regions2d.__init__( self, **kwargs )
        Subform.__init__( self )

    def check( self, info ):
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning

        warnings = []
        for idx,region in enumerate(self):
            regionWarnings = region.check( info )
            if regionWarnings: warnings.append( warning.Context('region index %i: %s'
                % (idx, region.moniker), regionWarnings ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

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
        This method call *fixDomains* on the *self* via the **regionsModule.Regions2d** class.
        """

        return regionsModule.Regions2d.fixDomains(self, domainMin, domainMax, fixToDomain)

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic for all regions and False otherwise.
        """

        for region in self :
            if( not( region.isIsotropic( ) ) ) : return( False )
        return( True )

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        Es = [ self[0][0].outerDomainValue ]
        for region in self: Es.extend( [ data.outerDomainValue for data in region[1:] ] )
        if( EMin is not None ) :
            if( EMin < ( 1.0 - 1e-15 ) * Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """

        # FIXME what if E lies on boundary between regions?
        for region in self:
            if( region.domainMax > E ) :
                return region.averageMu( E, accuracy )
        return self[-1].averageMu( E, accuracy )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as a :py:class:`Regions2d` with the P(mu) data
        returned as :py:class:`Xs_pdf_cdf1d` instances.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`Regions2d` instance.
        """

        _regions2d = Regions2d( axes = self.axes )
        for region in self : _regions2d.append( region.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( _regions2d )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

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
        """
        This method returns the class for representing linear point-wise 2-d data of *self*.
        """

        return( XYs2d )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2-d function) of an :py:class:`Regions2d` instance.
        """

        return( ( XYs2d, ) )

"""
# BRB need to implement
    def getAtEnergy( self, energy ) :
    def invert( self ) :
"""

class Form( baseModule.Form ) :
    """
    This class does not seem to be used.
    """

    moniker = 'angular'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        if( not( isinstance( angularSubform, Subform ) ) ) : raise TypeError( 'instance is not an angular subform' )
        baseModule.Form.__init__( self, label, productFrame, ( angularSubform, ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.angularSubform.convertUnits( unitMap )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

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
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        raise NotImplementedError( 'need to implement' )

    @classmethod
    def parseNodeUsingClass(cls, angularElement, xPath, linkData, **kwargs):
        """
        Parse *angularElement* into an instance *cls*.

        :param cls:             Form class to return. 
        :param angularElement:  Node to parse.
        :param xPath:           List containing xPath to current node, useful mostly for debugging.
        :param linkData:        dict that collects unresolved links.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *angularElement*.
        """

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
    r"""
    This class specifies that the :math:`P(\mu|E)` is for a two-body reaction and that the data are in the center-of-mass frame. 

    The following table list the primary members of this class:
    
    +-------------------+-----------------------------------------------------------+           
    | Member            | Description                                               |
    +===================+===========================================================+           
    | angularSubform    | Angular data.                                             |
    +-------------------+-----------------------------------------------------------+           
    
    :param label:               The label for this form.
    :param productFrame:        The frame the product data are specified in.
    :param angularSubform:      An :py:class:`Subform` representing the angular data.
    """

    moniker = 'angularTwoBody'
    subformAttributes = ( 'angularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularSubform ) :

        angularSubform.label = None                 # FIXME: Need to check the type of angularSubform.
        baseModule.Form.__init__( self, label, productFrame, ( angularSubform, ) )

    @property
    def data(self):
        """Return the data of *self*. Temporary method until I redo all the classes that have a data-like member."""

        return(self.angularSubform)

    @data.setter
    def data(self, data):
        r"""
        Sets the data of *self* to *data*. Temporary method until I redo all the classes that have a data-like member. Need to check data type.

        :param data:    The :math:`P(\mu|E)` data for *self.*
        """

        self.angularSubform = data

    @property
    def domainMin(self):
        """
        Returns the minimum projectile energy for *self*.
        """

        return self.angularSubform.domainMin

    def averageMu( self, E, accuracy = None ) :
        """
        This method returns the average value of mu for *self* evaluated at projectile energy E.

        :param E:           The energy of the projectile.
        :param accuracy:    The accurcay to calculate the average mu to.
        """
 
        return( self.angularSubform.averageMu( E, accuracy = accuracy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        kwargs['productFrame'] = self.productFrame
        return( calculateAverageProductData( self.angularSubform, style, indent, **kwargs ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.angularSubform.convertUnits( unitMap )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in.
        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    XYs1d instance for the energy spectrum.
        """

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
        if comKineticEnergy < 0.0:
            return energyModule.XYs1d(axes=energyModule.defaultAxes(energyUnit))
        productComEnergy = comKineticEnergy * residualMass / (productMass + residualMass)

        PofMu = self.angularSubform.evaluate(energyIn)
        if frame == xDataEnumsModule.Frame.centerOfMass:
            twoBodyCOMResolution = min(0.1, max(1e-6, kwargs.get('twoBodyCOMResolution', 1e-2)))
            data = [[productComEnergy* (1 - twoBodyCOMResolution), 0.0],
                    [productComEnergy, 1.0], 
                    [productComEnergy* (1 + twoBodyCOMResolution), 0.0]]
            return energyModule.XYs1d(data, axes=energyModule.defaultAxes(energyUnit)).normalize() * PofMu.integrate(domainMin=muMin, domainMax=muMax)

        projectileSpeed = math.sqrt(2.0 * energyIn / projectileMass)
        comSpeed = projectileMass / (projectileMass + targetMass) * projectileSpeed

        comSpeedProductEnergy = 0.5 * productMass * comSpeed * comSpeed
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
                if mu == 1 and abs(energyPrime - energyIn) < 1e-8 * energyIn:   # Attempt to handle forward elastic scattering exactly.
                    energyPrime = energyIn
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
        r"""
        This method returns the value of *self* at :math:`P(\mu|energy)`.

        :param energy:      This parameter is the projectile energy to evaluate *self* at.
        :param mu:          This parameter is the :math:`\mu` value to evaluate *self* at.
        :param phi:         This parameter is not used.
        :param frame:       The frame to evaluate *self* in.
        """

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
        """
        Need to figure out what this method is doing and document it.
        """

        item = self.angularSubform
        if( isinstance( item, ( linkModule.Link, linkModule.Link2 ) ) ) : links.append( [ item, item.link, item.path ] )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        This method call *fixDomains* on the *angularSubform* member.
        """

        return self.angularSubform.fixDomains(energyMin, energyMax, domainToFix)

    def getEnergyArray( self, EMin = None, EMax = None ) :
        """
        This method returns the list of projectile energies that the data of *self* is specific at.

        :parem EMin:        If *self* does not have an projectile energy values, the default minimum energy to return.
        :parem EMax:        If *self* does not have an projectile energy values, the default maximum energy to return.

        :return:            Python list of float.
        """

        return( self.angularSubform.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param reaction_suite:      The :py:class:`ReactionSuite` instance for this distribution.
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.
        :param muOut:               The outgoing mu range to integrate over.
        :param phiOut:              The outgoing phi range to integrate over.
        :param frame:               The frame the outgoing energy, mu and phi range specify.
        :param LegendreOrder:       The parameter is not used.

        :return:                    A float representing the value of the integration.
        """

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
        r"""
        This method returns a 2d angular instance that is the reflection of *self* (i.e., :math:`P(-\mu|E)`).
        """

        return( self.angularSubform.invert( ) )

    def isIsotropic( self ) :
        """
        This method returns True if *self* is isotropic for all projectile energies and False otherwise.
        """

        return( self.angularSubform.isIsotropic( ) )

    def isTwoBody( self ) :
        """
        This method always returns True since *self* is a two-body distribution.
        """

        return( True )

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an :py:class:`TwoBody` instance representing *self* with (xs, pdf, cdf) data for 
        representing the P(mu) data as needed for Monte Carlo transport.

    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        if( verbosity > 2 ) : print('%sGrouping %s' % ( indent, self.moniker ))

        subform = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( subform is None ) : return( None )

        form = TwoBody( style.label, self.productFrame, subform )
        return( form )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

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
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`XYs2d` instance.
        """

        return( self.angularSubform.toPointwise_withLinearXYs( **kwargs ) )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance *cls*.

        :param cls:         Form class to return. 
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
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
        r"""
        This function returns the tuple (mu1, mu2, isTangent) where mu1 and mu2 are the center-of-
        mass (COM) frame mus for an outgoing particle (product) that correspond to the lab frame mu *muLab*.
        The interaction is two-body, *COM_speed* is the COM mass speed and *productCOM_speed*
        is the product's speed in the COM frame.  If *COM_speed* <= productCOM_speed, 
        there is only one COM mu (mu1) for *muLab*. In this case, the returned tuple is (mu1, None, False).

        Otherwise, there are 0, 1 or 2 solutions as follows. If *muLab* is too large, there is no solution
        and the tuple (None, None, False) is returned. If *muLab* is tangent to the circle of the product's
        velolicty in the COM frame then (mu1, None, True) are returned. Otherwise (mu1, mu2, False) is
        returned with mu1 < mu2.

        :param muLab:               The desired :math:`\mu` value in the lab frame.
        :param COM_speed:           The center-of-mass speed of the system (i.e., projectile and target).
        :param productCOM_speed:    The speed of the outgoing particle in the center-of-mass frame.

        :return:                    (mu1, mu2, isTangent).
        """

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
        """
        This function determines all center-of-mass (COM) frame mu domains (i.e., [mu1, mu2] pairs) that lay 
        between the lab frame mus *muMinLab* and *muMaxLab*. The returned objects is a list that can contain 0, 
        1 or 2 pairs of mu values. Note if *muMinLab* >= *muMaxLab* an empty list is returned. If 
        *COM_speed* < *productCOM_speed* the list [[mu1, mu2]] is returned. Otherwise, there are 3 possible solutions:
        1) no COM mus corresonding to *muMaxLab* (and hence *muMinLab*) in which case an empty list is also returned, 
        2) only *muMaxLab* has corresonding COM mus and [[mu1, mu2]] is returned where mu1 < mu2 or 3) both *muMinLab*
        and *muMaxLab* have corresonding COM mus and [[mu1, mu2], [mu3, mu4]] where mu1 < mu2 < mu3 < mu4 is returned.
        """

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
    This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
    Average means over all outgoing angle and energy.
    
    :param style:   The style instance which the calculated data will belong to.
    :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
    :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
    """

    def calculateDepositionEnergyAtE( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average energy to the outgoing particle in the lab frame.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        a1x, a2y, Qp = parameters['a1x'], parameters['a2y'], parameters['Qp']
        dE = max( 0., E - Qp )
        return( a1x * E + a2y * dE + 2. * math.sqrt( a1x * a2y * E * dE ) * angularData.averageMu( E, accuracy = energyAccuracy ) )

    def calculateDepositionEnergyAtEForPhotoProjectile( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average energy to the outgoing photon in the lab frame
        when the projectile is a photon.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
        E__E_m2 = E / ( E + mass2 )
        dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
        return( 0.5 * E__E_m2**2 * massx + dE + E__E_m2 * math.sqrt( 2. * dE * massx ) * angularData.averageMu( E, accuracy = energyAccuracy ) )

    def calculateDepositionEnergyAtE_forInfiniteTargetMass( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average energy to the outgoing photon in the lab frame assuming theat the
        target mass can be assumed to be infinite. The returned value is alwasy *E*.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        return( E )

    class CalculateDepositionEnergyThicken :
        """
        For internal use. This class is needed by the :py:func:`fudgemath.thickenXYList` function.
        """

        def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

            self.data = data
            self.angular = angular
            self.func = func
            self.parameters = parameters
            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

        def evaluateAtX( self, E ) :
            """
            Returns the average energy to the outgoing particle for prjectile energy *E*.

            :param E:       Energy of the prjectile.
            """

            return( self.func( self.angular, E, self.parameters ) )

    def calculateDepositionMomentumAtE( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average momentum to the outgoing particle in the lab frame.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        mass1, massx, b1x, a2y, Qp = parameters['m1'], parameters['mx'], parameters['b1x'], parameters['a2y'], parameters['Qp']
        dE = max( 0., E - Qp )
        return( b1x * math.sqrt( 2. * E / mass1 ) + math.sqrt( 2. * massx * a2y * dE ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    def calculateDepositionMomentumAtEForPhotoProjectile( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average energy to the outgoing particle in the lab frame
        when the projectile is a photon.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
        E__E_m2 = E / ( E + mass2 )
        dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
        return( E__E_m2 * massx + math.sqrt( 2. * dE * massx ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    def calculateDepositionMomentumAtE_forInfiniteTargetMass( angularData, E, parameters ) :
        r"""
        For internal use. This function calculates the average energy to the outgoing photon in the lab frame assuming
        that the target mass can be assumed to be infinite.

        :param angularData:     The angular data :math:`P(\mu|E)`.
        :param E:               The projectile energy the average energy is evaluated at.
        :param parameters:      Parameters used to calculate the average energy.
        """

        return( math.sqrt( 2.0 * mass1 * E ) * angularData.averageMu( E, accuracy = momentumAccuracy ) )

    class CalculateDepositionMomentumThicken :
        """
        For internal use. This class is needed by the :py:func:`fudgemath.thickenXYList` function.
        """

        def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

            self.data = data
            self.angular = angular
            self.func = func
            self.parameters = parameters
            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

        def evaluateAtX( self, E ) :
            """
            Returns the average momentum to the outgoing particle for prjectile energy *E*.

            :param E:       Energy of the prjectile.
            """

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
