# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""Angular distribution classes with forms form, twoBodyForm and CoulombElasticForm."""

import math
import miscellaneous
import fudge
from fudge.core.utilities import brb
from fudge.core.math import fudgemath

from pqu import PQU as PQUModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import series1d as series1dModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule
from xData import link as linkModule

from fudge.gnd import abstractClasses as abstractClassesModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule

__metaclass__ = type

class XYs1d( XYsModule.XYs1d ) :

    def averageMu( self ) :

        allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                  standardsModule.interpolation.flatToken ]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations = allowedInterpolations )
        return( xys.integrateWithWeight_x( ) )

    def invert( self ) :

        reflected = self.copy( )
        pdf = [ [ -x, y ] for x, y in self ]
        pdf.reverse( )
        reflected.setData( pdf )
        return( reflected )

    def isIsotropic( self ) :

        return( self.rangeMin( ) == self.rangeMax( ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )


class Legendre( series1dModule.LegendreSeries ) :

    def averageMu( self ) :

        return( self.getCoefficientSafely( 1 ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class subform( baseModule.subform ) :
    """Abstract base class for angular subforms."""

    pass

class isotropic( subform ) :

    moniker = 'isotropic'

    def __init__( self ) :

        subform.__init__( self )

    def copy( self ):

        return( isotropic( ) )

    __copy__ = __deepcopy__ = copy

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMin( asPQU=True ).getUnitAsString() )

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

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        ptw = XYs2d( axes = XYs2d.defaultAxes( energyUnit = self.domainUnit( ) ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], value = self.domainMin( ), accuracy = 1e-6 ) )
        ptw.append( XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ], value = self.domainMax( ), accuracy = 1e-6 ) )
        return( ptw )

    def toXMLList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        return( isotropic( ) )

class forward( subform ) :

    moniker = 'forward'

    def __init__( self ) :

        subform.__init__( self )

    def copy( self ):

        return forward()

    __copy__ = __deepcopy__ = copy

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainUnit( ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def invert( self ) :

        return( self.toPointwise_withLinearXYs( invert = True ) )

    def isIsotropic( self ) :

        return( False )

    def averageMu( self, E, accuracy = None ) :

        return( 1. )

    def check( self, info ) :

        return []

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0, invert = False ) :

        eps = 1e-6
        ptw = XYs2d( XYs2d.defaultAxes( energyUnit = self.domainUnit( ) ) )
        if( invert ) :
            ptw.append( XYsModule.XYs1d( [ [ -1, 2 / eps ], [ eps - 1, 0. ] ], value = self.domainMin( ), accuracy = 1e-6 ) )
            ptw.append( XYsModule.XYs1d( [ [ -1, 2 / eps ], [ eps - 1, 0. ] ], value = self.domainMax( ), accuracy = 1e-6 ) )
        else :
            ptw.append( XYsModule.XYs1d( [ [ 1 - eps, 0. ], [ 1, 2 / eps ] ], value = self.domainMin( ), accuracy = 1e-6 ) )
            ptw.append( XYsModule.XYs1d( [ [ 1 - eps, 0. ], [ 1, 2 / eps ] ], value = self.domainMax( ), accuracy = 1e-6 ) )
        return( ptw )

    def toXMLList( self, indent = "", **kwargs ) :

        return( [ self.XMLStartTagString( indent = indent, emptyTag = True ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        raise 'hell'

class recoil( linkModule.link, subform ) :

    moniker = 'recoil'

    def __init__( self, link, root = None, path = None, relative = False ) :
        """
        Distributions for this particle can be computed from its recoil partner using kinematics.
        Only meant for angular 2-body reactions.
        The 'root' and 'path' usually do not need to be specified: they will be computed from the link.

        :param link: distribution form for recoil partner
        :type link: distribution.base.form
        :param root: filename where recoil partner is defined (usually None)
        :param path: xpath to recoil partner (usually will be computed from the link)
        :param relative: whether to write xpath as relative or absolute
        :return:
        """

        linkModule.link.__init__( self, link=link, root=root, path=path, relative=relative )
        subform.__init__( self )

    @property
    def partner( self ) :

        if( self.link is None ) : raise Exception( "Encountered unresolved link!" )
        return( self.link )

    def copy( self ):

        return recoil( self.partner )

    __copy__ = __deepcopy__ = copy

    def getNumericalDistribution( self ) :

        partnerForm = self.partner.invert( )
        return( partnerForm )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.partner.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def isIsotropic( self ) :

        return( self.partner.isIsotropic( ) )

    def averageMu( self, E, accuracy = None ) :

        return( -self.partner.averageMu( E, accuracy = accuracy ) )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if not isinstance( self.partner, twoBodyForm ) :
            warnings.append( warning.missingRecoilDistribution( self ) )

        return warnings

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.partner.invert( ).toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

class XYs2d( subform, multiD_XYsModule.XYs2d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )
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

        mode, function1, function2, frac = self.getBoundingSubFunctions( E )
        if( mode is None ) : raise Exception( 'angular has no data' )
        mu1 = function1.normalize( ).averageMu( )
        if( mode == '' ) :
            mu2 = function2.normalize( ).averageMu( )
            fraction = ( E - function1.value ) / ( function2.value - function1.value )
            mu1 = ( 1 - fraction ) * mu1 + fraction * mu2
            if( self.interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
                raise TypeError( 'Unsupported interpolation qualifier "%s"' % self.interpolationQualifier )
            elif( self.interpolationQualifier != standardsModule.interpolation.noneQualifierToken ) :
                raise TypeError( 'Unsupported interpolation qualifier "%s"' % self.interpolationQualifier )
        return( mu1 )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.value for data in self ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        for idx,function in enumerate(self):
            xys = function.toPointwise_withLinearXYs(1e-6)

            integral = xys.integrate(-1,1)
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQUModule.PQU(function.value,self.axes[-1].unit), idx,
                    integral, function ) )

            if xys.rangeMin() < 0.0:
                warnings.append( warning.negativeProbability( PQUModule.PQU(function.value,self.axes[-1].unit),
                    value = xys.rangeMin(), obj=function ) )

        return warnings

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, accuracy = accuracy, lowerEps = lowerEps,
            upperEps = upperEps, cls = XYs2d ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, Legendre ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', probabilityUnit = '', asLegendre = False ) :

        axes = axesModule.axes( rank = 3 )
        if( asLegendre ) :
            axes[0] = axesModule.axis( 'C_l(energy_in)', 0, '' )
            axes[1] = axesModule.axis( 'l', 1, '' )
        else :
            axes[0] = axesModule.axis( 'P(mu|energy_in)', 0, probabilityUnit )
            axes[1] = axesModule.axis( 'mu', 1, '' )
        axes[2] = axesModule.axis( 'energy_in', 2, energyUnit )
        return( axes )

class regions2d( subform, regionsModule.regions2d ) :

    def __init__( self, **kwargs ) :

        regionsModule.regions2d.__init__( self, **kwargs )
        subform.__init__( self )

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        for idx,region in enumerate(self):
            regionWarnings = region.check( info )
            if regionWarnings: warnings.append( warning.context('region index %i: %s'
                % (idx, region.moniker), regionWarnings ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        verbosity = kwargs.get( 'verbosity', 0 )
        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        EMax = kwargs['EMax']
        aveEnergies = []
        aveMomenta = []
        for i1, region in enumerate( self ) :
            kwargs['EMax'] = region.domainMax( )
            if( i1 == ( len( self ) - 1 ) ) : kwargs['EMax'] = EMax
            aveEnergy, aveMomentum = calculateAverageProductData( region, style, indent, **kwargs ) 
            kwargs['EMin'] = kwargs['EMax']
            aveEnergies.append( aveEnergy[0] )
            aveMomentum.append( aveMomentum[0] )
        return( aveEnergies, aveMomenta )

    def isIsotropic( self ) :

        for region in self :
            if( not( region.isIsotropic( ) ) ) : return( False )
        return( True )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ self[0][0].value ]
        for region in self: Es.extend( [ data.value for data in region[1:] ] )
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def averageMu( self, E, accuracy = None ) :

        # FIXME what if E lies on boundary between regions?
        for region in self:
            if region.domainMax() > E:
                return region.averageMu( E, accuracy )
        return self[-1].averageMu( E, accuracy )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        for region in self :
            if( region.interpolation != standardsModule.interpolation.linlinToken ) :
                print ( '    WARNING: ignoring interpolation "%s" for toPointwise_withLinearXYs' % region.interpolation )
        if( ( lowerEps < 0 ) or ( upperEps < 0 ) ) :
            raise ValueError( 'lowerEps = %s and upperEps = %s must be >= 0.' % ( lowerEps, upperEps ) )
        if( lowerEps == upperEps == 0 ) : raise ValueError( 'lowerEps and upperEps cannot both be 0.' )

        minEps = 5e-16
        if( lowerEps != 0 ) : lowerEps = max( minEps, lowerEps )
        if( upperEps != 0 ) : upperEps = max( minEps, upperEps )

        axes = XYs2d.defaultAxes(  )
        axes[2].unit = self.axes[2].unit
        linear = self.toLinearXYsClass( )( axes = axes )
        n1 = len( self )
        for i1, region in enumerate( self ) :
            _region = region.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            if( i1 < n1 - 1 ) :
                if( lowerEps > 0  ) : _region[-1].value = _region[-1].value * ( 1 - lowerEps )
            if( i1 > 0 ) :
                if( upperEps > 0  ) : _region[0].value = _region[0].value * ( 1 + upperEps )
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
    """
    This class only seems to be used in site_packages/legacy/toENDF6/productData/distributions/uncorrelated.py
    to write GND structure to ENDF6. This class should be deprecated.
    """

    moniker = 'angular'
    subformAttributes = ( 'angularSubform', )

    def __init__( self, label, productFrame, angularSubform, makeCopy = True ) :

        if( not( isinstance( angularSubform, subform ) ) ) : raise TypeError( 'instance is not an angular subform' )
        if( makeCopy ) : angularSubform = angularSubform.copy( )

        baseModule.form.__init__( self, label, productFrame, ( angularSubform, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise 'FIXME, when am I called'

        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        product = kwargs['product']
        reactionSuite = kwargs['reactionSuite']

        if( product.name != 'gamma' ) :
            raise Exception( 'For form %s, calculateAverageProductData is only for gammas, not %s' % ( self.moniker, product.name ) )

        depData = []
        angularSubform = self.angularSubform

        Es = angularSubform.getEnergyArray( kwargs['EMin'], kwargs['EMax'] )
        if( ( 'discrete' in product.attributes ) or ( 'primary' in product.attributes ) ) :
            massRatio = 0.
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
            else :
                Eg = product.attributes['primary'].getValueAs( energyUnit )
                mass1 = reactionSuite.projectile.getMass( massUnit )
                mass2 = reactionSuite.target.getMass( massUnit )
                massRatio = mass2 / ( mass1 + mass2 )
            depEnergy = [ [ E, Eg + massRatio * E ] for E in Es ]
            depMomentum = [ [ E, ( Eg + massRatio * E ) * angularSubform.averageMu( E, accuracy = 0.1 * momentumAccuracy ) ] for E in Es ]
        else :
            raise Exception( 'Unsupported gamma; gamma must be "discrete" or "primary"' )

        axes = energyDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( energyDepositionModule.XYs1d( data = depEnergy, axes = axes,
                label = style.label, accuracy = energyAccuracy ) )
        axes = momentumDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( momentumDepositionModule.XYs1d( data = depMomentum, axes = axes,
                label = style.label, accuracy = momentumAccuracy ) )

        return( depData )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        raise Exception( 'need to implement' )

#    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.gnd import miscellaneous as miscellaneousModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        product = tempInfo['product']
        if( product.name != 'gamma' ) : raise Exception( 'For form %s, process is only for gammas, not %s' % ( self.moniker, product.name ) )

        newForms = []
        angularSubform = self.angularSubform
        crossSection = tempInfo['crossSection']
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            projectileName, productName = processInfo.getProjectileName( ), product.name
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
                TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( processInfo, projectileName, productName, Eg, crossSection,
                   angularSubform, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            elif( 'primary' in product.attributes ) :
                projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
                mass1, mass2 = projectile.getMass( massUnit ), target.getMass( massUnit )
                massRatio = mass2 / ( mass1 + mass2 )
                ELevel = product.attributes['primary'].getValueAs( energyUnit )
                TM_1, TM_E = transferMatricesModule.primaryGammaAngularData( processInfo, projectileName, productName, ELevel, massRatio, crossSection,
                   angularSubform, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            else :
                raise Exception( 'Unsupported gamma; gamma must be "discrete" or "primary"' )
            miscellaneousModule.TMs2Form( processInfo, tempInfo, newForms, TM_1, TM_E, crossSection.axes )

        return( newForms )

    @staticmethod
    def parseXMLNode( angularElement, xPath, linkData ) :
        """Translate <angular> element from xml."""

        xPath.append( angularElement.tag )
        subformClass = {
                isotropic.moniker   : isotropic,
                forward.moniker     : forward,
                XYs2d.moniker       : XYs2d,
                regions2d.moniker   : regions2d,
            }.get( angularElement.tag )
        if( subformClass is None ) : raise Exception( "encountered unknown angular subform: %s" % angularElement.tag )
        angularSubform = subformClass.parseXMLNode( angularElement, xPath, linkData )
        xPath.pop( )
        return( angularSubform )

class twoBodyForm( baseModule.form ) :

    moniker = 'angularTwoBody'
    subformAttributes = ( 'angularSubform', )

    def __init__( self, label, productFrame, angularSubform, makeCopy = True ) :

        angularSubform.label = None
        if( makeCopy ) : angularSubform = angularSubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, ) )

    def averageMu( self, E, accuracy = None ) :
 
        return( self.angularSubform.averageMu( E, accuracy = accuracy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        kwargs['productFrame'] = self.productFrame
        return( calculateAverageProductData( self.angularSubform, style, indent, **kwargs ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.angularSubform.getEnergyArray( EMin = EMin, EMax = EMax ) )

    def invert( self ) :

        return( self.angularSubform.invert( ) )

    def isIsotropic( self ) :

        return( self.angularSubform.isIsotropic( ) )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']
        outputChannel = tempInfo['reaction'].outputChannel

        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )

        crossSection = tempInfo['crossSection']
        angularSubform = self.angularSubform
        if( isinstance( angularSubform, recoil ) ) :
            if( isinstance( angularSubform.link, CoulombExpansionForm ) ) : return( None )
            angularSubform = angularSubform.getNumericalDistribution( )

        Q = outputChannel.Q['eval']
        Q = Q.getValue( Q.domainMin( ), unit = tempInfo['incidentEnergyUnit'] )
        residual = outputChannel.products[1]
        if( tempInfo['productIndex'] == '1' ) : residual = outputChannel.products[0]
        residualMass = tempInfo['masses']['Residual']
        tempInfo['masses']['Residual'] = residual.getMass( tempInfo['massUnit'] )
        TM_1, TM_E = transferMatricesModule.twoBodyTransferMatrix( style, tempInfo, self.productFrame, crossSection, angularSubform, 
                Q, comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
        tempInfo['masses']['Residual'] = residualMass
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.angularSubform.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    @staticmethod
    def parseXMLNode( angularElement, xPath, linkData ) :
        """
        Translate <angularTwoBody> element from xml. Returns instance of angular.twoBodyForm.
        """

        xPath.append( angularElement.tag )

        subformElement = angularElement[0]
        subformClass = {    isotropic.moniker   : isotropic,
                            recoil.moniker      : recoil,
                            XYs2d.moniker       : XYs2d,
                            regions2d.moniker   : regions2d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )
        angularForm = twoBodyForm( angularElement.get( 'label' ), angularElement.get( 'productFrame' ),
                                             angularSubform, makeCopy=False )
        xPath.pop( )
        return( angularForm )

class CoulombExpansionForm( baseModule.form ) :
    """Charged-particle Coulomb scattering, stored as three separate terms (nuclear + real/imaginary interference)"""

    moniker = 'CoulombExpansion'
    subformAttributes = ( 'nuclear_term', 'interferenceReal_term', 'interferenceImaginary_term' )

    def __init__( self, label, productFrame, nuclear_term, interferenceReal_term, interferenceImaginary_term, makeCopy = True ) :

        if( makeCopy ) :
            nuclear_term = nuclear_term.copy()
            interferenceReal_term = interferenceReal_term.copy()
            interferenceImaginary_term = interferenceImaginary_term.copy()
        if( isinstance( nuclear_term, XYs2d ) ) :
            assert isinstance( interferenceReal_term, XYs2d ) and isinstance( interferenceImaginary_term, XYs2d )
        elif( isinstance( nuclear_term, regions2d) ) :
            assert isinstance( interferenceReal_term, regions2d ) and isinstance( interferenceImaginary_term, regions2d )
        else :
            raise Exception( 'Only XYs2d and regions2d are currently supported, not %s' % brb.getType( nuclear_term ) )
        baseModule.form.__init__( self, label, productFrame, ( nuclear_term, interferenceReal_term, interferenceImaginary_term ) )
        self.nuclear_term.label = 'nuclear_term'
        self.interferenceReal_term.label = 'interferenceReal_term'
        self.interferenceImaginary_term.label = 'interferenceImaginary_term'

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        print '    processSnMultiGroup not implemented for distribution form %s.' % self.moniker
        return( None )
    
    def toXMLList( self, indent = "", **kwargs ):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s" productFrame="%s">' % ( indent, self.moniker, self.label, self.productFrame ) ]
        for term in ( 'nuclear_term', 'interferenceReal_term', 'interferenceImaginary_term' ) :
            xmlString += getattr( self, term ).toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        kwargs = {}
        for subformElement in element:
            subformClass = {    XYs2d.moniker       : XYs2d,
                                regions2d.moniker   : regions2d,
                }.get( subformElement.tag )
            if( subformClass is None ) : raise Exception( "encountered unknown CoulombElastic subform: %s" % subformElement.tag )
            kwargs[ subformElement.get('label') ] = subformClass.parseXMLNode( subformElement, xPath, linkData )
        kwargs['makeCopy'] = False
        form = CoulombExpansionForm( element.get( 'label' ), element.get('productFrame'), **kwargs )
        xPath.pop()
        return( form )

class NuclearPlusCoulombInterferenceForm( baseModule.form ) :
    """Charged-particle elastic scattering (not including pure Coulomb portion), stored as a single term (nuclear + interference)"""

    moniker = 'NuclearPlusCoulombInterference'
    subformAttributes = ( 'angularSubform', )

    def __init__( self, label, productFrame, angularSubform, makeCopy = True ) :

        if( not( isinstance( angularSubform, ( XYs2d, regions2d ) ) ) ) : raise TypeError( 'instance is not an angular subform' )
        if( makeCopy ) : angularSubform = angularSubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise CoulombDepositionNotSupported( "Cannot compute average product data for %s distribution" % self.moniker )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        print '    processSnMultiGroup not implemented for distribution form %s.' % self.moniker
        return( None )
    
    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s" productFrame="%s">' % ( indent, self.moniker, self.label, self.productFrame ) ]
        xmlString += self.angularSubform.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        # nuclearPlusCoulombInterference is divided into three sections: nuclear, real and imag. interference
        # each may be LegendrePointwise or LegendrePiecewise

        xPath.append( element.tag )
        angularSubformClass = {
            XYs2d.moniker : XYs2d 
        }.get( element[0].tag )
        angularSubform = angularSubformClass.parseXMLNode( element[0], xPath, linkData )
        npci = NuclearPlusCoulombInterferenceForm( element.get( 'label' ), element.get('productFrame'),
                                                        angularSubform = angularSubform, makeCopy=False )
        xPath.pop()
        return npci

class CoulombDepositionNotSupported( Exception ):
    """ Custom Exception, returned when calculatedDepositionData() called for CoulombExpansion or NuclearPlusCoulomb """
    pass

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
    momentumDepositionUnit = energyUnit + '/c'
    massUnit = energyUnit + '/c**2'
    multiplicity = kwargs['multiplicity']
    productMass = kwargs['product'].getMass( massUnit )
    energyAccuracy = kwargs['energyAccuracy']
    momentumAccuracy = kwargs['momentumAccuracy']
    reactionSuite = kwargs['reactionSuite']
    reaction = kwargs['reaction']
    outputChannel = kwargs['outputChannel']
    productIndex = kwargs['productIndex']

    mass1 = reactionSuite.projectile.getMass( massUnit )
    mass2 = reactionSuite.target.getMass( massUnit )
    massx = outputChannel.products[0].getMass( massUnit )
    massy = outputChannel.products[1].getMass( massUnit )

    if( productIndex == '1' ) : massx, massy = massy, massx
    m12, mxy = mass1 + mass2, massx + massy
    b1x, a1x, a2y = mass1 * massx / m12, mass1 * massx / ( m12 * m12 ), mass2 * massy / ( m12 * mxy )
    Qm = m12 - mxy
    Q = reaction.getQ( energyUnit, final = False, groundStateQ = True )
    Qp = -float( Q ) * m12 / mass2              # This is the threshold in the COM frame.

    if( mass1 == 0. ) :                         # Photo as projectile
        energyFunc = calculateDepositionEnergyAtEForPhotoProjectile
        momentumFunc = calculateDepositionMomentumAtEForPhotoProjectile
        parameters = { 'm2' : mass2, 'mx' : massx, 'my' : massy, 'Q' : Q }
    else :
        energyFunc = calculateDepositionEnergyAtE
        momentumFunc = calculateDepositionMomentumAtE
        parameters = { 'm1' : mass1, 'mx' : massx, 'a1x' : a1x, 'a2y' : a2y, 'b1x' : b1x, 'Qp' : Qp }

    Es = self.getEnergyArray( kwargs['EMin'], kwargs['EMax'] )
    aveEnergy = [ [ E, energyFunc( self, E, parameters ) ] for E in Es ]
    aveEnergy = fudgemath.thickenXYList( aveEnergy, calculateDepositionEnergyThicken( aveEnergy, self, energyFunc, 
            parameters, energyAccuracy, 1e-10 ) )

    aveMomentum = [ [ E, momentumFunc( self, E, parameters ) ] for E in Es ]
    aveMomentum = fudgemath.thickenXYList( aveMomentum, calculateDepositionMomentumThicken( aveMomentum, self, momentumFunc, 
            parameters, momentumAccuracy, 1e-10 ) )

    return( [ aveEnergy ], [ aveMomentum ] )
