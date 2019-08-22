# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import math

from pqu import PQU 

from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import regions as regionsModule
from xData import link as linkModule
from xData import ancestry as ancestryModule

from fudge.gnd import abstractClasses as abstractClassesModule
from . import base as baseModule
from fudge.gnd import tokens as tokensModule
from fudge.gnd import styles as stylesModule

__metaclass__ = type

lowerEps = 1e-8
upperEps = 1e-8

#
# crossSection forms.
#
class baseCrossSectionForm( abstractClassesModule.form ) :

    __genre = baseModule.crossSectionToken

class pointwise( baseCrossSectionForm, XYsModule.XYs ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYsModule.XYs.__init__( self, **kwargs )

    def heat( self, ancestors, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False, 
            doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) :
        """
        Returns a linear version of the cross section heated to 'temperature'. If the current temperature of the cross section
        is greater than 'temperature' a raise is executed.
        If lowerlimit is None, it is set to 'oneOverV' except when the reaction is determined to be a threshold reaction, 
        then it is set to 'threshold'. Any cross section with domainMin greater than 2.5e-4 eV is determined to be a
        threshold reaction. If upperlimit is None it is set to 'constant'.
        If heatBelowThreshold is False, then EMin is set to the larger of EMin and self's domainMin.

        For more information on EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin and heatAllEDomain
        see the module crossSectionAdjustForHeatedTarget.
        """

        from crossSectionAdjustForHeatedTarget import heat
        from pqu import PQU

        gotT = currentTemperature = ancestors.findAttributeInAncestry( 'getTemperature' )( self.label )
        if( gotT.isTemperature( ) ) : gotT = PQU.PQU( gotT.getValueAs( 'K' ), 'K' ) * PQU.PQU( '1 k' )
        gotT = gotT.getValueAs( self.domainUnit( ) )

        if( isinstance( temperature, str ) ) : temperature = PQU.PQU( temperature )
        temperature_ = temperature
        if( temperature.isTemperature( ) ) :
            temperature_ = PQU.PQU( temperature.getValueAs( 'K' ), 'K' ) * PQU.PQU( '1 k' )
        wantedT = temperature_.getValueAs( self.domainUnit( ) )

        dT = wantedT - gotT
        if( abs( dT ) <= 1e-2 * wantedT ) : dT = 0.
        if( dT < 0 ) : raise Exception( 'Current temperature "%s" (%.4e) higher than desired temperature "%s" (%.4e)' % ( currentTemperature, gotT, temperature, wantedT ) ) 

        projectile, target = ancestors.findAttributeInAncestry( 'projectile' ), ancestors.findAttributeInAncestry( 'target' )
        if( projectile.getMass( 'amu' ) == 0 ) : raise Exception( 'Heating with gamma as projectile not supported.' )
        massRatio = target.getMass( 'amu' ) / projectile.getMass( 'amu' )

        heated = unheated = self
        if( not( unheated.isInterpolationLinear( ) ) ) : unheated = self.toPointwise_withLinearXYs( lowerEps, upperEps )
        if( ( len( unheated ) > 0 ) and ( dT > 0. ) ) :
            if( isinstance( EMin, str ) ) : EMin = PQU.PQU( EMin )
            EMin_ = EMin
            if( not( isinstance( EMin, float ) ) ) : EMin_ = EMin.getValueAs( self.domainUnit( ) )
            if( not( heatBelowThreshold ) ) : EMin_ = max( unheated.domainMin( ), EMin_ )
            if( lowerlimit is None ) :
                lowerlimit = "oneOverV"
                domainMin = PQU.PQU( 1e-5, 'eV' ).getValueAs( self.domainUnit( ) )
                if( domainMin < 0.04 * unheated.domainMin( ) ) : lowerlimit = "threshold"
            if( upperlimit is None ) : upperlimit = 'constant'

            heated = heat.crossSectionAdjustForHeatedTarget( massRatio, dT, EMin_, unheated, lowerlimit = lowerlimit, \
                upperlimit = upperlimit, interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                heatAllEDomain = heatAllEDomain )
        accuracy = max( interpolationAccuracy, self.getAccuracy( ) )
        return( pointwise( data = heated, axes = unheated.axes, accuracy = accuracy ) )

    def process( self, component_, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing import miscellaneousModule

        if( processInfo.verbosity >= 20 ) : print '%sprocessing pointwise cross section' % verbosityIndent
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        forms = []

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            crossSection = miscellaneousModule.groupOneFunctionAndFlux( projectile, processInfo, self )
            tempInfo['groupedCrossSectionNorm'] = crossSection

            groupedFlux = tempInfo['groupedFlux']
            crossSectionGrouped_ = [ crossSection[i] / groupedFlux[i] for i in xrange( len( crossSection ) ) ]

            axes = self.axes.copy( standAlone = True )
            forms.append( grouped( axes, crossSectionGrouped_ ) )

        return( forms )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        ptw = XYsModule.XYs.toPointwise_withLinearXYs( self, lowerEps = lowerEps, upperEps = upperEps, cls = pointwise )
        ptw.label = self.label
        return( ptw )

    def toXMLList( self, indent = "", **kwargs ) :

        return( XYsModule.XYs.toXMLList( self, indent, **kwargs ) )

class piecewise( baseCrossSectionForm, regionsModule.regions ) :
    """
    This class stores a cross section in two or more regions, which may have different interpolations.
    Each region must contain at least two points. Each pair of adjacent regions must overlap at exactly one point.
    """

    def __init__( self, **kwargs ) :

        if( 'dimension' not in kwargs ) : kwargs['dimension'] = 1
        if( kwargs['dimension'] != 1 ) : raise ValueError( 'Dimension = %s != 1' % ( kwargs['dimension'] ) )
        regionsModule.regions.__init__( self, **kwargs )

    def toLinearXYsClass( self ) :

        return( pointwise )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        accuracy = 1e-3
        if( len( self ) > 0 ) : accuracy = self[0].getAccuracy( )
        xys = regionsModule.regions.toPointwise_withLinearXYs( self, accuracy, lowerEps, upperEps )
        return( pointwise( data = xys, axes = xys.axes, accuracy = xys.getAccuracy( ) ) )

class grouped( abstractClassesModule.multiGroup ) :

    __genre = baseModule.crossSectionToken

    def __init__( self, axes, groupData ) :

        abstractClassesModule.multiGroup.__init__( self, axes, groupData )

class resonanceLink( linkModule.link ) :

    moniker = 'resonanceRegion'
    __genre = baseModule.crossSectionToken

class resonancesWithBackground( baseCrossSectionForm ) :
    """
    This class stores cross sections that include resonances along with a background contribution.
    Contains a link to the resonances, and the 'background' which could be pointwise or piecewise data.
    The full pointwise cross section can be obtained by first reconstructing resonances
    and then adding the background contribution (users should use the reactionSuite.reconstructResonances method).
    """

    moniker = tokensModule.resonancesWithBackgroundFormToken

    def __init__( self, label, tabulatedData, resonanceLink ) :
        
        self.link = resonanceLink
        self.tabulatedData = tabulatedData
        self.tabulatedData.attribute = self.tabulatedData.label = None
        self.tabulatedData.ancestor = self

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def findEntity( self, entityName, attribute = None, value = None ):

        if entityName in (pointwise.moniker, piecewise.moniker):
            return self.tabulatedData
        return ancestryModule.ancestry.findEntity( self, entityName, attribute, value )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.domain( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.tabulatedData.domainUnit( ) )

    @property
    def label( self ) :

        return( self.__label )
    
    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        _component = self.findClassInAncestry( component )
        reconstructed = _component.getStyleOfClass( stylesModule.crossSectionReconstructed )
        if( reconstructed is not None ) :
            return( reconstructed.toPointwise_withLinearXYs( lowerEps, upperEps ) )     # Return a copy.
        raise Exception( 'resonancesWithBackground cross section has not been reconstructed via reactionSuite.reconstructResonances' )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        xmlList = [ indent + '<%s%s>' % ( self.moniker, attributeStr ), self.link.toXML( indent2, **kwargs )  ]
        xmlList += self.tabulatedData.toXMLList( indent2, **kwargs )
        xmlList[-1] += '</%s>' % self.moniker
        return( xmlList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        xlink = resonanceLink.parseXMLNode( element.find('resonanceRegion') )
        form = element[1]
        formClass = {pointwise.moniker: pointwise,
                piecewise.moniker: piecewise,
                weightedPointwise.moniker: weightedPointwise,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown cross section form: %s" % form.tag)
        resWithBack = resonancesWithBackground( element.get( 'label' ), formClass.parseXMLNode(form,xPath,linkData), xlink )
        xPath.pop()
        return resWithBack
    
class reference( linkModule.link, baseCrossSectionForm ) :
    """This cross section form consists of a reference to another cross section."""

    moniker = tokensModule.referenceFormToken

    def __init__( self, link = None, root = None, path = None, label = None ) :

        linkModule.link.__init__(self, link = link, root = root, path = path, label = label )
        baseCrossSectionForm.__init__( self )

    @property
    def crossSection( self ) :

        if self.link is None : raise Exception( "Unresolved link!" )
        return self.link

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domain( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.crossSection.domainUnit( ) )

    def getReference( self ) :

        return( self.crossSection )

    def setReference( self, crossSection ) :

        if( not( isinstance( crossSection, (component, None.__class__) ) ) ) :
            raise Exception( 'crossSection argument must be a cross section component not type %s' % type( crossSection ) )
        self.crossSection = crossSection

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        return( self.crossSection.toPointwise_withLinearXYs( lowerEps, upperEps ) )

    def check( self ) :

        from fudge.gnd import warning
        warnings = []
        if self.getRootAncestor() != self.getReference().getRootAncestor():
            warnings.append( warning.badCrossSectionReference() )
        return warnings

class weightedPointwise( baseCrossSectionForm ) :
    """
    This class is used to store one cross section as a weighted multiple of another cross section.
    This is typically used to express the 'production' cross section for radioactive products of a reaction.
    The production cross section can be given as an energy-dependent multiple of the reaction cross section.

    The weightedPointwise class contains a link to another cross section, and also contains a set of
    energy-dependent weights that should be multiplied by the other cross section.
    """

    moniker = tokensModule.weightedPointwiseFormToken

    def __init__( self, label, crossSection, weights ) :

        self.setReference( crossSection )
        self.weights = weights

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.weights.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.weights.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.weights.domain( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.crossSection.domainUnit( ) )

    def setReference( self, crossSection ):

        self.crossSection = crossSection

    @property
    def label( self ) :

        return( self.__label )
    
    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        xys = self.crossSection.toPointwise_withLinearXYs( lowerEps, upperEps ) * self.weights
        return( pointwise( data = xys, axes = xys.axes, accuracy = xys.getAccuracy( ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        xmlList = [ '%s<%s%s xlink:type="simple" xlink:href="%s">' % \
                ( indent, attributeStr, self.form, self.crossSection.toXLink( ) ) ]
        xmlList += self.weights.toXMLList( indent2, **kwargs )
        xmlList[-1] += '</%s>' % self.form
        return( xmlList )

    @staticmethod
    def parseXMLNode( WPelement, xPath, linkData, **kwargs ) :

        xlink = linkModule.link.parseXMLNode(WPelement).path
        weights = XYsModule.XYs.parseXMLNode( WPelement[0], xPath, linkData )
        ref = weightedPointwise( None, weights )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

def defaultAxes( energyUnit = 'eV', crossSectionUnit = 'b' ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( baseModule.crossSectionToken, 0, crossSectionUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

def chargeParticle_changeInterpolationSubFunction( threshold, x, x1, y1, x2, y2 ) :

    B = math.log( x2 * y2 / ( x1 * y1 ) ) / ( 1. / math.sqrt( x1 - threshold ) - 1. / math.sqrt( x2 - threshold ) )
    A = x1 * y1 * math.exp( B / math.sqrt( x1 - threshold ) )
    y = A * math.exp( - B / math.sqrt( x - threshold ) ) / x
    return( y )

#
# crossSection component
#
class component( abstractClassesModule.component ) :

    moniker = baseModule.crossSectionToken
    __genre = baseModule.crossSectionToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self,
                ( pointwise, piecewise, grouped, resonanceLink, resonancesWithBackground, reference, weightedPointwise ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.evaluated.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.evaluated.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.evaluated.domainUnit( ) )

    def hasLinearForm( self ) :

        for form in self :
            if( isinstance( form, pointwise ) ) :
                if( form.isInterpolationLinear( ) ) : return( form )
        return( None )

    def heat( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, 
            heatAllPoints = False, doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) : 
        """
        Returns the result of self.toPointwise_withLinearXYs( ).heat( ... ). See method pointwise.heat for more information.
        """

        return( self.toPointwise_withLinearXYs( ).heat( self, temperature, EMin, lowerlimit, upperlimit, 
                interpolationAccuracy, heatAllPoints, doNotThin, heatBelowThreshold, heatAllEDomain ) )

    def check( self, info ) :
        """
        Check cross section data for correct threshold, negative cross sections, etc.
        Returns a list of any warnings encountered during checking.
        
        :param dict info:  A Python dictionary containing the parameters that control the cross section checking.
        :keyword boolean CoulombReaction: True if target and projectile are both charged particles, or if two or more products are charged particles.
        :keyword float Q: a parameter of `info`: if Q is positive (and CoulombReaction=False), cross section must start at 1e-5 eV, otherwise, cross section threshold must agree with Q-value
        :keyword float kinematicFactor: a parameter of `info`: equal to ( (targetMass + projectileMass) / targetMass )
        :keyword string dThreshold: a parameter of `info`: allowable threshold energy mismatch as a string suitable for the PQU class
        :keyword float crossSectionEnergyMax: a parameter of `info`: the cross section must extend up to limit (usually 20 MeV)
        """
        from fudge.gnd import warning
        warnings = []

        # FIXME: domain is giving incorrect answers for threshold reactions with resonance contributions!
        lower, upper = self.domain( asPQU = True )
        thresh = PQU.PQU( 0, 'eV' )
        if 'Q' in info:
            thresh = PQU.PQU( -info['Q'] * info['kinematicFactor'], 'eV' )

            if not info['CoulombChannel']:
                # if Q is positive, cross section must start at 1e-5 eV
                # otherwise, cross section threshold must agree with Q-value:
                if thresh.value>0:
                    if abs(thresh-lower) > PQU.PQU( info['dThreshold'] ):
                        warnings.append( warning.threshold_mismatch( lower, thresh, self ) )
                elif lower != PQU.PQU(1e-5,'eV'):
                    # ignore 2nd,3rd,4th-chance fission (they have Q>0 but are still threshold reactions):
                    from fudge.gnd import channels
                    parent = self.getAncestor( )
                    if (not hasattr( parent, 'outputChannel' ) or ( parent.outputChannel.getFissionGenre()
                            in (None, channels.fissionGenreTotal, channels.fissionGenreFirstChance) ) ):
                        warnings.append( warning.threshold_mismatch( lower, PQU.PQU(1e-5, 'eV'), self ) )
            else:
                # charged-particle reaction generally doesn't 'turn on' right at threshold due to Coulomb barrier.
                # In this case, just ensure the cross section stays zero until at or above threshold:
                if (lower < thresh):
                    warnings.append( warning.Coulomb_threshold_mismatch( thresh, lower, self ) )

        # cross section must extend up to limit (usually 20 MeV):
        if upper < PQU.PQU( info['crossSectionEnergyMax'] ):
            warnings.append( warning.gapInCrossSection( upper, PQU.PQU( info['crossSectionEnergyMax'] ),
                self ) )

        evaluatedCrossSection = self.evaluated
        if( isinstance( evaluatedCrossSection, resonancesWithBackground ) ) :
            # ensure no gaps between resonance and background:
            resDomain = info['reactionSuite'].resonances.domain( asPQU = True )
            bckDomain = evaluatedCrossSection.tabulatedData.domain( asPQU = True )
            if bckDomain[0] > resDomain[1]:
                warnings.append( warning.gapInCrossSection(resDomain[1],bckDomain[0], self ) )
            linearCrossSection = self[info['reconstructedStyle']]

        elif( isinstance( evaluatedCrossSection, reference )
                and isinstance( evaluatedCrossSection.link.evaluated, resonancesWithBackground) ):
            linearCrossSection = evaluatedCrossSection.link[info['reconstructedStyle']]

        else:
            # can it be converted to pointwise linear?
            try:
                linearCrossSection = evaluatedCrossSection.toPointwise_withLinearXYs(0,1e-8)
            except Exception as e:
                warnings.append( warning.ExceptionRaised( e, self ) )
                if info['failOnException']: raise
                return warnings

        # test for negative values, and for non-zero cross section at threshold
        if linearCrossSection.rangeMin() < 0:
            for i,(en,xsc) in enumerate(linearCrossSection):
                if xsc < 0:
                    warnings.append( warning.negativeCrossSection( PQU.PQU(en, linearCrossSection.axes[-1].unit), i, self ) )

        if thresh.value>0:
            if linearCrossSection[0][1] != 0:
                warnings.append( warning.nonZero_crossSection_at_threshold( linearCrossSection[0][1], self ) )

        return warnings

    def process( self, processInfo, tempInfo, verbosityIndent ) :
        """
        Process self into a form suitable for a specific application.
        """

        crossSection = self.toPointwise_withLinearXYs( )
        forms = crossSection.process( self, processInfo, tempInfo, verbosityIndent )
        for form in forms : self.add( form )
        return( crossSection )
        
    def findEntity( self, entityName, attribute = None, value = None ):

        for form in self :
            if( form.moniker == entityName ) : return( form )
        return ancestryModule.ancestry.findEntity( self, entityName, attribute, value )

    def evaluateWithUncertainty( self, E, useCovariance = True, covariance = None, covarianceSuite = None ) :
        """

        :param E:
        :param useCovariance: use this to override covarance usage
        :type useCovariance: bool
        :param covariance:
        :param covarianceSuite:
        :return:
        """
        import fudge.gnd.covariances.base as covModule

        if( not isinstance( E, PQU.PQU ) ) : raise TypeError( "E must be an PQU.PQU instance")

        ptwise = self.hasLinearForm()
        if ptwise is None: ptwise = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)

        meanValue = ptwise.evaluate( E.getValueAs( ptwise.domainUnit( ) ) )

        # We might be done
        if not useCovariance: return PQU.PQU( meanValue, unit=ptwise.rangeUnit() )

        # Get that the covariance goes with the data.
        if covariance is None:
            if hasattr(self.evaluated, 'uncertainties') and self.evaluated.uncertainties:
                if hasattr(self.evaluated.uncertainties[0].data,'link'):
                    covariance = self.evaluated.uncertainties[0].data.link.forms['eval']
                elif covarianceSuite is not None:
                    covariance = self.evaluated.uncertainties[0].data.follow(startNode=covarianceSuite).forms['eval']
                else: return PQU.PQU( meanValue, unit = ptwise.rangeUnit() )
            else: return PQU.PQU( meanValue, unit = ptwise.rangeUnit() )
        else:
            if covariance is not( isinstance( covariance, covModule.covarianceMatrix ) ):
                raise TypeError( 'covariance must be of type covarianceMatrix, got %s'%str(type(covariance)))

        # Compute uncertainty on convolution given mean value and covariance.
        try:
            uncertainty = covariance.getUncertaintyVector( self.evaluated, relative = False )\
                                                       .evaluate( E.getValueAs( ptwise.domainUnit( ) ) )
        except:
            print "WARNING: could not get uncertainty in evaluateWithUncertainty for %s" % str(covariance.__class__)
            return PQU.PQU( meanValue, unit = ptwise.rangeUnit() )

        # We haven't changed units on uncertainty or the meanValue, so we can reconstruct the physical quantity easily.
        return( PQU.PQU( meanValue, unit = ptwise.rangeUnit(), uncertainty = uncertainty ) )

    def integrateTwoFunctionsWithUncertainty( self, f2, domainMin = None, domainMax = None, useCovariance = True, covariance = None, covarianceSuite = None, normalize = False ) :
        '''
        Computes the spectrum (i.e. weighted) integral of self with the spectrum (i.e. weight)
        specified by ``spectrum``  optionally from Emin to Emax.  If the covariance is
        provided, the uncertainty on the spectrum integral will be computed.

        :param f2: spectrum over which to average
        :type f2: XYs instance
        :param domainMin: Lower integration limit.
            If None (the default), then the lower limit of self's domain is used
        :type domainMin: PQU or None
        :param domainMax: Upper integration limit.
            If None (the default), then the upper limit of self's domain is used
        :type domainMax: PQU or None
        :param useCovariance: use this to override covarance usage
        :type useCovariance: bool
        :param covariance: covariance to use when computing uncertainty on the spectral average.
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :param normalize: if set to True, normalize integral of the spectrum over the interval so we are doing a spectrum average
        :type normalize: bool
        :rtype: PQU


        :How does it work?:
            Given a  weighting spectrum :math:`\phi(E)`, we intend to average the cross section in self as follows:

            .. math::
                < \sigma > = \int_{E_{min}}^{E_{max}} dE \phi(E) \sigma(E)

            To compute the uncertainty, we resort to the basis function expansion of the covariance
            as follows:

            .. math::
                \Delta^2\sigma(E,E') = \sum_{ij} \Delta^2\sigma_{ij} B_i(E) B_j(E')

            In the ENDF format, all covariances are assumed to be grouped in energy so the
            basis functions are simple window functions.   With this,

            .. math::
                \delta <\sigma> = \sqrt{ \int dE dE' \phi(E)\phi(E') \Delta^2\sigma(E,E') }
                                = \sqrt{ \sum_{ij} \Delta^2\sigma_{ij} <B_i> <B_j> }

        '''
        # Check that the inputs are of the correct type
        if not isinstance( f2, XYsModule.XYs ): raise TypeError( "spectrum must be an XYs instance")

        # Convert the cross section to pointwise
        ptwise = self.hasLinearForm()
        if ptwise is None: ptwise = self.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)

        # Check domains
        if domainMin is None: domainMin = PQU.PQU( max( ptwise.domainMin(), f2.domainMin() ), ptwise.domainUnit() )
        elif not isinstance( domainMin, PQU.PQU ): raise TypeError( "domainMin must be an PQU.PQU instance")
        if domainMax is None: domainMax = PQU.PQU( min( ptwise.domainMax(), f2.domainMax() ), ptwise.domainUnit() )
        elif not isinstance( domainMax, PQU.PQU ): raise TypeError( "domainMax must be an PQU.PQU instance")

        # Convolve spectrum and self to get mean value
        meanValue = ptwise.integrateTwoFunctions(f2,domainMin,domainMax)

        # Normalize the mean if we're averaging
        if normalize:
            norm = f2.integrate(domainMin,domainMax)
            meanValue = meanValue/norm

        # We might be done
        if not useCovariance: return meanValue # already a PQU

        # Get that the covariance goes with the data.
        import fudge.gnd.covariances.base as covModule
        if covariance is None:
            if hasattr(self.evaluated, 'uncertainties') and self.evaluated.uncertainties:
                if hasattr(self.evaluated.uncertainties[0].data,'link'):
                    covariance = self.evaluated.uncertainties[0].data.link
                elif covarianceSuite is not None:
                    covariance = self.evaluated.uncertainties[0].data.follow(startNode=covarianceSuite)
                else: return meanValue
            else: return meanValue # already a PQU
        else:
            if covariance is not( isinstance( covariance, covModule.covarianceMatrix ) ):
                raise TypeError( 'covariance must be of type covarianceMatrix, got %s'%str(type(covariance)))
        try:
            theCovariance = covariance.forms['eval'].toCovarianceMatrix() # may be redundant, but at least we'll get a usable data type
        except:
            print "WARNING: could not get covariance in integrateTwoFunctionsWithUncertainty for form %s" % str(covariance.forms['eval'].__class__)
            return meanValue

        # Compute weighting vector from spectrum and possibly the cross section (if covariance is relative)
        grid = theCovariance.matrix.axes[-1].grid
        gridUnit = theCovariance.matrix.axes[-1].unit
        covGroupBdries = list(grid)
        if theCovariance.type == 'absolute':
            phi = f2.group( covGroupBdries, norm=None )
        elif theCovariance.type == 'relative':
            phi = f2.group( covGroupBdries, ptwise, norm=None )
        else: raise ValueError( "Unknown covariance type: %s"%str(type(covariance)))

        # Compute the uncertainty itself
        import numpy
        phi = numpy.matrix( phi )
        theArray = theCovariance.matrix.array.constructArray()
        coco = ( phi * theArray * phi.T )[0,0]
        if coco < 0.0: print "WARNING: covariance of spectrum integral is %s < 0.0"%str(coco)
        uncertainty = PQU.PQU( math.sqrt( max( coco, 0.0 ) ), unit=meanValue.unit )
        if normalize: uncertainty = uncertainty/norm

        # it worked, return result
        return PQU.PQU( meanValue.value, unit=meanValue.unit, uncertainty=uncertainty.getValueAs(meanValue.unit) )

def parseXMLNode( crossSectionElement, xPath, linkData ):
    """
    Reads an xml <crossSection> element into fudge, including all cross section forms (pointwise, piecewise, etc.)
    contained inside the crossSection.
    """

    xPath.append( crossSectionElement.tag )
    crossSectionComponent = component( )

    for form in crossSectionElement :
        formClass = { 
                pointwise.moniker :   pointwise,
                piecewise.moniker :   piecewise,
                reference.moniker :   reference,
                resonancesWithBackground.moniker :   resonancesWithBackground,
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown cross section form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath = xPath, linkData = linkData )
        crossSectionComponent.add( newForm )

    xPath.pop( )
    return( crossSectionComponent )
