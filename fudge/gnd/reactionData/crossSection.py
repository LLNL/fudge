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

import math

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import axes as axesModule
from xData import link as linkModule
from xData import XYs as XYsModule
from xData import regions as regionsModule
from xData import gridded as griddedModule

from fudge.processing import group as groupModule

from fudge.gnd import abstractClasses as abstractClassesModule
from fudge.gnd import tokens as tokensModule
from fudge.gnd import styles as stylesModule

from . import base as baseModule

__metaclass__ = type

lowerEps = 1e-8
upperEps = 1e-8

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( baseModule.crossSectionToken, 0, 'b' )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

#
# crossSection forms.
#
class baseCrossSectionForm( abstractClassesModule.form ) :

    pass

class XYs1d( baseCrossSectionForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYsModule.XYs1d.__init__( self, **kwargs )

    def changeInterpolation( self, interpolation, accuracy, lowerEps = 0, upperEps = 0, cls = None ) :

        def func( x, parameters ) :

            self, T, subParameters = parameters
            x1, y1, x2, y2, f1, f2 = subParameters
            if( ( x is None ) or ( x < x1 ) or ( x > x2 ) ) :
                index = self.lowerIndexBoundingX( x )
                x1, y1 = self[index]
                if( ( index + 1 ) == len( self ) ) : return( y1 )
                x2, y2 = self[index+1]
                f1 = 1 / math.sqrt( x1 - T )
                f2 = 1 / math.sqrt( x2 - T )
                parameters[2] = [ x1, y1, x2, y2, f1, f2 ]

            if( x == x1 ) : return( y1 )
            if( x == x2 ) : return( y2 )
            f = 1 / math.sqrt( x - T )
            alpha = ( f - f1 ) / ( f2 - f1 )
            return( math.pow( x1 * y1, 1 - alpha ) * math.pow( x2 * y2, alpha ) / x )

        if( cls is None ) : cls =XYs1d 
        if( self.interpolation == standardsModule.interpolation.chargedParticleToken ) :
            if( interpolation != standardsModule.interpolation.linlinToken ) : raise TypeError( 'Only "%s" interpolation for conversion from %s' %
                    ( standardsModule.interpolation.linlinToken, self.interpolation ) )
            temp = self.cloneToInterpolation( interpolation )
            parameters = [ temp, 0, [ None, None, None, None, None, None ] ]
            return( temp.applyFunction( func, parameters, accuracy = accuracy, biSectionMax = -1, checkForRoots = False ) )

        return( XYsModule.XYs1d.changeInterpolation( self, interpolation, accuracy, lowerEps = lowerEps,
                upperEps = upperEps, cls = cls ) )

    def heat( self, currentTemperature, newTemperature, massRatio, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, 
                heatAllPoints = False, doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True,
                setThresholdToZero = False ) :
        """
        Returns a linear version of the cross section heated to 'newTemperature'. If the current temperature of the 
        cross section, given by 'currentTemperature', is greater than 'newTemperature' a raise is executed.
        If lowerlimit is None, it is set to 'oneOverV' except when the reaction is determined to be a threshold reaction, 
        then it is set to 'threshold'. Any cross section with domainMin greater than 2.5e-4 eV is determined to be a
        threshold reaction. If upperlimit is None it is set to 'constant'.
        If heatBelowThreshold is False, then EMin is set to the larger of EMin and self's domainMin.

        The unit of 'currentTemperature' and 'newTemperature' must be the same as the self's domain unit.

        For more information on EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin and heatAllEDomain
        see the module crossSectionAdjustForHeatedTarget.
        """

        from crossSectionAdjustForHeatedTarget import heat

        dT = newTemperature - currentTemperature
        if( abs( dT ) <= 1e-2 * newTemperature ) : dT = 0.
        if( dT < 0 ) : raise Exception( 'Current temperature "%s" (in energy units) higher than desired temperature "%s"' %
                ( currentTemperature, newTemperature ) ) 

        heated = unheated = self
        if( not( unheated.isInterpolationLinear( ) ) ) :
            unheated = self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = upperEps )
        if( ( len( unheated ) > 0 ) and ( dT > 0. ) ) :
            if( not( heatBelowThreshold ) ) : EMin = max( unheated.domainMin, EMin )
            if( lowerlimit is None ) :
                lowerlimit = "oneOverV"
# BRB6 Hardwired.
                domainMin = PQUModule.PQU( 1e-5, 'eV' ).getValueAs( self.domainUnit )
                if( domainMin < 0.04 * unheated.domainMin ) : lowerlimit = "threshold"
            if( upperlimit is None ) : upperlimit = 'constant'

            heated = heat.crossSectionAdjustForHeatedTarget( massRatio, dT, EMin, unheated, lowerlimit = lowerlimit,
                upperlimit = upperlimit, interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                heatAllEDomain = heatAllEDomain )
        if( ( unheated[0][1] == 0.0 ) and setThresholdToZero ) : heated[0][1] = 0.0
        return( XYs1d( data = heated, axes = unheated.axes ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print '%sMulti-grouping XYs1d cross section' % indent

        crossSectionGrouped = miscellaneousModule.groupOneFunctionAndFlux( style, tempInfo, self )
        return( groupModule.toMultiGroup1d( gridded1d, style, tempInfo, self.axes, crossSectionGrouped ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        ptw = XYsModule.XYs1d.toPointwise_withLinearXYs( self, cls = XYs1d, hh = True, **kwargs )
        ptw.label = self.label
        return( ptw )

    def toXMLList( self, indent = "", **kwargs ) :

        def dataToString( self, dummy, indent = '', **kwargs ) :

            maxXFormat = "%.17e"

            valuesPerLine = int( kwargs.get( 'valuesPerLine', 100 ) )
            if( ( valuesPerLine % 2 ) == 1 ) : valuesPerLine += 1
            xFormat = kwargs.get( "xFormat", "%14.8e" )

            XMLList = []
            Line = []
            i1 = 1
            for value in self :
                if( ( i1 % 2 ) != 0 ) :
                    xString = xFormat % value
                    if( i1 > 1 ) :
                        if( xString == priorXString ) : 
                            if( xFormat == maxXFormat ) :
                                raise ValueError( 'x-value strings identical: "%s" and "%s"' % ( priorXString, xString ) )
                            kwargs['xFormat'] = maxXFormat
                            return( dataToString( self, dummy, indent = indent, **kwargs ) )
                    Line.append( xString )
                    priorXString = xString
                else :
                    Line.append( "%12.6e" % value )
                if( ( i1 % valuesPerLine ) == 0 ) :
                    XMLList.append( indent + ' '.join( Line ) )
                    Line = []
                i1 += 1
            if( len( Line ) > 0 ) : XMLList.append( indent + ' '.join( Line ) )
            return( XMLList )

        kwargs['dataToString'] = dataToString
        kwargs['dataToStringParent'] = None
        return( XYsModule.XYs1d.toXMLList( self, indent, **kwargs ) )

class regions1d( baseCrossSectionForm, regionsModule.regions1d ) :
    """
    This class stores a cross section in two or more regions, which may have different interpolations.
    Each region must contain at least two points. Each pair of adjacent regions must overlap at exactly one point.
    """

    def __init__( self, **kwargs ) :

        regionsModule.regions1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        return( linear.processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        xys = regionsModule.regions1d.toPointwise_withLinearXYs( self, **kwargs )
        return( XYs1d( data = xys, axes = xys.axes ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class gridded1d( baseCrossSectionForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )

class resonanceLink( linkModule.link ) :

    moniker = 'resonanceRegion'

class resonancesWithBackground( baseCrossSectionForm ) :
    """
    This class stores cross sections that include resonances along with a background contribution.
    Contains a link to the resonances, and the 'background' which could be XYs1d or regions1d data.
    The full XYs1d cross section can be obtained by first reconstructing resonances
    and then adding the background contribution (users should use the reactionSuite.reconstructResonances method).
    """

    moniker = tokensModule.resonancesWithBackgroundFormToken
    ancestryMembers = ( 'tabulatedData', )

    def __init__( self, label, tabulatedData, resonanceLink ) :

        baseCrossSectionForm.__init__( self )

        self.link = resonanceLink

        self.tabulatedData = tabulatedData
        self.tabulatedData.attribute = None
        self.tabulatedData.label = None
        self.tabulatedData.setAncestor( self )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    @property
    def domainMin( self ) :

        return( self.tabulatedData.domainMin )

    @property
    def domainMax( self ) :

        return( self.tabulatedData.domainMax )

    @property
    def domainUnit( self ) :

        return( self.tabulatedData.domainUnit )

    def findEntity( self, entityName, attribute = None, value = None ):

        if entityName in ( XYs1d.moniker, regions1d.moniker ) :
            return self.tabulatedData
        return ancestryModule.ancestry.findEntity( self, entityName, attribute, value )

    @property
    def label( self ) :

        return( self.__label )
    
    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.tabulatedData.convertUnits( unitMap )

# BRB FIXME This is a kludge until we fix production/crossSection/reference
# For example, see n + U236 in ENDF/B-7.1.
# <productions>
#    <production label="0" outputChannel="U235" ENDF_MT="16">
#      <crossSection>
#        <reference label="eval" ...
    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.getAncestor( ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        _component = self.findClassInAncestry( component )
        reconstructed = _component.getStyleOfClass( stylesModule.crossSectionReconstructed )
        if( reconstructed is not None ) :
            return( reconstructed.toPointwise_withLinearXYs( **kwargs ) )     # Return a copy.
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
        xlink = resonanceLink.parseXMLNode( element.find('resonanceRegion'), xPath, linkData )
        form = element[1]
        formClass = {   XYs1d.moniker               : XYs1d,
                        regions1d.moniker           : regions1d,
                        weightedPointwise.moniker   : weightedPointwise,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown cross section form: %s" % form.tag)
        resWithBack = resonancesWithBackground( element.get( 'label' ), formClass.parseXMLNode(form,xPath,linkData), xlink )
        xPath.pop()
        return resWithBack
    
class reference( linkModule.link, baseCrossSectionForm ) :
    """This cross section form consists of a reference to another cross section."""

    moniker = tokensModule.referenceFormToken

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.link.__init__(self, link = link, root = root, path = path, label = label, relative = relative )
        baseCrossSectionForm.__init__( self )

    @property
    def crossSection( self ) :

        if self.link is None : raise Exception( "Unresolved link!" )
        return self.link

    @property
    def domainMin( self ) :

        return( self.crossSection.domainMin )

    @property
    def domainMax( self ) :

        return( self.crossSection.domainMax )

    @property
    def domainUnit( self ) :

        return( self.crossSection.domainUnit )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def getReference( self ) :

        return( self.crossSection )

    def setReference( self, crossSection ) :

        if( not( isinstance( crossSection, (component, None.__class__) ) ) ) :
            raise TypeError( 'crossSection argument must be a cross section component not type %s' % type( crossSection ) )
        self.crossSection = crossSection

    def processMultiGroup( self, style, tempInfo, indent ) :

        addToComponent = tempInfo.get( 'addToComponent', None )
        tempInfo['addToComponent'] = False
        multiGroup = self.crossSection.processMultiGroup( style, tempInfo, indent )
        tempInfo['addToComponent'] = addToComponent
        if( addToComponent is None ) : del tempInfo['addToComponent']
        return( multiGroup )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.crossSection.toPointwise_withLinearXYs( **kwargs ) )

    def check( self ) :

        from fudge.gnd import warning
        warnings = []
        if self.getRootAncestor() != self.getReference().getRootAncestor():
            warnings.append( warning.badCrossSectionReference() )
        return warnings

class CoulombElasticReference( reference ) :
    """Special type of link: points to dCrossSection_dOmega/CoulombElastic"""

    moniker = "CoulombElasticReference"

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

        baseCrossSectionForm.__init__( self )
        self.setReference( crossSection )
        self.weights = weights

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    @property
    def domainMin( self ) :

        return( self.weights.domainMin )

    @property
    def domainMax( self ) :

        return( self.weights.domainMax )

    @property
    def domainUnit( self ) :

        return( self.crossSection.domainUnit )

    def setReference( self, crossSection ):

        self.crossSection = crossSection

    @property
    def label( self ) :

        return( self.__label )
    
    def toPointwise_withLinearXYs( self, **kwargs ) :

        xys = self.crossSection.toPointwise_withLinearXYs( **kwargs ) * self.weights
        return( XYs1d( data = xys, axes = xys.axes ) )

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
        weights = XYsModule.XYs1d.parseXMLNode( WPelement[0], xPath, linkData )
        ref = weightedPointwise( None, weights )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

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

    def __init__( self ) :

        abstractClassesModule.component.__init__( self,
                ( XYs1d, regions1d, gridded1d, resonanceLink, resonancesWithBackground, reference, weightedPointwise ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainMin( self ) :

        return( self.evaluated.domainMin )

    @property
    def domainMax( self ) :

        return( self.evaluated.domainMax )

    @property
    def domainUnit( self ) :

        return( self.evaluated.domainUnit )

    def hasLinearForm( self ) :

        for form in self :
            if( isinstance( form, XYs1d ) ) :
                if( form.isInterpolationLinear( ) ) : return( form )
        return( None )

    def heat( self, style, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, 
            heatAllPoints = False, doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False,
            addToSuite = False ) : 
        """
        Returns the result of self.toPointwise_withLinearXYs( ).heat( ... ). See method XYs1d.heat for more information.
        If setThresholdToZero is True and self's cross section at the first point is 0., then the heated cross section's
        first value will also be 0.
        """

        styles = self.findAttributeInAncestry( 'styles' )
        currentstyle = styles[style.derivedFrom]
        if( not( isinstance( currentstyle, ( stylesModule.evaluated, stylesModule.heated ) ) ) ) :
            TypeError( 'Form to heat is not heatable form: its style moniker is "%s"' % style.label )
        currentTemperature = PQUModule.PQU( currentstyle.temperature.getValueAs( 'K' ), 'K' ) * PQUModule.PQU( '1 k' )
        currentTemperature = currentTemperature.getValueAs( self.domainUnit )

        newTemperature = PQUModule.PQU( style.temperature.getValueAs( 'K' ), 'K' ) * PQUModule.PQU( '1 k' )
        newTemperature = newTemperature.getValueAs( self.domainUnit )

        reactionSuite = self.getRootAncestor( )
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.getMass( 'amu' )
        if( projectileMass == 0 ) : raise ValueError( 'Heating with gamma as projectile not supported.' )

        target = reactionSuite.PoPs[reactionSuite.target]
        massRatio = target.getMass( 'amu' ) / projectileMass

        linear = self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        heated = linear.heat( currentTemperature, newTemperature, massRatio, EMin, lowerlimit, upperlimit, interpolationAccuracy, 
                heatAllPoints, doNotThin, heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero )
        heated.label = style.label
        if( addToSuite ) : self.add( heated )
        return( heated )

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
        lower = PQUModule.PQU( self.domainMin, self.domainUnit )
        upper = PQUModule.PQU( self.domainMax, self.domainUnit )
# BRB6 hardwired
        thresh = PQUModule.PQU( 0, 'eV' )
        if 'Q' in info:
            thresh = -info['Q'] * info['kinematicFactor']
            if not isinstance(info['Q'], PQUModule.PQU): thresh = PQUModule.PQU(thresh, 'eV')

            if not info['CoulombChannel']:
                # if Q is positive, cross section must start at 1e-5 eV
                # otherwise, cross section threshold must agree with Q-value:
                if thresh.value>0:
                    if abs(thresh-lower) > PQUModule.PQU( info['dThreshold'] ):
                        warnings.append( warning.threshold_mismatch( lower, thresh, self ) )
                elif lower != PQUModule.PQU(1e-5,'eV'):
                    # ignore 2nd,3rd,4th-chance fission (they have Q>0 but are still threshold reactions):
                    from fudge.gnd import channels
                    parent = self.getAncestor( )
                    if (not hasattr( parent, 'outputChannel' ) or ( parent.outputChannel.getFissionGenre()
                            in (None, channels.fissionGenreTotal, channels.fissionGenreFirstChance) ) ):
                        warnings.append( warning.threshold_mismatch( lower, PQUModule.PQU(1e-5, 'eV'), self ) )
            else:
                # charged-particle reaction generally doesn't 'turn on' right at threshold due to Coulomb barrier.
                # In this case, just ensure the cross section stays zero until at or above threshold:
                if (lower < thresh):
                    warnings.append( warning.Coulomb_threshold_mismatch( lower, thresh, self ) )

        # cross section must extend up to limit (usually 20 MeV):
        if upper < PQUModule.PQU( info['crossSectionEnergyMax'] ):
            warnings.append( warning.gapInCrossSection( upper, PQUModule.PQU( info['crossSectionEnergyMax'] ),
                self ) )

        evaluatedCrossSection = self.evaluated
        if( isinstance( evaluatedCrossSection, resonancesWithBackground ) ) :
            # ensure no gaps between resonance and background:
            resonances = info['reactionSuite'].resonances
            if resonances.resolved is not None:
                resDomainMax = resonances.resolved.domainMax
            if resonances.unresolved is not None and resonances.unresolved.reconstructCrossSection:
                resDomainMax = resonances.unresolved.domainMax
            bckDomainMin = evaluatedCrossSection.tabulatedData.domainMin
            if bckDomainMin > resDomainMax:
                warnings.append( warning.gapInCrossSection(resDomainMax,bckDomainMin, self ) )
            linearCrossSection = self[info['reconstructedStyleName']]

        elif( isinstance( evaluatedCrossSection, reference )
            and isinstance( evaluatedCrossSection.link, resonancesWithBackground ) ) :
                linearCrossSection = evaluatedCrossSection.link.ancestor[info['reconstructedStyleName']]

        elif isinstance( evaluatedCrossSection, CoulombElasticReference ):
            return warnings     # already checked as part of dCrossSection_dOmega

        else:
            # can it be converted to XYs1d linear?
            try:
                linearCrossSection = evaluatedCrossSection.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
            except Exception as e:
                warnings.append( warning.ExceptionRaised( e, self ) )
                if info['failOnException']: raise
                return warnings

        # test for negative values, and for non-zero cross section at threshold
        if( linearCrossSection.rangeMin < 0 ) :
            for i,(en,xsc) in enumerate(linearCrossSection):
                if xsc < 0:
                    warnings.append( warning.negativeCrossSection( PQUModule.PQU(en, linearCrossSection.axes[-1].unit), i, self ) )

        if thresh.value>0:
            if linearCrossSection[0][1] != 0:
                warnings.append( warning.nonZero_crossSection_at_threshold( linearCrossSection[0][1], self ) )

        return warnings

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

        if( not isinstance( E, PQUModule.PQU ) ) : raise TypeError( "E must be an PQUModule.PQU instance")

        ptwise = self.hasLinearForm()
        if ptwise is None: ptwise = self.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )

        EinDomainUnit = E.getValueAs( ptwise.domainUnit )
        if EinDomainUnit < ptwise.domainMin or EinDomainUnit > ptwise.domainMax: meanValue = 0.0
        else: meanValue = ptwise.evaluate( EinDomainUnit )

        # We might be done
        if not useCovariance: return PQUModule.PQU( meanValue, unit=ptwise.rangeUnit )

        # Get that the covariance goes with the data.
        covariance = self.getMatchingCovariance(covariance,covarianceSuite)

        if covariance is None: return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit )
        else:
            try:
                return( PQUModule.PQU(
                    meanValue,
                    unit = ptwise.rangeUnit,
                    uncertainty = covariance.getUncertaintyVector( self.evaluated, relative = False ).evaluate(EinDomainUnit) ) )
            except IndexError as err: # FIXME: a kludge until cov routines working again, try it on 27Al(n,tot)
                print "WARNING: Could not extract uncertianty, got error %s"%str(err)
                return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit )

#        if covariance is None:
#            if hasattr(self.evaluated, 'uncertainties') and self.evaluated.uncertainties:
#                if hasattr(self.evaluated.uncertainties[0].data,'link'):
#                    covariance = self.evaluated.uncertainties[0].data.link['eval']
#                elif covarianceSuite is not None:
#                    covariance = self.evaluated.uncertainties[0].data.follow(startNode=covarianceSuite)['eval']
#                else: return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#            else: return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#        else:
#            if covariance is not( isinstance( covariance, covModule.covarianceMatrix ) ):
#                raise TypeError( 'covariance must be of type covarianceMatrix, got %s'%str(type(covariance)))

#        # Compute uncertainty on convolution given mean value and covariance.
#        try:
#            uncertainty = covariance.getUncertaintyVector( self.evaluated, relative = False ).evaluate( EinDomainUnit )
#        except:
#            print "WARNING: could not get uncertainty in evaluateWithUncertainty for %s" % str(covariance.__class__)
#            return PQUModule.PQU( meanValue, unit = ptwise.rangeUnit() )
#        if uncertainty is None:

#        # We haven't changed units on uncertainty or the meanValue, so we can reconstruct the physical quantity easily.
#        return( PQUModule.PQU( meanValue, unit = ptwise.rangeUnit(), uncertainty = uncertainty.evaluate(EinDomainUnit) ) )

    def getMatchingCovariance(self, covariance = None, covarianceSuite = None):
        '''
        Retrieve the uncertainty that goes with this cross section, if we can find it
        :return:
        '''
        import fudge.gnd.covariances.base as covModule

        # Get that the covariance goes with the data.
        if covariance is None:
            if hasattr(self.evaluated, 'uncertainties') and self.evaluated.uncertainties:
                if hasattr(self.evaluated.uncertainties[0].data,'link') and self.evaluated.uncertainties[0].data.link is not None:
                    covariance = self.evaluated.uncertainties[0].data.link['eval']
                elif covarianceSuite is not None:
                    covariance = self.evaluated.uncertainties[0].data.follow(startNode=covarianceSuite)['eval']
        else:
            if covariance is not( isinstance( covariance, covModule.covarianceMatrix ) ):
                raise TypeError( 'covariance must be of type covarianceMatrix, got %s'%str(type(covariance)))

        return covariance

    def integrateTwoFunctionsWithUncertainty( self, f2, domainMin = None, domainMax = None, useCovariance = True, covariance = None, covarianceSuite = None, normalize = False ) :
        '''
        Computes the spectrum (i.e. weighted) integral of self with the spectrum (i.e. weight)
        specified by ``spectrum``  optionally from Emin to Emax.  If the covariance is
        provided, the uncertainty on the spectrum integral will be computed.

        :param f2: spectrum over which to average
        :type f2: XYs1d instance
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
        if not isinstance( f2, XYsModule.XYs1d ): raise TypeError( "spectrum must be an XYs1d instance")

        # Convert the cross section toXYs1d 
        ptwise = self.hasLinearForm()
        if ptwise is None: ptwise = self.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )

        # Check domains
        if( domainMin is None ) :
            domainMin = PQUModule.PQU( max( ptwise.domainMin, f2.domainMin ), ptwise.domainUnit )
        elif( not isinstance( domainMin, PQUModule.PQU ) ) :
            raise TypeError( "domainMin must be an PQUModule.PQU instance")

        if( domainMax is None ) :
            domainMax = PQUModule.PQU( min( ptwise.domainMax, f2.domainMax ), ptwise.domainUnit )
        elif( not isinstance( domainMax, PQUModule.PQU ) ) :
            raise TypeError( "domainMax must be an PQUModule.PQU instance")

        # Convolve spectrum and self to get mean value
        meanValue = ptwise.integrateTwoFunctions(f2,domainMin,domainMax)

        # Normalize the mean if we're averaging
        if normalize and meanValue.value != 0.0:
            norm = f2.integrate(domainMin,domainMax)
            if norm.value == 0.0: raise ValueError('zero norm (%s) while integrating function with %s '%(str(norm),str(self.ancestor.ancestor)))
            meanValue = meanValue/norm

        # We might be done
        if not useCovariance: return meanValue # already a PQU

        # Get that the covariance goes with the data.
        covariance = self.getMatchingCovariance(covariance,covarianceSuite)
        if covariance is None: return meanValue

        try:
            theCovariance = covariance.toCovarianceMatrix() # may be redundant, but at least we'll get a usable data type
        except Exception, err:
            print "WARNING: could not get covariance in integrateTwoFunctionsWithUncertainty for form %s, got error message '%s'" % (str(covariance.__class__),err)
            return meanValue

        # Compute weighting vector from spectrum and possibly the cross section (if covariance is relative)
        grid = theCovariance.matrix.axes[-1].values.values
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
        try:
            coco = ( phi * theArray * phi.T )[0,0]
            if coco < 0.0: print "WARNING: covariance of spectrum integral is %s < 0.0"%str(coco)
            uncertainty = PQUModule.PQU( math.sqrt( max( coco, 0.0 ) ), unit=meanValue.unit )
            if normalize and meanValue.value != 0.0: uncertainty = uncertainty/norm

            # it worked, return result
            return PQUModule.PQU( meanValue.value, unit=meanValue.unit, uncertainty=uncertainty.getValueAs(meanValue.unit) )
        except Exception, err:
            print "WARNING: could not compute uncertainty in integrateTwoFunctionsWithUncertainty for form %s, got error message '%s'" % (str(covariance.__class__),err)
            return meanValue

def parseXMLNode( crossSectionElement, xPath, linkData ):
    """
    Reads an xml <crossSection> element into fudge, including all cross section forms (XYs1d, regions1d, etc.)
    contained inside the crossSection.
    """

    xPath.append( crossSectionElement.tag )
    crossSectionComponent = component( )

    for form in crossSectionElement :
        formClass = { 
                XYs1d.moniker                       : XYs1d,
                regions1d.moniker                   : regions1d,
                gridded1d.moniker                   : gridded1d,
                reference.moniker                   : reference,
                CoulombElasticReference.moniker     : CoulombElasticReference,
                resonancesWithBackground.moniker    : resonancesWithBackground,
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown cross section form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath = xPath, linkData = linkData )
        crossSectionComponent.add( newForm )

    xPath.pop( )
    return( crossSectionComponent )
