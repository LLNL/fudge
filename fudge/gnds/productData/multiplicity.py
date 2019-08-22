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

try :
    from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
    floatToShortestString = pointwiseXY_CModule.floatToShortestString
except :
    from pqu import PQU as PQUModule
    floatToShortestString = PQUModule.floatToShortestString

from PoPs import IDs as IDsPoPsModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import constant as constantModule
from xData import XYs as XYsModule
from xData import series1d as series1dModule
from xData import regions as regionsModule
from xData import gridded as griddedModule
from xData import link as linkModule

from fudge.core.utilities import brb
from fudge.gnds import tokens as tokensModule
from fudge.gnds import abstractClasses as abstractClassesModule

from fudge.gnds.reactions import base as reactionBaseModule

__metaclass__ = type

energyDependentToken = 'energyDependent'

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( 'multiplicity', 0, '' )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class baseMultiplicityForm( abstractClassesModule.form ) :

    pass

class unspecified( baseMultiplicityForm ) :

    moniker = 'unspecified'

    def __init__( self, label ) :

        baseMultiplicityForm.__init__( self )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def convertUnits( self, unitMap ) :

        pass

    def evaluate( self, E ) :

        return( 0. )

    def getXMLAttribute( self ) :

        return( self.moniker )

    @property
    def label( self ) :

        return( self.__label )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        raise Exception( 'Linear-linear XYs1d data for multiplicity form %s does not make sense' % self.moniker )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s/>' % ( indent, self.moniker, attributeStr ) ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        val = cls( element.get('label') )
        xPath.pop()
        return val

class constant1d( baseMultiplicityForm, constantModule.constant1d ) :

    def __init__( self, multiplicity, domainMin, domainMax, axes, label = None ) :

        baseMultiplicityForm.__init__( self )
        constantModule.constant1d.__init__( self, multiplicity, domainMin, domainMax, axes = axes, label = label )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def getXMLAttribute( self ) :

        integer = int( self.constant )
        if( self.constant == integer ) : return( integer )
        return( floatToShortestString( self.constant ) )

    def isConstant( self ) :

        return( True )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """This method returns the multiplicity as linear-linear XYs1d data which spans self's domain."""

        return( XYs1d( data = [ [ self.domainMin, self.constant ], [ self.domainMax, self.constant ] ], axes = self.axes ) )

class XYs1d( baseMultiplicityForm, XYsModule.XYs1d ) :

    def __init__( self, **kwargs ) :

        baseMultiplicityForm.__init__( self )
        XYsModule.XYs1d.__init__( self, **kwargs )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.domainMin, self.domainMax )

    def evaluate( self, E ) :

        if( E < self[0][0] ) : E = self[0][0]
        if( E > self[-1][0] ) : E = self[-1][0]
        return( XYsModule.XYs1d.evaluate( self, E ) )

    def getXMLAttribute( self ) :

        return( energyDependentToken )

    def isConstant( self ) :

        return self.rangeMin == self.rangeMax

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        if( tempInfo['verbosity'] > 2 ) : print('%sMulti-grouping XYs1d mutiplicity' % indent)

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class regions1d( baseMultiplicityForm, regionsModule.regions1d ) :

    def __init__( self, **kwargs ) :

        baseMultiplicityForm.__init__( self )
        regionsModule.regions1d.__init__( self, **kwargs )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ self[0][0][0], self[-1][-1][0] ] )

    def evaluate( self, E ) :

        if( E < self[0][0][0] ) : E = self[0][0][0]
        for region in self :
            if( E < region[-1][0] ) : break
        return( region.evaluate( E ) )

    def getXMLAttribute( self ) :

        return( energyDependentToken )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps =  1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """See regionsXYs.toPointwise_withLinearXYs on the use of lowerEps, upperEps."""

        kwargs['removeOverAdjustedPoints'] = True
        xys = regionsModule.regions1d.toPointwise_withLinearXYs( self, **kwargs )
        return( XYs1d( data = xys, axes = xys.axes ) )

class reference( linkModule.link, baseMultiplicityForm ) :

    moniker = tokensModule.referenceFormToken

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.link.__init__( self, link = link, root = root, path = path, label = label, relative = relative )
        baseMultiplicityForm.__init__( self )

    @property
    def reference( self ):
        if self.link is None:
            raise Exception("Unresolved link!")
        return self.link

    @property
    def domainMin( self ) :

        return( self.reference.domainMin )

    @property
    def domainMax( self ) :

        return( self.reference.domainMax )

    @property
    def domainUnit( self ) :

        return( self.reference.domainUnit )

    def convertUnits( self, unitMap ) :

        pass

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.reference.getEnergyLimits( EMin, EMax ) )

    def evaluate( self, E ) :

        return( self.reference.evaluate( E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.reference.toPointwise_withLinearXYs( **kwargs ) )

    def getXMLAttribute( self ) :

        return( self.reference.getXMLAttribute( ) )

class polynomial( baseMultiplicityForm, series1dModule.polynomial1d ) :

    def __init__( self, **kwargs ) :

        baseMultiplicityForm.__init__( self )
        series1dModule.polynomial1d.__init__( self, **kwargs )

    @property
    def domainUnit( self ) :

        return( self.findClassInAncestry( reactionBaseModule.base_reaction ).domainUnit )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def getXMLAttribute( self ) :

        return( energyDependentToken )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        xys = series1dModule.polynomial1d.toPointwise_withLinearXYs( self, accuracy = 1e-6 )
        return( XYs1d( data = xys, axes = self.axes ) )

class partialProduction( baseMultiplicityForm ) :

    moniker = tokensModule.partialProductionFormToken

    def __init__( self, label, data = 1 ) :

        baseMultiplicityForm.__init__( self )
        self.data = data

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def evaluate( self, E ) :

        return( self.data )

    def getXMLAttribute( self ) :

        return( tokensModule.partialProductionFormToken )

    @property
    def label( self ) :

        return( self.__label )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        brb.objectoutline( self.data )
        raise 'Hell'

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.data ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        label = element.get( 'label' )
        mult = float( element.get( 'value' ) )
        if( mult.is_integer( ) ) : mult = int( mult )
        multiplicity = partialProduction( label, mult )
        xPath.pop( )
        return( multiplicity )

class gridded1d( baseMultiplicityForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        baseMultiplicityForm.__init__( self )
        griddedModule.gridded1d.__init__( self, **kwargs )

    def getXMLAttribute( self ) :

        return( energyDependentToken )

class branching1d( baseMultiplicityForm ) :

    moniker = 'branching1d'

    def __init__( self, label ) :

        baseMultiplicityForm.__init__( self )

        self.__label = label

    @property
    def label( self ) :

        return( self.__label )

    def convertUnits( self, unitMap ) :

        pass

    def processMultiGroup( self, style, tempInfo, indent ) :

        def countGammas( product, probability ) :

            count = 0.0
            for decayMode in product.decayData.decayModes :
                branchingRatio = decayMode.probability[0].value * probability
                count += branchingRatio

                decayPath = decayMode.decayPath[0]
                residual = decayPath.products[0].pid
                if( residual == IDsPoPsModule.photon ) : residual = decayPath.products[1].pid
                count += countGammas( PoPs[residual], branchingRatio )

            return( count )

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']
        crossSection = tempInfo['crossSection']
        domainMin = crossSection.domainMin
        domainMax = crossSection.domainMax
        
        PoPs = tempInfo['reactionSuite'].PoPs

        distribution = self.ancestor.ancestor.distribution[0]
        nuclide = PoPs[distribution.pids.initial]

        multipliticy = countGammas( nuclide, 1.0 )
        multipliticy = XYs1d( data = [ [ domainMin, multipliticy ], [ domainMax, multipliticy ] ], axes = defaultAxes( tempInfo['incidentEnergyUnit'] ) )

        return( multipliticy.processMultiGroup( style, tempInfo, indent ) )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s/>' % ( indent, self.moniker, attributeStr ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        multiplicity = branching1d( element.get( 'label' ) )

        xPath.pop( )
        return( multiplicity )

def parseXMLNode( multElement, xPath, linkData ) :
    """Translate an xml <multiplicity> element into fudge."""

    xPath.append( multElement.tag )

    multiplicityComponent = component( )

    formClasses = {}
    for formClass in [ unspecified, constant1d, XYs1d, regions1d, reference, polynomial, partialProduction, gridded1d, branching1d ] :
        formClasses[formClass.moniker] = formClasses
    for form in multElement :
        formClass = formClasses.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown multiplicity form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath, linkData )
        multiplicityComponent.add( newForm )

    xPath.pop( )

    return( multiplicityComponent )

#
# multiplicity component
#
class component( abstractClassesModule.component ) :

    moniker = 'multiplicity'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( unspecified, constant1d, XYs1d, regions1d, reference, polynomial,
                partialProduction, gridded1d, branching1d ) )

    def isConstant( self ) :

        return( isinstance( self.evaluated, constant1d ) )

    def getConstant( self ) :

        if( self.isConstant()  ) : return( self.evaluated.evaluate( 0 ) )
        raise Exception( 'multiplicity type = "%s" is not single valued' % self.evaluated.moniker )

    @property
    def domainMin( self ) :

        return( self.evaluated.domainMin )

    @property
    def domainMax( self ) :

        return( self.evaluated.domainMax )

    @property
    def domainUnit( self ) :

        return( self.evaluated.domainUnit )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.evaluated.getEnergyLimits( EMin, EMax ) )
        
    def evaluate( self, E ) :

        return( self.evaluated.evaluate( E ) )

    def check( self, info ):

        from fudge.gnds import warning
        warnings = []

        if self.isConstant():
            if self.getConstant() < 1:
                warnings.append( warning.negativeMultiplicity( self.getConstant(), self ) )
        else:   # energy-dependent mult.
            for form in self :
                domain = None
                if hasattr(form, 'domain') :
                    domain = form.domainMin, form.domainMax
                if( ( domain is not None ) and ( domain != info['crossSectionDomain'] ) ) :
                    # photon multiplicities can start above cross section domainMin, but domainMax should agree:
                    if( self.ancestor.id == IDsPoPsModule.photon ) :
                        startRatio = domain[0] / info['crossSectionDomain'][0]
                        endRatio = domain[1] / info['crossSectionDomain'][1]
                        if (startRatio < 1 - standardsModule.floats.epsilon or endRatio < 1 - standardsModule.floats.epsilon
                            or endRatio > 1 + standardsModule.floats.epsilon):
                            warnings.append(warning.domain_mismatch(
                                *(domain + info['crossSectionDomain']), obj=self))
                    # For all other products, both domainMin and domainMax should agree within epsilon
                    else:
                        for idx, (e1, e2) in enumerate(zip(domain, info['crossSectionDomain'])):
                            ratio = e1 / e2
                            if (ratio < 1 - standardsModule.floats.epsilon or ratio > 1 + standardsModule.floats.epsilon):
                                if ( idx==0 and domain[0] > info['crossSectionDomain'][0] and form.evaluate( domain[0] ) == 0 ):
                                    continue    # multiplicity starting after xsc is ok if first multiplicity point == 0
                                warnings.append(warning.domain_mismatch( *(domain + info['crossSectionDomain']), obj=self))
                                break

                if hasattr(form, 'rangeMin') and form.rangeMin < 0:
                    warnings.append( warning.negativeMultiplicity( form.rangeMin, obj=form ) )

        return warnings

    def getXMLAttribute( self ) :

        return( self.evaluated.getXMLAttribute( ) )
