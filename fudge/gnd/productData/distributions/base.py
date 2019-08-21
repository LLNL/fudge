# <<BEGIN-copyright>>
# <<END-copyright>>

""" base classes for distributions """

import math

from fudge.core.ancestry import ancestry
from fudge.core.math import *   #math classes, fudgemath, fudge2dGrouping
from fudge.legacy.converting import endfFormats
from fudge.core.utilities import fudgeZA
import fudge
from fudge.processing import processingInfo
from fudge.processing.montecarlo import fudge2dEqualProbableBinning
from fudge.processing.deterministic import transferMatrices
import miscellaneous

__metaclass__ = type

#
# Standard genre can only be composited from one of the following standard forms.
#

angularTwoBodyGenre = 'angularTwoBody'
unknownGenre = 'unknown'
referenceGenre = 'reference'
NBodyGenre = 'NBody'

noneComponentToken = 'none'
unknownComponentToken = 'unknown'
angularComponentToken = 'angular'
energyComponentToken = 'energy'
energyAngularComponentToken = 'energyAngular'
angularEnergyComponentToken = 'angularEnergy'
LegendreComponentToken = 'Legendre'
LegendreEnergyConservationComponentToken = 'LegendreEnergyConservation'
referenceComponentToken = 'reference'       # 'reference', should be replaced with 'link':
CoulombElasticComponentToken = 'CoulombElastic'
LLNLAngularEnergyComponentToken = 'LLNLAngularEnergy'
LLNL_withAngularComponentToken = 'LLNLAngular_angularEnergy'
uncorrelatedComponentToken = 'uncorrelated'

noneFormToken = 'none'
constantFormToken = 'constant'  # angular.py overrides with isotropicFormToken
isotropicFormToken = 'isotropic'
linearFormToken = 'linear'
pointwiseFormToken = 'pointwise'
piecewiseFormToken = 'piecewise'
semiPiecewiseFormToken = 'semiPiecewise'
equalProbableBinsFormToken = 'equalProbableBins'
nuclearPlusCoulombInterferenceFormToken = 'nuclearPlusCoulombInterference'
mixedRangesFormToken = 'mixedRanges'
recoilFormToken = 'recoil'
groupedFormToken = 'grouped'
LegendrePointwiseFormToken = 'LegendrePointwise'
LegendrePiecewiseFormToken = 'LegendrePiecewise'
LegendrePiecewise_NuclearFormToken = 'LegendrePiecewise_Nuclear'
LegendrePiecewise_CoulombInterferenceRealFormToken = 'LegendrePiecewise_CoulombInterferenceReal'
LegendrePiecewise_CoulombInterferenceImaginaryFormToken = 'LegendrePiecewise_CoulombInterferenceImaginary'
LLNLLegendrePointwiseFormToken = 'LLNLLegendrePointwise'
NBodyPhaseSpaceFormToken = 'NBodyPhaseSpace'
generalEvaporationFormToken = 'generalEvaporation'
simpleMaxwellianFissionFormToken = 'simpleMaxwellianFission'
evaporationFormToken = 'evaporation'
WattFormToken = 'Watt'
MadlandNixFormToken = 'MadlandNix'
weightedFunctionalsFormToken = 'weightedFunctionals'
KalbachMannFormToken = 'KalbachMann'

class distribution( ancestry ) :

    def __init__( self, nativeData, moniker = 'distributions' ) :

        ancestry.__init__( self, moniker, None )
        self.nativeData = nativeData
        self.components = {}

    def __len__( self ) :
        "Returns the number of components representing this data."

        return( len( self.components ) )

    def __getitem__( self, name ) :
        "Returns the component named name."

        return( self.components[name] )

    def addComponent( self, component ) :

        self.components[component.moniker] = component
        component.setParent( self )

    def getComponentNames( self ) :
        """Returns a list of the name of each component currently in self."""

        return( self.components.keys( ) )

    def getNativeData( self ) :

        return( self.components[self.nativeData].getNativeData( ) )

    def getNativeDataToken( self ) :

        return( self.nativeData )

    def calculateDepositionData( self, processInfo, tempInfo, name = None ) :

        if( name is None ) : name = self.nativeData
        if( name == noneFormToken ) : return( [] )
        return( self.components[name].calculateDepositionData( processInfo, tempInfo ) )

    def process( self, processInfo, tempInfo, verbosityIndent, addToDistribution = True ) :

        components = self.components[self.nativeData].process( processInfo, tempInfo, verbosityIndent )
        if( addToDistribution ) :
            for component in components :
                if( component.moniker in self.components ) :
                    self.components[component.moniker].addForm( component[component.nativeData] )
                else :
                    self.addComponent( component )
        return( components )

    def setNativeData( self, nativeData ) :

        self.nativeData = nativeData

    def getGropuedChiAndDepEnergyAtLowestE( self, groupBoundaries ) :
        """This is a special routine to handle delayed fission neutrons for LLNL_SnCOut."""

        try :
            return( self.components[LegendreComponentToken].getGropuedChiAndDepEnergyAtLowestE( groupBoundaries ) )
        except :
            return( None )

    def toLLNLSnCOut( self, toOtherData, moreOtherData, name, addToFissionMatrix ) :

        self.components[LegendreComponentToken].toLLNLSnCOut( toOtherData, moreOtherData, name, addToFissionMatrix )
        self.components[LegendreEnergyConservationComponentToken].toLLNLSnCOut( toOtherData, moreOtherData, name )

    def addAttributesToGNDString( self ) :

        return( '' )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.components[self.nativeData].toPointwiseLinear( accuracy, lowerEps, upperEps ) )

    def toXMLList( self, indent = '' ) :

        xmlString = [ '%s<%s nativeData="%s"%s>' % ( indent, self.moniker, self.nativeData, self.addAttributesToGNDString( ) ) ]
        names = sorted( self.components.keys( ) )
        if( angularComponentToken in names ) :
            names.remove( angularComponentToken )
            names.insert( 0, angularComponentToken )
        for name in names : xmlString += self.components[name].toXMLList( indent = indent + '  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        if( self.nativeData == noneComponentToken ) :
            if( MT in [ 452, 455, 456 ] ) : return
            targetInfo['MF6LCTs'].append( None )
            particleID = targetInfo[ 'zapID' ]    # get ZAP for the outgoing particle
            if( particleID == 'gamma' ) :
                ZAP = 0
            else :
                Z, A, suffix, ZAP = fudgeZA.gndNameToZ_A_Suffix( particleID )
            nPoints, multiplicityList = targetInfo['multiplicity'].endfMultiplicityList( targetInfo )
            endfMFList[6][MT] += [ endfFormats.endfContLine( ZAP, targetInfo['particleMass'] / targetInfo['neutronMass'], 0, 0, 1, nPoints ) ]
            endfMFList[6][MT] += multiplicityList
            return
        if( hasattr( self.components[self.nativeData], 'toENDF6' ) ) :
            self.components[self.nativeData].toENDF6( MT, endfMFList, flags, targetInfo )
        elif self.nativeData == referenceComponentToken:
            pass
        else :
            print 'WARNING: Distribution, no toENDF6 for nativeData = %s' % self.nativeData

class component( ancestry ) :

    def __init__( self, moniker, nativeData ) :

        ancestry.__init__( self, moniker, None )
        self.forms = {}
        self.nativeData = nativeData

    def __len__( self ) :

        return( len( self.forms ) )

    def __getitem__( self, name ) :
        "Returns the formed named name."

        return( self.forms[name] )

    def addForm( self, form ) :

        self.forms[form.moniker] = form
        form.component = self.moniker
        form.setParent( self )

    def getFormNames( self ) :
        """Returns a list of the name of each form currently in self."""

        return( self.forms.keys( ) )

    def getNativeData( self ) :

        return( self.forms[self.nativeData] )

    def getNativeDataToken( self ) :

        return( self.nativeData )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        raise Exception( '"calculateDepositionData" not implemented for component = "%s"' % self.moniker )
        
    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( '"process" not implemented for component = "%s"' % self.moniker )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.forms[self.nativeData].toPointwiseLinear( accuracy, lowerEps, upperEps ) )

    def toXMLList( self, indent = '' ) :

        xmlString = [ '%s<%s nativeData="%s">' % ( indent, self.moniker, self.nativeData ) ]
        for name in self.forms : xmlString += self.forms[name].toXMLList( indent = indent + '  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class form( ancestry ) :

    def __init__( self, moniker ) :

        ancestry.__init__( self, moniker, None )

    def getName( self ) :

        return( self.moniker )

class distributionFormWithRegions( form ) :

    def __init__( self, name ) :

        form.__init__( self, name )
        self.data = []

    def __len__( self ) :

        return( len( self.data ) )

    def __getitem__( self, index ) :

        return( self.data[index] )

    def append( self, region ) :

        if( len( self.data ) != region.index ) : raise Exception( 'region index = %d is invalid: must be %d' % ( region.index, len( self.data ) ) )
        if( len( self.data ) > 0 ) :
            xMin1, xMax1 = self.data[-1].domain( )
            xMin2, xMax2 = region.domain( )
            if( xMax1 != xMin2 ) : raise Exception( 'domains = %e and %e not the same' % ( xMax1, xMin2 ) )
        self.data.append( region )

    def toXMLList( self, indent = ""  ) :

        xmlString = []
        for region in self : xmlString += region.toXMLList( indent )
        return( xmlString )

class equalProbableBinsFormBase( form ) :

    def __init__( self ) :

        form.__init__( self, equalProbableBinsFormToken )
        self.numberOfBins = len( data[0][1] ) - 1

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlString[-1] += '<equalProbableBins2d bins="%d">' % ( self.numberOfBins )
        xmlString += self.axes.toXMLList( indent = indent2 )
        for indexE, energyMu in enumerate( self.data ) :
            energy, mu = energyMu
            muString = endl1dmathmisc.list1dToXMLEqualProbableBins1dString( mu )
            xmlString.append( '%s<energy value="%s" index="%d">%s</energy>' % ( indent2, fudge.gnd.miscellaneous.floatToString( energy ), indexE, muString ) )
        xmlString[-1] += '</equalProbableBins2d></%s>' % self.moniker
        return( xmlString )

class referenceComponent( component ) :

    def __init__( self, referenceInstance ) :
        """Link to another part of the file using link.  'reference' should use xPath syntax."""

        component.__init__( self, referenceComponentToken, noneFormToken )
        self.setReference( referenceInstance )

    def setReference( self, referenceInstance ) :

        self.referenceInstance = referenceInstance

    def calculateDepositionData( self, processInfo, tempInfo ) :

        return( self.referenceInstance.calculateDepositionData( processInfo, tempInfo ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.referenceInstance.process( processInfo, tempInfo, verbosityIndent, addToDistribution = False ) )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.referenceInstance.toPointwiseLinear( accuracy, lowerEps, upperEps ) )

    def toXMLList( self, indent = "" ) :

        return( [ '%s<%s xlink:type="simple" xlink:href="%s"/>' % ( indent, self.moniker, self.referenceInstance.toXLink( ) ) ] )

    @staticmethod
    def parseXMLNode( referenceElement, linkData={} ):
        from fudge.gnd import link
        xlink = link.parseXMLNode( referenceElement ).path
        ref = referenceComponent( None )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

class unknownComponent( component ) :
    """This class is needed for particles that have multiplicity data but no distribution data."""

    def __init__( self ) :

        component.__init__( self, unknownComponentToken, noneFormToken )

    @staticmethod
    def parseXMLNode( element, linkData={} ): return unknownComponent()

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        from fudge.legacy.converting import gndToENDF6
        from fudge.core.math.xData import axes
        gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 0, axes.labToken, [] )

def parseXMLNode( distributionsElement, linkData={} ):
    """Translate a <distributions> element from xml."""
    dist = distribution( distributionsElement.get('nativeData') )
    for component in distributionsElement:
        parserClass = {
                angularComponentToken:          fudge.gnd.productData.distributions.angular,
                energyComponentToken:           fudge.gnd.productData.distributions.energy,
                energyAngularComponentToken:    fudge.gnd.productData.distributions.energyAngular,
                angularEnergyComponentToken:    fudge.gnd.productData.distributions.angularEnergy,
                LLNLAngularEnergyComponentToken: fudge.gnd.productData.distributions.angularEnergy.LLNLComponent,
                LLNL_withAngularComponentToken: fudge.gnd.productData.distributions.angularEnergy.LLNL_withAngularComponent,
                LegendreComponentToken:         fudge.gnd.productData.distributions.Legendre.component,
                LegendreEnergyConservationComponentToken:
                                                fudge.gnd.productData.distributions.Legendre.energyConservationComponent,
                uncorrelatedComponentToken:     fudge.gnd.productData.distributions.uncorrelated,
                CoulombElasticComponentToken:   fudge.gnd.productData.distributions.angular.CoulombElasticComponent,
                referenceComponentToken:        referenceComponent,
                unknownComponentToken:          unknownComponent,
                }.get(component.tag)
        if parserClass is None:
            raise Exception(" can't handle %s distributions yet" % component.tag)
        try:
            newComponent = parserClass.parseXMLNode( component, linkData )
        except Exception as e:
            raise Exception, "%s%s" % (dist.toXLink(), e)
        dist.addComponent( newComponent )
    return dist
