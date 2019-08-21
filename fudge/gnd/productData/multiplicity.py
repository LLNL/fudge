# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.legacy.converting import endfFormats
from fudge.core.math.xData import axes, regions, XYs
from fudge.core.math.xData import polynomial as polynomial_
from fudge.gnd import tokens

__metaclass__ = type

energyDependent = 'energyDependent'

multiplicityForms = [ tokens.constantFormToken, tokens.partialProductionFormToken, tokens.pointwiseFormToken, tokens.piecewiseFormToken, \
    tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken, tokens.referenceFormToken, tokens.weightedReferenceFormToken, \
    tokens.polynomialFormToken, tokens.unknownFormToken ]

#
# multiplicity genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.multiplicityToken

    def __init__( self, nativeData ) :

        baseClasses.componentBase.__init__( self, multiplicityForms, nativeData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        keys = self.forms.keys( )
        for form in keys :
            if( form == tokens.pointwiseFormToken ) :
                ps = self.forms[form].process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

    def isConstant( self ) :

        return( self.nativeData in [ tokens.constantFormToken, tokens.partialProductionFormToken ] )

    def getConstant( self ) :

        if( self.isConstant()  ) : return( self.getValue( 0 ) )
        raise Exception( 'multiplicity type = %s does not have a single value' % self.nativeData )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.forms[self.nativeData].domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.forms[self.nativeData].domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.forms[self.nativeData].getEnergyLimits( EMin, EMax ) )
        
    def getValue( self, E ) :

        return( self.forms[self.nativeData].getValue( E ) )

    def check( self, info ):

        from fudge.gnd import warning
        warnings = []

        if self.isConstant():
            if self.getConstant() < 1:
                warnings.append( warning.negativeMultiplicity( self.getConstant(), self.toXLink() ) )
        else:   # energy-dependent mult.
            for form in self.forms.values():
                if hasattr(form, 'getDomain') and form.getDomain() != info['crossSectionDomain']:
                    # domains aren't identical, but may still be 'mutual' if 1st multiplicity point is 0:
                    if ( (form.getDomain()[1] != info['crossSectionDomain'][1])  or
                            (form.getValue( form.getDomain()[0]) != 0) ):
                        warnings.append( warning.domain_mismatch( 
                            *(form.getDomain() + info['crossSectionDomain']), obj=form ) )
                if hasattr(form, 'yMin') and form.yMin() < 0:
                    warnings.append( warning.negativeMultiplicity( form.yMin(), obj=form ) )
                elif isinstance(form, weightedReference) and form.weights.yMin() < 0:
                    warnings.append( warning.negativeMultiplicity( form.weights.yMin(), obj=form ) )

        return warnings

    def getXMLAttribute( self ) :

        return( self.forms[self.nativeData].getXMLAttribute( ) )

    def toPointwiseLinear( self, lowerEps = 1.e-8, upperEps = 1.e-8, template = None ) :
        """This method calls the toPointwiseLinear method for the 'nativeData' to return the multiplicity as linear-linear pointwise form.
        If nativeData is constantFormToken then template (typically a cross section component) is needed to get EMin, EMax and axes information."""

        if( self.nativeData in [ tokens.constantFormToken, tokens.polynomialFormToken ] ) :
            return( self.forms[self.nativeData].toPointwiseLinear( lowerEps, upperEps, template ) )
        else :
            return( self.forms[self.nativeData].toPointwiseLinear( lowerEps, upperEps ) )

    def toXMLList( self, indent = "" ) :

        if( self.getNativeDataToken( ) in [ tokens.constantFormToken, tokens.partialProductionFormToken ] ) : return( [] )
        return( baseClasses.componentBase.toXMLList( self, indent ) )

    def endfMultiplicityList( self, tempInfo ) :

        if( hasattr( self.getFormByToken( self.getNativeDataToken( ) ), 'endfMultiplicityList' ) ) :
            return( self.getFormByToken( self.getNativeDataToken( ) ).endfMultiplicityList( tempInfo ) )
        else :
            raise Exception( 'Unsupport multiplicity nativeData = %s in method endfMultiplicityList' % self.nativeData )

class unknown( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.unknownFormToken

    def __init__( self ) :

        baseClasses.formBase.__init__( self )

    def getValue( self, E ) :

        return( 0. )

    def getXMLAttribute( self ) :

        return( tokens.unknownFormToken )

    def toPointwiseLinear( self, lowerEps, upperEps ) :

        raise Exception( 'Linear-linear pointwise data for multiplicity form %s does not make sense' % tokens.unknownFormToken )

    def toXMLList( self, indent ) :

        return( [] )

class constant( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.constantFormToken

    def __init__( self, value ) :

        baseClasses.formBase.__init__( self )
        self.value = value

    def endfMultiplicityList( self, tempInfo ) :

        endfMult = [ endfFormats.endfInterpolationLine( [ 2, 2 ] ) ]
        endfMult.append( endfFormats.endfDataLine( [ tempInfo['EMin'], self.value, tempInfo['EMax'], self.value ] ) )
        return( 2, endfMult )

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def getValue( self, E ) :

        return( self.value )

    def getXMLAttribute( self ) :

        return( self.value )

    def toPointwiseLinear( self, lowerEps, upperEps, template  ) :
        """This method returns a the multiplicity as linear-linear pointwise data which spans template's domain."""

        axes_ = pointwise.defaultAxes( energyUnit = template.getDomainUnit( ) )
        return( pointwise( axes_, [ [ template.domainMin( ), self.value ], [ template.domainMax( ), self.value ] ], 1e-12 ) )

    def toXMLList( self, indent ) :

        return( [] )

class pointwise( base.XYPointwiseFormBase ) :

    genre = component.genre

    def __init__( self, axes, data, accuracy, **kwargs ) :

        base.XYPointwiseFormBase.__init__( self, axes, data, accuracy, **kwargs )
        self.toForms = toForms = { tokens.groupedFormToken : grouped, tokens.groupedWithCrossSectionFormToken : groupedWithCrossSection }

    def endfMultiplicityList( self, tempInfo ) :

        nPoints = len( self )
        endfMult = [ endfFormats.endfInterpolationLine( [ nPoints, endfFormats.twoAxesToENDFInterpolation( self.axes, 0 ) ] ) ]
        endfMult += endfFormats.endfNdDataList( self )
        return( nPoints, endfMult )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.getDomain( ) )

    def getValue( self, E ) :

        if( E < self[0][0] ) : E = self[0][0]
        if( E > self[-1][0] ) : E = self[-1][0]
        return( XYs.XYs.getValue( self, E ) )

    def toENDF6List( self ) :

        return( ( 1, ) + self.endfMultiplicityList( None ) )

    def getXMLAttribute( self ) :

        return( energyDependent )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, multiplicityName = base.multiplicityToken, 
        multiplicityInterpolation = axes.linearToken ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, frame = axes.labToken, \
            interpolation = axes.interpolationXY( energyInterpolation, multiplicityInterpolation ) )
        axes_[1] = axes.axis( multiplicityName, 1, '', frame = axes.labToken )
        return( axes_ )

class piecewise( baseClasses.formBase, regions.regionsXYs ) :

    genre = component.genre
    form = tokens.piecewiseFormToken
    moniker = tokens.piecewiseFormToken

    def __init__( self, axes ) :

        baseClasses.formBase.__init__( self )
        regions.regionsXYs.__init__( self, axes )

    def endfMultiplicityList( self, tempInfo ) :

        multData = self.data
        nPoints = len( multData )
        endfMult = [ endfFormats.endfInterpolationLine( [ nPoints, 2 ] ) ]
        endfMult += endfFormats.endfNdDataList( multData.data )
        return( nPoints, endfMult )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ self[0][0][0], self[-1][-1][0] ] )

    def getValue( self, E ) :

        if( E < self[0][0][0] ) : E = self[0][0][0]
        for region in self :
            if( E < region[-1][0] ) : break
        return( region.getValue( E ) )

    def getXMLAttribute( self ) :

        return( energyDependent )

    def toENDF6List( self ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        interpolationFlatData, multiplicityFlatData = [], []
        counter = 0
        lastX, lastY = None, None
        for region in self :
            ENDFInterpolation = gndToENDF6.axesToEndfInterpolationFlag( region.axes )
            data = region.copyDataToXYs( xUnit = 'eV' )
            if( lastX is not None ) :
                if( lastY == data[0][1] ) : data = data[1:]
            counter += len( data )
            interpolationFlatData.append( counter )
            interpolationFlatData.append( ENDFInterpolation )
            for xy in data : multiplicityFlatData += xy
            lastX, lastY = data[-1]
        return( interpolationFlatData, counter, endfFormats.endfDataList( multiplicityFlatData ) )

    def toPointwiseLinear( self, lowerEps, upperEps ) :
        """See regionsXYs.toPointwiseLinear on the use of lowerEps, upperEps."""

        axes_ = axes.defaultAxes( labelsUnits = { 0 : [ self.axes[0].getLabel( ), self.axes[1].getUnit( ) ], 
                                                  1 : [ self.axes[0].getLabel( ), self.axes[1].getUnit( ) ] } )
        accuracy = 1e-6
        if( len( self ) > 0 ) : accuracy = self[0].getAccuracy( )
        xys = regions.regionsXYs.toPointwiseLinear( self, accuracy, lowerEps, upperEps, removeOverAdjustedPoints = True )
        return( pointwise( xys.axes, xys, accuracy = xys.getAccuracy( ) ) )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedFormBase.__init__( self, axes, data )

    def getXMLAttribute( self ) :

        return( energyDependent )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes, data )

    def getXMLAttribute( self ) :

        return( energyDependent )

class reference( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.referenceFormToken

    def __init__( self, referenceInstance ) :
        
        baseClasses.formBase.__init__( self )
        self.setReference( referenceInstance )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.referenceInstance.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.referenceInstance.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.referenceInstance.multiplicity.getEnergyLimits( EMin, EMax ) )

    def getValue( self, E ) :

        return( self.referenceInstance.multiplicity.getValue( E ) )

    def setReference( self, referenceInstance ) :

        self.referenceInstance = referenceInstance

    def toPointwiseLinear( self, lowerEps, upperEps ) :

        return( self.referenceInstance.multiplicity.toPointwiseLinear( lowerEps, upperEps ) )

    def toXMLList( self, indent = "" ) :

        return( [ '%s<%s xlink:type="simple" xlink:href="%s"/>' % ( indent, self.form, self.referenceInstance.toXLink( ) ) ] )

    def getXMLAttribute( self ) :

        return( self.referenceInstance.multiplicity.getXMLAttribute( ) )

    @staticmethod
    def parseXMLNode( form, linkData={} ):
        from fudge.gnd import link
        """ translate <reference> element from xml """
        xlink = link.parseXMLNode( form ).path
        ref = reference(None)
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

class weightedReference( reference ):

    genre = component.genre
    tag = tokens.weightedReferenceFormToken

    def __init__( self, referenceInstance, weights ) :
        """Multiplicity = reference * weights, where reference is an xlink, weights is an XYs instance."""

        from fudge.core.utilities import brb
        reference.__init__( self, referenceInstance )
        self.form = tokens.weightedReferenceFormToken
        if( not( isinstance( weights, XYs.XYs ) ) ) : raise Exception( 'weights must be an instance of XYs.XYs: it is %s' % brb.getType( weights ) )
        weights.tag = tokens.pointwiseFormToken
        self.weights = weights

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.referenceInstance.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.referenceInstance.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyLimits( self, EMin, EMax ) :

        return( self.weights.getDomain( ) )

    def toPointwiseLinear( self, lowerEps, upperEps ) :

        referenceMultiplicity = self.referenceInstance.multiplicity.toPointwiseLinear( lowerEps, upperEps )
        weights = self.weights.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )
        return( referenceMultiplicity.xSlice( xMin = weights.domainMin( ) ) * weights ) 

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        indent3 = indent2 + '  '
        xml = [ '%s<%s>' % (indent, self.form) ]
        xml.append( '%s<reference xlink:type="simple" xlink:href="%s"/>' % (indent2, self.referenceInstance.toXLink() ) )
        xml.append( '%s<weights nativeData="%s">' % ( indent2, self.weights.tag ) )
        xml += self.weights.toXMLList( indent = indent3, pairsPerLine = 200 )
        xml[-1] += '</weights></%s>' % self.form
        return( xml )

    def toENDF6( self, flags, targetInfo ) :
        "Convert back to Legendre form when going back to ENDF:"

        from fudge.core.math.xData import LegendreSeries
        from fudge.gnd.productData import distributions
        axes_ = distributions.angular.LegendrePointwise.defaultAxes( self.weights.axes[1].frame, axes.linearToken, axes.linearToken )
        form = distributions.angular.LegendrePointwise( axes_ )
        weights = self.weights.copyDataToXYs( xUnit = 'eV' )
        for i, enlc in enumerate( weights ) :
            energy, LegendreCoeffs = enlc
            form[i] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), [ LegendreCoeffs ], index = i, value = energy )
        return( form )

    @staticmethod
    def parseXMLNode( WRelement, linkData={} ):
        from fudge.gnd import link
        xlink = link.parseXMLNode( WRelement[0] ).path
        if WRelement[1].get('nativeData')!=tokens.pointwiseFormToken:
            raise Exception("parsing %s multiplicity weights not yet supported" % WRelement[1].get('nativeData'))
        weights = XYs.XYs.parseXMLNode( WRelement[1][0] )
        ref = weightedReference( None, weights )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

class polynomial( baseClasses.formBase, polynomial_.polynomial ) :

    genre = component.genre
    form = tokens.polynomialFormToken
    tag = tokens.polynomialFormToken

    def __init__( self, axes, coefficients ) :

        baseClasses.formBase.__init__( self )
        polynomial_.polynomial.__init__( self, axes, coefficients )

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyLimits( self, EMin, EMax ) :

        return( [ EMin, EMax ] )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        endfMFList[1][MT]  = [ endfFormats.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 1, 0, 0 ) ]
        endfMFList[1][MT] += [ endfFormats.endfHeadLine( 0, 0, 0, 0, len( self ), 0 ) ]
        endfMFList[1][MT] += endfFormats.endfDataList( self.coefficients ) + [ endfFormats.endfSENDLineNumber( ) ]

    def getXMLAttribute( self ) :

        return( energyDependent )

    def toPointwiseLinear( self, lowerEps, upperEps, template, accuracy = 1e-6 ) :

        xys = polynomial_.polynomial.toPointwiseLinear( self, template.domainMin( ), template.domainMax( ), accuracy )
        return( pointwise( self.axes, xys, accuracy = accuracy ) )

    @staticmethod
    def defaultAxes( energyName = 'energy_in', energyUnit = 'eV', multiplicityName = base.multiplicityToken ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( energyName, 0, energyUnit, frame = axes.labToken, interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[1] = axes.axis( multiplicityName, 1, '', frame = axes.labToken )
        return( axes_ )

class partialProduction( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.partialProductionFormToken

    def __init__( self, data = 1 ) :

        baseClasses.formBase.__init__( self )
        self.data = data

    def getValue( self, E ) :

        return( self.data )

    def getXMLAttribute( self ) :

        return( tokens.partialProductionFormToken )

    def toPointwiseLinear( self, lowerEps, upperEps ) :

        from fudge.core.utilities import brb
        brb.objectoutline( self.data )
        raise 'Hell'

    def toXMLList( self, indent ) :

        return( [] )

def parseXMLNode( multElement, linkData={} ):
    """ translate an xml <multiplicity> element into fudge """

    nativeData = multElement.get('nativeData')
    mult = component( nativeData )
    if nativeData==tokens.unknownFormToken: mult.addForm( unknown() )
    for form in multElement:
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.piecewiseFormToken: piecewise,
                tokens.referenceFormToken: reference,
                tokens.weightedReferenceFormToken: weightedReference,
                tokens.polynomialFormToken: polynomial,
                tokens.groupedFormToken: grouped,
                tokens.groupedWithCrossSectionFormToken: groupedWithCrossSection,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown multiplicity form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception ("/%s %s" % (form.tag, e))
        mult.addForm( newForm )
        newForm.setParent( mult )
    return mult
