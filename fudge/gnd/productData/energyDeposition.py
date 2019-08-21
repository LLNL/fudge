# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

energyDepositionForms = [ tokens.pointwiseFormToken, tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken ]

#
# energyDeposition genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.energyDepositionToken

    def __init__( self, nativeData ) :

        baseClasses.componentBase.__init__( self, energyDepositionForms, nativeData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        keys = self.forms.keys( )
        for form in keys :
            if( form == tokens.pointwiseFormToken ) :
                ps = self.forms[form].process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedFormBase.__init__( self, axes, data )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes, data )

class pointwise( base.XYPointwiseFormBase ) :

    genre = component.genre

    def __init__( self, axes, data, accuracy, **kwargs ) :

        base.XYPointwiseFormBase.__init__( self, axes, data, accuracy, **kwargs )
        self.toForms = { tokens.groupedFormToken : grouped, tokens.groupedWithCrossSectionFormToken : groupedWithCrossSection }

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, energyDepositionName = 'energyDeposition',
        energyDepositionInterpolation = axes.linearToken, energyDepositionUnit = None ) :

        if( energyDepositionUnit is None ) : energyDepositionUnit = energyUnit
        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, frame = axes.labToken, \
            interpolation = axes.interpolationXY( energyInterpolation, energyDepositionInterpolation ) )
        axes_[1] = axes.axis( energyDepositionName, 1, energyDepositionUnit, frame = axes.labToken )
        return( axes_ )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <depositionEnergy> component element into fudge, including all forms in the component."""
    dec = component( element.get('nativeData') )
    for form in element:
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.groupedFormToken: grouped,
                tokens.groupedWithCrossSectionFormToken: groupedWithCrossSection,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown depositionEnergy form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/depositionEnergy/%s: %s" % (form.tag, e)
        dec.addForm( newForm )
    return dec
