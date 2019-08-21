# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

momentumDepositionForms = [ tokens.constantFormToken, tokens.partialProductionFormToken, tokens.pointwiseFormToken, tokens.groupedFormToken ]

#
# momentumDeposition genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.momentumDepositionToken

    def __init__( self, nativeData ) :

        baseClasses.componentBase.__init__( self, momentumDepositionForms, nativeData )

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

class pointwise( base.XYPointwiseFormBase ) :

    genre = component.genre

    def __init__( self, axes, data, accuracy, **kwargs ) :

        base.XYPointwiseFormBase.__init__( self, axes, data, accuracy, **kwargs )
        self.toForms = { tokens.groupedFormToken : grouped }

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, momentumDepositionName = 'momentumDeposition',
        momentumDepositionInterpolation = axes.linearToken, momentumDepositionUnit = 'eV/c' ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, frame = axes.labToken, \
            interpolation = axes.interpolationXY( energyInterpolation, momentumDepositionInterpolation ) )
        axes_[1] = axes.axis( momentumDepositionName, 1, momentumDepositionUnit, frame = axes.labToken )
        return( axes_ )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <depositionMomentum> component element into fudge, including all forms in the component."""
    dmc = component( element.get('nativeData') )
    for form in element:
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.groupedFormToken: grouped,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown depositionMomentum form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/depositionMomentum/%s: %s" % (form.tag, e)
        dmc.addForm( newForm )
    return dmc
