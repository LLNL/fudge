# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.gnd import miscellaneous
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

energyDependent = 'energyDependent'

QForms = [ tokens.notApplicableFormToken, tokens.constantFormToken, tokens.pointwiseFormToken, tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken ]

#
# Q genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.QToken

    def __init__( self, form = None ) :

        baseClasses.componentBase.__init__( self, QForms )
        if( form is not None ) : self.addForm( form )

    def getConstantAs( self, unit ) :

        if( self.nativeData in [ tokens.constantFormToken ]  ) : return( self.getValue( 0, unit ) )
        raise Exception( 'Q type = %s does not have a single value' % self.nativeData )

    def getValue( self, E, unit ) :

        return( self.forms[self.nativeData].getValue( E, unit ) )

    def getXMLAttribute( self ) :

        return( self.forms[self.nativeData].getXMLAttribute( ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        keys = self.forms.keys( )
        for form in keys :
            if( form == tokens.pointwiseFormToken ) :
                ps = self.forms[form].process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

class constant( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.constantFormToken

    def __init__( self, data ) :

        baseClasses.formBase.__init__( self )
        self.data = data

    def getValue( self, E, unit ) :

        return( self.data.getValueAs( unit ) )

    def getXMLAttribute( self ) :

        return( self.data )

    def toXMLList( self, indent ) :

        return( [] )

class pointwise( baseClasses.formBase, XYs.XYs ) :

    genre = component.genre
    form = tokens.pointwiseFormToken
    tag = tokens.pointwiseFormToken

    def __init__( self, axes, data, accuracy ) :

        baseClasses.formBase.__init__( self )
        XYs.XYs.__init__( self, axes, data, accuracy, isPrimaryXData = True )
        self.toForms = { tokens.groupedFormToken : grouped }

    def process( self, processInfo, tempInfo ) :

        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        Q_E = self.getFormByToken( tokens.pointwiseFormToken )
        grouped = miscellaneous.makeGrouped( self, processInfo, tempInfo, Q_E )
        self.addForm( grouped( grouped ) )
        grouped = miscellaneous.makeGrouped( self, processInfo, tempInfo, Q_E, normType = 'groupedFlux' )
        self.addForm( groupedWithCrossSection( grouped ) )

    def getXMLAttribute( self ) :

        return( energyDependent )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, QUnit = 'eV', QInterpolation = axes.linearToken ) :

        axes_ = axes.axes( )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, frame = axes.labToken, interpolation = axes.interpolationXY( energyInterpolation, QInterpolation ) )
        axes_[1] = axes.axis( base.QToken, 1, QUnit, frame = axes.labToken )
        return( axes_ )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedFormBase.__init__( self, axes, data )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes, data )

class notApplicable( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.notApplicableFormToken

    def __init__( self ) :

        baseClasses.formBase.__init__( self )
        self.data = 'N/A'

    def getValue( self, E, unit ) :

        return( 0. )

    def getXMLAttribute( self ) :

        return( self.data )

    def toXMLList( self, indent ) :

        return( [] )
