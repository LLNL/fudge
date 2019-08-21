# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

energyDepositionToken = 'depositionEnergy'
momentumDepositionToken = 'depositionMomentum'
multiplicityToken = 'multiplicity'

class XYPointwiseFormBase( baseClasses.formBase, XYs.XYs ) :

    form = tokens.pointwiseFormToken
    tag = tokens.pointwiseFormToken

    def __init__( self, axes_, data, accuracy, **kwargs ) :

        baseClasses.formBase.__init__( self )
        kwargs['isPrimaryXData'] = True
        XYs.XYs.__init__( self, axes_, data, accuracy, **kwargs )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing import miscellaneous

        forms = []

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
            norm = tempInfo['groupedCrossSectionNorm']
            crossSection = tempInfo['crossSection']
            axes_ = self.axes.copy( )
            axes_[0].interpolation = axes.interpolationXY( axes.linearToken, axes.flatToken )
            grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
            forms = [ self.toForms[ tokens.groupedFormToken ]( axes_, grouped ) ]
            if( tokens.groupedWithCrossSectionFormToken in self.toForms ) :
                norm = tempInfo['groupedFlux']
                axes_ = axes_.copy( )
                axes_[1].label = "%s %s" % ( axes_[1].getLabel( ), crossSection.axes[1].getLabel( ) )
                if( axes_[1].getUnit( ) == '' ) :
                    axes_[1].unit = "%s" % ( crossSection.axes[1].getUnit( ) )
                else :
                    axes_[1].unit = "%s * %s" % ( axes_[1].getUnit( ), crossSection.axes[1].getUnit( ) )
                grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
                forms.append( self.toForms[tokens.groupedWithCrossSectionFormToken]( axes_, grouped ) )
        return( forms )

    def toXMLList( self, indent = "" ) :

        return( XYs.XYs.toXMLList( self, indent = indent, incrementalIndent = '  ', pairsPerLine = 100 ) )
