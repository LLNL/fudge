# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.gnd import miscellaneous, tokens
from fudge.core.math.xData import axes, XYs

__metaclass__ = type

availableMomentumForms = [ tokens.groupedFormToken ]

#
# availableMomentum genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.availableMomentumToken

    def __init__( self ) :

        baseClasses.componentBase.__init__( self, availableMomentumForms )

    def makeGrouped( self, processInfo, tempInfo ) :

        energyUnit = tempInfo['crossSection'].getDomainUnit( )
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        massInE = tempInfo['reactionSuite'].projectile.getMass( energyUnit + '/c**2' )
        xMin, xMax = tempInfo['crossSection'].getDomain( )
        momentum = calculateMomentumPoints( massInE, xMin, xMax, energyUnit )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, momentum )
        self.addForm( grouped( axes_, grouped_ ) )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedFormBase.__init__( self, axes_, groupData )

def calculateMomentumPoints( massInE, EMin, EMax, energyUnit, accuracy = 1e-4 ) :
    """This is a temporary function (hopefully) to calculate momentum vs E.
    What is really needed is a function like rriint in mcfgen."""

    import math
    momentum = []

    def insertPointIfNeeded( Ep1, Ep2 ) :

        E1, p1 = Ep1
        E2, p2 = Ep2
        s = ( p2 - p1 ) / ( E2 - E1 )
        Em = massInE / ( 2. * s * s )
        pc = math.sqrt( 2 * massInE * Em )
        pi = s * ( Em - E1 ) + p1
        if( abs( pi / pc  - 1 ) > accuracy ) :
            m = [ Em, pc ]
            momentum.append( m )
            insertPointIfNeeded( Ep1, m )
            insertPointIfNeeded( m, Ep2 )

    momentum.append( [ EMin, math.sqrt( 2 * massInE * EMin ) ] )
    momentum.append( [ EMax, math.sqrt( 2 * massInE * EMax ) ] )
    insertPointIfNeeded( momentum[0], momentum[1] )
    momentum.sort( )
    axes_ = axes.defaultAxes( labelsUnits = { 0 : [ 'energy_in', energyUnit ], 1 : [ 'momentum', energyUnit + '/c' ] } )
    return( XYs.XYs( axes_, momentum, accuracy = accuracy ) )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <availableMomentum> component element into fudge, including all forms in the component."""
    amc = component()
    for form in element:
        formClass = {
                tokens.groupedFormToken: grouped,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown availableMomentum form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/availableMomentum/%s: %s" % (form.tag, e)
        amc.addForm( newForm )
    amc.nativeData = element.get('nativeData')
    return amc
