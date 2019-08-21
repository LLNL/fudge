# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.gnd import miscellaneous, tokens
from fudge.core.math.xData import axes, XYs

__metaclass__ = type

availableEnergyForms = [ tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken ]

#
# availableEnergy genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.availableEnergyToken

    def __init__( self ) :

        baseClasses.componentBase.__init__( self, availableEnergyForms )

    def makeGrouped( self, processInfo, tempInfo ) :

        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        energyUnit = tempInfo['crossSection'].getDomainUnit( )
        axes_ = axes.defaultAxes( labelsUnits = { 0 : [ 'energy_in', energyUnit ], 1 : [ 'energy', energyUnit ] } )
        xMin, xMax = tempInfo['crossSection'].getDomain( )
        energy = XYs.XYs( axes_, [ [ xMin, 0. ], [ xMax, xMax ] ], accuracy = 1e-8 )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, energy )
        self.addForm( grouped( axes_, grouped_ ) )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, energy, normType = 'groupedFlux' )
        self.addForm( groupedWithCrossSection( axes_, grouped_ ) )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedFormBase.__init__( self, axes_, groupData )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes_, groupData )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <availableEnergy> component element into fudge, including all forms in the component."""
    aec = component()
    for form in element:
        formClass = {
                tokens.groupedFormToken: grouped,
                tokens.groupedWithCrossSectionFormToken: groupedWithCrossSection,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown availableEnergy form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/availableEnergy/%s: %s" % (form.tag, e)
        aec.addForm( newForm )
    aec.nativeData = element.get('nativeData')
    return aec
