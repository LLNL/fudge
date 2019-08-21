# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the reaction class.
"""

QNotApplicableToken = 'notApplicable'

from fudge.legacy.converting import endfFormats, endf_endl

import fudge
from . import base, reaction

__metaclass__ = type

class partialGammaProduction( reaction.reaction ) :
    """The origin of gammas may be unknown, especially at higher incident energies where many reaction channels
    could be responsible for gamma production. This class contains all gammas whose origin is unknown.
    The outputChannel contains only a gamma product. This is similar to C=55 in ENDL, and to MT=3 in ENDF."""

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, URRProbabilityTable = None, documentation = None, attributes = {} ) :

        reaction.reaction.__init__( self, outputChannel, label, ENDF_MT, crossSection, URRProbabilityTable, documentation, attributes )
        self.moniker = base.partialGammaProductionToken

    def __str__( self ):
        return self.name

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ) :

        return( False )

    def check( self, info ) :

        from fudge.gnd import warning
        warnings = []

        crossSectionWarnings = self.crossSection.check( info )
        if crossSectionWarnings:
            warnings.append( warning.context("Cross section:", crossSectionWarnings) )
        return warnings


def parseXMLNode( reactionElement, linkData={} ):
    try:
        crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( reactionElement[0], linkData )
    except Exception as e:
        raise Exception, '/reactionSuite/partialGammaProduction[@label="%s"]%s' % (reactionElement.get('label'),e)
    try:
        outputChannel = fudge.gnd.channels.parseXMLNode( reactionElement[1], linkData )
    except Exception as e:
        raise Exception, '/reactionSuite/partialGammaProduction[@label="%s"]%s' % (reactionElement.get('label'),e)
    MT = int( reactionElement.get('ENDF_MT') )
    attributes = {'date': reactionElement.get('date')}
    reac = partialGammaProduction( outputChannel, reactionElement.get('label'), MT, crossSection, attributes = attributes )
    if reactionElement.find('documentations'):
        for doc in reactionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )
    return reac
