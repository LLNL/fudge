# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the production class, for storing production cross section
of radioactive products
"""

from fudge.core.utilities import fudgeExceptions
from pqu import physicalQuantityWithUncertainty
from fudge.legacy.converting import endfFormats, endf_endl

import fudge
from . import base

__metaclass__ = type

class production( base.base_reaction ):
    """ Production reactions are another special case of the <reaction> class, used to store the cross section for
    producing specific radioactive products. Requires a cross section and product id, but no distributions.
    As with the summedReaction, 'outputChannel' should be a channels.simpleOutputChannel instance """

    def __init__( self, product, label, ENDF_MT, crossSection = None, Q = None, documentation = None, attributes = {} ) :

        base.base_reaction.__init__( self, base.productionToken, label, ENDF_MT, crossSection, documentation, attributes )
        self.product = product
        self.Q = Q

    def __str__( self ):
        return "%s production" % self.product

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def check( self, info ) :
        warnings = self.__checkCrossSection__( info )
        return warnings

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

        MT = int( self.attributes['ENDF_MT'] )
        if( flags['verbosity'] >= 10 ) : print '%sproduction: %s' % ( verbosityIndent, self.product )
        if( MT not in targetInfo['MF8'] ) : targetInfo['MF8'][MT] = []
        targetInfo['MF8'][MT].append( self )

def parseXMLNode( productionElement, linkData={} ):
    try:
        crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( productionElement.find('crossSection'), linkData )
    except Exception as e:
        raise Exception, '/reactionSuite/production[@label="%s"]%s' % (productionElement.get('label'),e)
    product = linkData['particles'].getParticle( productionElement.get('product') )
    MT = int( productionElement.get('ENDF_MT') )
    Q = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( productionElement.get('Q') )
    Q = fudge.gnd.channelData.Q.component( fudge.gnd.channelData.Q.constant( Q ) )
    attributes = {'date': productionElement.get('date')}
    reac = production( product, productionElement.get('label'), MT, crossSection, Q, attributes = attributes )
    if productionElement.find('documentations'):
        for doc in productionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )
    return reac

