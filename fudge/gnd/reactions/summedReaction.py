# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the summedReaction class
"""

from pqu import physicalQuantityWithUncertainty
from fudge.legacy.converting import endfFormats, endf_endl

import fudge
from . import base

__metaclass__ = type

class summedReaction( base.base_reaction ):
    """ special case for sums of other reactions, such as (n,total), where outgoing particle descriptions are omitted.
    The 'outputChannel' should be a channels.simpleOutputChannel instance """

    def __init__( self, name, label, ENDF_MT, crossSection = None, Q = None, summands = None, documentation = None, attributes = {} ):

        base.base_reaction.__init__( self, base.summedReactionToken, label, ENDF_MT, crossSection, documentation, attributes )
        self.name = name
        self.Q = Q
        self.summands = summands or []

    def __str__( self ):
        return self.name

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ) :

        return( False )

    def check( self, info ) :
        warnings = self.__checkCrossSection__( info )
        return warnings

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :
        from fudge.legacy.converting import gndToENDF6

        MT = int( self.attributes['ENDF_MT'] )
        if flags['verbosity'] >= 10 : print '%ssummed reaction: %s' % (verbosityIndent, self.name)
        targetInfo['Q'] = self.getQ( 'eV', groundStateQ = True )
        crossSection = self.getCrossSection( )
        formToken = crossSection.nativeData
        targetInfo['EMin'], targetInfo['EMax'] = crossSection.getDomain( )
        if self.summands:
            level = min( [ max( [p.getLevelAsFloat('eV') for p in reac.link.outputChannel])
                for reac in self.summands] )
        else: level = 0
        crossSection.toENDF6( MT, endfMFList, targetInfo, level, LR=0 )

def parseXMLNode( reactionElement, linkData={} ):
    try:
        summands = [fudge.gnd.link.parseXMLNode( link, linkData ) for link in reactionElement.findall('summand')]
        linkData['unresolvedLinks'] += [(sum.path, sum) for sum in summands]
        crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( reactionElement.find('crossSection'), linkData )
    except Exception as e:
        raise Exception, '/reactionSuite/summedReaction[@label="%s"]%s' % (reactionElement.get('label'),e)
    name = reactionElement.get('name')
    MT = int( reactionElement.get('ENDF_MT') )
    Q = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( reactionElement.get('Q') )
    Q = fudge.gnd.channelData.Q.component( fudge.gnd.channelData.Q.constant( Q ) )
    attributes = {'date': reactionElement.get('date')}
    reac = summedReaction( name, reactionElement.get('label'), MT, crossSection, Q, summands, attributes = attributes )
    if reactionElement.find('documentations'):
        for doc in reactionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )
    return reac
