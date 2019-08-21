# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the fissionComponent class
"""

from fudge.legacy.converting import endfFormats, endf_endl

import fudge
from . import base, reaction

__metaclass__ = type

class fissionComponent( reaction.reaction ) :
    """ special case of reaction, for storing 1st, 2nd, 3rd, etc fission chances when we also have total """

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, URRProbabilityTable = None, documentation = None, attributes = {} ) :
        """Creates a new reaction object of type outputChannel (that is, use 2-body outputChannel for 2-body reaction)."""

        reaction.reaction.__init__( self, outputChannel, label, ENDF_MT, crossSection, URRProbabilityTable, documentation, attributes )
        self.moniker = base.fissionComponentToken

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def check( self, info ) :
        warnings = self.__checkCrossSection__( info )
        if not info['crossSectionOnly']:
            warnings += self.__checkOutputChannel__( info )
        return warnings
