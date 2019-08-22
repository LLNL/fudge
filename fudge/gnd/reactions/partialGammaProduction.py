# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, documentation = None, attributes = {} ) :

        reaction.reaction.__init__( self, outputChannel, label, ENDF_MT, crossSection, documentation, attributes )
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


def parseXMLNode( reactionElement, xPath=[], linkData={} ):
    xPath.append( (reactionElement.tag, reactionElement.get('label')) )
    crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( reactionElement[0], xPath, linkData )
    outputChannel = fudge.gnd.channels.parseXMLNode( reactionElement[1], xPath, linkData )
    MT = int( reactionElement.get('ENDF_MT') )
    attributes = {'date': reactionElement.get('date')}
    reac = partialGammaProduction( outputChannel, reactionElement.get('label'), MT, crossSection, attributes = attributes )
    if reactionElement.find('documentations'):
        for doc in reactionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )
    xPath.pop()
    return reac
