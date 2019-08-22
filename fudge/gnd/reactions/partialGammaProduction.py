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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
This module contains the reaction class.
"""

import fudge
from . import base, reaction

__metaclass__ = type

class partialGammaProduction( reaction.reaction ) :
    """
    The origin of gammas may be unknown, especially at higher incident energies where many reaction channels
    could be responsible for gamma production. This class contains all gammas whose origin is unknown.
    The outputChannel contains only a gamma product. This is similar to C=55 in ENDL, and to MT=3 in ENDF.
    """

    moniker = 'partialGammaProduction'

    def __init__( self, outputChannel, label, ENDF_MT, documentation = None, date = None ) :

        reaction.reaction.__init__( self, outputChannel, label, ENDF_MT, documentation, date = date )

    def __str__( self ):
        return ("%s production" % self.outputChannel)

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


def parseXMLNode( reactionElement, xPath, linkData ) :

    xPath.append( '%s[@label="%s"]' % ( reactionElement.tag, reactionElement.get( 'label' ) ) )
    crossSectionComponent = fudge.gnd.reactionData.crossSection.parseXMLNode( element.find( 'crossSection' ), xPath, linkData )
    outputChannel = fudge.gnd.channels.parseXMLNode( reactionElement[1], xPath, linkData )
    MT = int( reactionElement.get( 'ENDF_MT' ) )
    date = reactionElement.get( 'date' )
    reac = partialGammaProduction( outputChannel, reactionElement.get( 'label' ), MT, date = date )
    for crossSection in crossSectionComponent : reac.crossSection.add( crossSection )
    if( reactionElement.find( 'documentations' ) ) :
        for doc in reactionElement.find( 'documentations' ) :
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode( doc ) )
    xPath.pop( )
    return( reac )
