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
This module contains the production class, for storing production cross section
of radioactive products
"""

from fudge.core.utilities import fudgeExceptions

import fudge
from . import base as baseModule

__metaclass__ = type

class production( baseModule.base_reaction ):
    """
    Production reactions are another special case of the <reaction> class, used to store the cross section for
    producing specific radioactive products. Requires a cross section and product id, but no distributions.
    As with the summedReaction, 'outputChannel' should be a channels.simpleOutputChannel instance.
    """

    moniker = 'production'

    def __init__( self, outputChannel, label, ENDF_MT, documentation = None, date = None ) :
        """Creates a new production reaction object."""

        baseModule.base_reaction.__init__( self, label, ENDF_MT, documentation, date = date )

        if( not isinstance( outputChannel, fudge.gnd.channels.channel ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'Input channel not instance of class channel.' )
        outputChannel.setAncestor( self )
        self.data = {}
        self.outputChannel = outputChannel

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def check( self, info ) :

        warnings = self.__checkCrossSection__( info )
        return warnings

    def getQ( self, unit, final = True, groundStateQ = False ) :

        return( self.outputChannel.getConstantQAs( unit, final = final ) )

    # FIXME next method nearly identical to reactions.calculateDepositionData. Define in base_reaction?
    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :
        """
        Calculate average energy deposited to the outgoing product.

        :param processInfo: dictionary, contains at least 'verbosity' (int), 'incidentEnergyUnit' (string) and 'reactionSuite'
        :tempInfo: dictionary for storing intermediary data.
        :param verbosityIndent: string, indentation level
        :return:
        """

        if( processInfo.verbosity >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        tempInfo['outputChannel'] = self
        tempInfo['EMin'], tempInfo['EMax'] = self.domain( )
        tempInfo['incidentEnergyUnit'] = self.crossSection.domainUnit( )
        for productIndex, product in enumerate( self.outputChannel ) :
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productToken'] = product.name
            tempInfo['productLabel'] = product.label
            product.calculateDepositionData( processInfo, tempInfo, verbosityIndent + processInfo.verbosityIndentStep )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )

        crossSectionComponent = fudge.gnd.reactionData.crossSection.parseXMLNode( element.find( 'crossSection' ), xPath, linkData )
        outputChannel = fudge.gnd.channels.parseXMLNode( element.find('outputChannel'), xPath, linkData )
        MT = int( element.get('ENDF_MT') )
        date = element.get( 'date' )
        reac = production( outputChannel, element.get('label'), MT, date = date )
        for crossSection in crossSectionComponent : reac.crossSection.add( crossSection )
        if element.find('documentations'):
            for doc in element.find('documentations'):
                reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )

        xPath.pop()
        return reac
