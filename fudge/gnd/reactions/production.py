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
    """
    Production reactions are another special case of the <reaction> class, used to store the cross section for
    producing specific radioactive products. Requires a cross section and product id, but no distributions.
    As with the summedReaction, 'outputChannel' should be a channels.simpleOutputChannel instance.
    """

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, documentation = None, attributes = {} ) :
        """Creates a new production reaction object."""

        base.base_reaction.__init__( self, base.productionToken, label, ENDF_MT, crossSection, documentation, attributes )

        if( not isinstance( outputChannel, fudge.gnd.channels.channel ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'Input channel not instance of class channel.' )
        outputChannel.setParent( self )
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

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

        from fudge.gnd.reactionData import crossSection
        from fudge.gnd.productData.distributions import base as distribution_base

        if( flags['verbosity'] >= 10 ) : print '%sproduction: %s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        MT = int( self.attributes['ENDF_MT'] )
        numberOfProductionGammas = 0
        for product in self.outputChannel :
            if( product.getName( ) == 'gamma' ) :
                if( product.getDistributionNativeData( ) not in [ distribution_base.noneComponentToken, distribution_base.unknownComponentToken ] ) :
                    numberOfProductionGammas += 1
        if( numberOfProductionGammas > 0 ) :
            outputChannel = self.outputChannel
            if( len( outputChannel ) != 1 ) : raise Exception( 'production only supported for one product: have %d' % len( outputChannel ) )
            productName = outputChannel[0].getName( )
            if( productName != 'gamma' ) : raise Exception( 'production only supported for gamma as product: product is "%s"' % productName )

            MF = 13
            if( isinstance( self.getCrossSection( ).getNativeData( ), crossSection.reference ) ) : MF = 12
            if( MT not in targetInfo['production_gammas'] ) : targetInfo['production_gammas'][MT] = [ MF ]
            if( MF != targetInfo['production_gammas'][MT][0] ) : raise Exception( 'Prior instances were MF = %d, which cannot be mixed with MF = %d' % 
                ( targetInfo['production_gammas'][MT][0], MF ) )
            targetInfo['production_gammas'][MT].append( self )
        else :
            if( MT not in targetInfo['MF8'] ) : targetInfo['MF8'][MT] = []
            targetInfo['MF8'][MT].append( self )

def parseXMLNode( element, xPath=[], linkData={} ):

    xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label')) )

    crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( element.find('crossSection'), xPath, linkData )
    outputChannel = fudge.gnd.channels.parseXMLNode( element.find('outputChannel'), xPath, linkData )
    MT = int( element.get('ENDF_MT') )
    attributes = {'date': element.get('date')}
    reac = production( outputChannel, element.get('label'), MT, crossSection, attributes = attributes )
    if element.find('documentations'):
        for doc in element.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )

    xPath.pop()
    return reac
