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

        from fudge.gnd import warning
        warnings = self.__checkCrossSection__( info )

        # does self.crossSection equal the sum over summand reaction cross sections?
        sum_ = self.summands[0].link.crossSection.toPointwise_withLinearXYs()
        for i in range(1,len(self.summands)):
            sum_, current = sum_.mutualify( 1e-8, 1e-8, 0,
                    self.summands[i].link.crossSection.toPointwise_withLinearXYs(), 1e-8, 1e-8, 0 )
            sum_ += current
        quotedXsec = self.crossSection.toPointwise_withLinearXYs()
        if sum_.getDomain() != quotedXsec.getDomain():
            warnings.append( warning.summedCrossSectionDomainMismatch( obj=self ) )
        else:
            diff = quotedXsec - sum_
            diff.setSafeDivide(True)
            relativeDiff = abs( diff ) / quotedXsec
            if relativeDiff.yMax() > info['crossSectionMaxDiff']:
                warnings.append( warning.summedCrossSectionMismatch( relativeDiff.yMax()*100, obj=self ) )

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

def parseXMLNode( reactionElement, xPath=[], linkData={} ):
    xPath.append( '%s[@label="%s"]' % (reactionElement.tag, reactionElement.get('label')) )

    summands = [fudge.gnd.link.parseXMLNode( link, linkData ) for link in reactionElement.findall('summand')]
    linkData['unresolvedLinks'] += [(sum.path, sum) for sum in summands]
    crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( reactionElement.find('crossSection'),
            xPath, linkData )
    name = reactionElement.get('name')
    MT = int( reactionElement.get('ENDF_MT') )
    Q = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( reactionElement.get('Q') )
    Q = fudge.gnd.channelData.Q.component( fudge.gnd.channelData.Q.constant( Q ) )
    attributes = {'date': reactionElement.get('date')}
    reac = summedReaction( name, reactionElement.get('label'), MT, crossSection, Q, summands, attributes = attributes )
    if reactionElement.find('documentations'):
        for doc in reactionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )
    xPath.pop()
    return reac
