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

from fudge.core.utilities import brb

from fudge.core.utilities import fudgeExceptions
from pqu import PQU

import fudge
from . import base as baseModule
from fudge.gnd.channelData import Q as QModule

from fudge.gnd.reactionData import availableEnergy as availableEnergyModule
from fudge.gnd.reactionData import availableMomentum as availableMomentumModule

__metaclass__ = type

class reaction( baseModule.base_reaction ) :
    """This is the class for a normal gnd reaction."""

    moniker = 'reaction'

    def __init__( self, outputChannel, label, ENDF_MT, documentation = None, date = None, EFL = None ) :
        """
        Creates a new reaction object. Reaction is two-body or uncorrelated-body, depending on
        the outputChannel type. This class is only meant to be used for 'distinct' reactions (distinct reactions
        do not overlap any other reactions; altogether they sum up to total).
        To store a sum over these distinct reactions, use the summedReaction class instead.
        """

        baseModule.base_reaction.__init__( self, label, ENDF_MT, documentation, date = date, EFL = EFL )

# ?????? BRB, I do not think this next line will work anymore. I think outputChannel must always be an instance of fudge.gnd.channels.channel.
        if( type( outputChannel ) == type( '' ) ) : outputChannel = fudge.gnd.channels.channel( 'output', outputChannel )
        if( not isinstance( outputChannel, fudge.gnd.channels.channel ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'Input channel not instance of class channel.' )
        outputChannel.setAncestor( self )
        self.outputChannel = outputChannel

        self.availableEnergy = availableEnergyModule.component( )
        self.availableEnergy.setAncestor( self )

        self.availableMomentum = availableMomentumModule.component( )
        self.availableMomentum.setAncestor( self )

    def __eq__(self, other):
        if (not baseModule.isGNDReaction( other )): return False
        selfParent, otherParent = self.getReactionSuite( ), other.getReactionSuite( )
        return( ( selfParent.projectile == otherParent.projectile ) and ( selfParent.target == otherParent.target ) 
            and ( self.outputChannel == other.outputChannel ) )
    
    def __cmp__( self, other ) :
        """Test if self is <, == or > other."""

        if( not baseModule.isGNDReaction( other ) ) : raise fudgeExceptions.FUDGE_Exception( "Other not an reaction object." )
        selfParent, otherParent = self.getReactionSuite( ), other.getReactionSuite( )
        if( selfParent.projectile < otherParent.projectile ) : return( -1 )
        if( selfParent.projectile > otherParent.projectile ) : return(  1 )
        if( selfParent.target < otherParent.target ) : return( -1 )
        if( selfParent.target > otherParent.target ) : return(  1 )
        if( self.outputChannel < other.outputChannel ) : return( -1 )
        if( self.outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def __str__( self ) :

        return( str( self.outputChannel ) )

    def getQ( self, unit, final = True, groundStateQ = False ) :        # ????? unit for Q is needed, badly.
        """Returns the Q-value for this input/output channel. It will be converted to a float if possible, otherwise a string value is returned."""

        def getLevel( particle, unit ) :

            if( hasattr( particle, 'getLevelAsFloat' ) ) : return( particle.getLevelAsFloat( unit = unit ) )
            return( 0. )

        def getFinalProductsLevelMultiplicity( particle, final ) :
            """Returns a list of [ productName, multiplicity ] for the final products of particle."""

            if( ( final == False ) or ( particle.decayChannel is None ) ) :
                if( particle.particle.name == 'gamma' ) :
                    pms = [ [ particle, 0., 0. ] ]
                else :
                    pms = [ [ particle, particle.multiplicity.getConstant( ), getLevel( particle, unit ) ] ]
            else :
                pms = []
                for product in particle.decayChannel :
                    pms += getFinalProductsLevelMultiplicity( product, final )
            return( pms )

        return( self.outputChannel.getConstantQAs( unit, final = final ) )
        Q = self.Q
        if( not( Q is None ) ) :
                if( isinstance( Q, PQU.PQU ) ) :
                    Q = Q.getValueAs( unit )
                else :
                    pass                                        # ????? This needs work. For example, check for energy dependent Q.
        else :
            pm = {}
            for particle in [ self.getReactionSuite( ).projectile, self.getReactionSuite( ).target ] :
                if( particle.name not in pm ) : pm[particle.name] = [ particle, 0, 0. ]
                pm[particle.name][1] += 1
                pm[particle.name][2] += getLevel( particle, unit )
            for particle in self.outputChannel :
                pms = getFinalProductsLevelMultiplicity( particle, final = final )
                for p, m, level in pms :
                    if( p.particle.name not in pm ) : pm[p.particle.name] = [ p, 0, 0. ]
                    pm[p.particle.name][1] -= m
                    pm[p.particle.name][2] -= level
            Q = 0.
            levelQ = 0.
            try :
                for p in pm :
                    Q += pm[p][1] * pm[p][0].getMass( 'MeV/c**2' )  # Unit 'MeV' should not be hard-wired.
                    levelQ += pm[p][2]
            except :
                print p, pm[p]
                raise
            if( groundStateQ ) : levelQ = 0.
            Q = Q + levelQ
        return( Q )

    def getThreshold( self, unit ) :

        Q = self.getQ( unit = unit, final = False )
        if( Q >= 0. ) : return( 0. )
        return( -Q * ( 1. + self.getReactionSuite( ).projectile.getMass( 'amu' ) / self.getReactionSuite( ).target.getMass( 'amu' ) ) )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False, 
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) :

        return( self.crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin, 
            heatBelowThreshold, heatAllEDomain ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :
        """
        Calculate average energy deposited to the outgoing product.

        :param processInfo: processingInfo.processInfo instance with style = gnd.styles.averageProductData instance
        :tempInfo: dictionary for storing intermediary data.
        :param verbosityIndent: string, indentation level
        :return:
        """

        if( isinstance( self.outputChannel, fudge.gnd.channels.sumOfRemainingOutputChannels ) ) : return    # BRB FIXME - Why is this?
        if( processInfo.verbosity >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        tempInfo['outputChannel'] = self
        tempInfo['EMin'], tempInfo['EMax'] = self.domain( )
        tempInfo['incidentEnergyUnit'] = self.crossSection.domainUnit( )
        for productIndex, product in enumerate( self.outputChannel ) : 
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productToken'] = product.name
            tempInfo['productLabel'] = product.label
            product.calculateDepositionData( processInfo, tempInfo, verbosityIndent + processInfo.verbosityIndentStep )

    def check( self, info ) :

        warnings = self.__checkCrossSection__( info )
        if not info['crossSectionOnly']:
            warnings += self.__checkOutputChannel__( info )
        return warnings

    def isBasicReaction( self ) :

        return( True )

    def isCompleteReaction( self ) :

        return( True )

    def processSn( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.gnd import xParticle

#        if( isinstance( self.outputChannel, fudge.gnd.channels.sumOfRemainingOutputChannels ) ) : return
        crossSection = self.crossSection
# The next lines are probably not needed as grouping has been updated to handle grouping with domain is outside group boundaries. BRB Sept. 2014
#        if( 'LLNL_Pn' in processInfo['styles'] ) :
#            if( processInfo['particles'][tempInfo['reactionSuite'].projectile.name].groups[-1] <= crossSection.domainMin( ) ) : return

        tempInfo['outputChannel'] = self
        if( processInfo.verbosity >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        if( 'LLNL_Pn' in processInfo['styles'] ) :
            fluxl0 = processInfo['flux']['data'][0]
            tempInfo['groupedFlux'] = fluxl0.groupOneFunction( processInfo.getParticleGroups( projectile ) )
        tempInfo['crossSection'] = crossSection.processSn( processInfo, tempInfo, verbosityIndent + '    ' )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            try :
                availableEnergy = self.getData( fudge.gnd.reactionData.base.availableEnergyToken )
            except :
                availableEnergy = fudge.gnd.reactionData.availableEnergy.component( )
                self.addData( availableEnergy )
            availableEnergy.makeGrouped( processInfo, tempInfo )

            try :
                availableMomentum = self.getData( fudge.gnd.reactionData.base.availableMomentumToken )
            except :
                availableMomentum = fudge.gnd.reactionData.availableMomentum.component( )
                self.addData( availableMomentum )
            availableMomentum.makeGrouped( processInfo, tempInfo )

            if( isinstance( self.outputChannel.Q, QModule.pointwise ) ) :
                self.getData( fudge.gnd.channelData.energyDependentQToken ).makeGrouped( processInfo, tempInfo )    # ????? This is not right.

            tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) +  \
                self.outputChannel.toString( simpleString = True )

        processInfo['workFile'].append( 'c%s' % self.label )
        for productIndex, product in enumerate( self.outputChannel ) : 
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productToken'] = product.name
            tempInfo['productGroupToken'] = tempInfo['productToken']
            if( isinstance( product.particle, xParticle.nuclearLevel ) ) : tempInfo['productGroupToken'] = product.particle.groundState.name
            tempInfo['productLabel'] = product.label
            product.processSn( processInfo, tempInfo, verbosityIndent + '    ' )

        del processInfo['workFile'][-1]

    def toString( self, indent = '' ) :

        s, p = indent, ''
        for particle in self.outputChannel :
            s += p + particle.toString( )
            p = ' + '
        return( s + '\n' )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
        """Translate a <reaction> element from xml into a reaction class instance."""

        xPath.append( '%s[@label="%s"]' % ( element.tag, element.get( 'label' ) ) )

        crossSectionComponent = fudge.gnd.reactionData.crossSection.parseXMLNode( 
                element.find( 'crossSection' ), xPath, linkData )

        outputChannel = fudge.gnd.channels.parseXMLNode( element.find( 'outputChannel' ), xPath, linkData )
        outputChannel.fissionGenre = element.get( 'fissionGenre' )

        reac = cls( outputChannel=outputChannel, label=element.get('label'), ENDF_MT=int(element.get('ENDF_MT')),
                date=element.get('date') )
        if 'process' in element:
            reac.process = element.get('process')
        for crossSection in crossSectionComponent : reac.crossSection.add( crossSection )

        if( element.find( 'documentations' ) ) :
            for doc in element.find( 'documentations' ) :
                reac.addDocumentation( fudge.gnd.documentation.documentation.parseXMLNode( doc, xPath, linkData ) )

        xPath.pop( )
        return( reac )
