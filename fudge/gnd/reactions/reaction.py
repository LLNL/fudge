# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the reaction class.
"""

from pqu import PQU

import fudge
from fudge.core.utilities import brb
from fudge.core.utilities import fudgeExceptions

from fudge.gnd.channelData import Q as QModule
from fudge.gnd.reactionData import availableEnergy as availableEnergyModule
from fudge.gnd.reactionData import availableMomentum as availableMomentumModule

from . import base as baseModule

__metaclass__ = type

class reaction( baseModule.base_reaction ) :
    """This is the class for a normal gnd reaction."""

    moniker = 'reaction'

    def __init__( self, outputChannel, label, ENDF_MT, documentation = None, date = None, EFL = None ) :
        """
        Creates a new reaction object. Reaction is two-body or uncorrelated-body, depending on
        the outputChannel type. This class is only meant to be used for 'distinct' reactions (distinct reactions
        do not overlap any other reactions; altogether they sum up to total).
        To store a sum over these distinct reactions, use the product class.
        """

        baseModule.base_reaction.__init__( self, label, outputChannel, ENDF_MT, documentation, date = date, EFL = EFL )

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
        """
        Returns the Q-value for this input/output channel. It will be converted to a float if possible, 
        otherwise a string value is returned.
        """

        def getLevel( particle, unit ) :

            if( hasattr( particle, 'getLevelAsFloat' ) ) : return( particle.getLevelAsFloat( unit = unit ) )
            return( 0. )

        def getFinalProductsLevelMultiplicity( particle, final ) :
            """Returns a list of [ productName, multiplicity ] for the final products of particle."""

            if( ( final == False ) or ( particle.outputChannel is None ) ) :
                if( particle.particle.name == 'gamma' ) :
                    pms = [ [ particle, 0., 0. ] ]
                else :
                    pms = [ [ particle, particle.multiplicity.getConstant( ), getLevel( particle, unit ) ] ]
            else :
                pms = []
                for product in particle.outputChannel :
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

    def isBasicReaction( self ) :

        return( True )

    def isCompleteReaction( self ) :

        return( True )

    def toString( self, indent = '' ) :

        s, p = indent, ''
        for particle in self.outputChannel :
            s += p + particle.toString( )
            p = ' + '
        return( s + '\n' )
