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

import base as baseModule
from fudge.gnd import miscellaneous as miscellaneousModule

from fudge.gnd import tokens as tokensModule
from fudge.gnd import abstractClasses as abstractClassesModule

from xData import axes as axesModule
from xData import XYs as XYsModule

__metaclass__ = type

class baseAvailableMomentumForm( abstractClassesModule.form ) :

    __genre = baseModule.availableMomentumToken

#
# availableMomentum forms
#
class grouped( baseAvailableMomentumForm, abstractClassesModule.multiGroup ) :

    def __init__( self, label, axes, data ) :

        baseMomentumDepositionForm.__init__( self )
        abstractClassesModule.multiGroup.__init__( self, label, axes, data )

#
# availableMomentum forms
#
class component( abstractClassesModule.component ) :

    __genre = baseModule.availableMomentumToken
    moniker = baseModule.availableMomentumToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( grouped, ) )

    def makeGrouped( self, processInfo, tempInfo ) :

        energyUnit = tempInfo['crossSection'].domainUnit( )
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        massInE = tempInfo['reactionSuite'].projectile.getMass( energyUnit + '/c**2' )
        domainMin, domainMax = tempInfo['crossSection'].domain( )
        momentum = calculateMomentumPoints( massInE, domainMin, domainMax, energyUnit )
        axes, _grouped = miscellaneousModule.makeGrouped( self, processInfo, tempInfo, momentum )
        self.addForm( grouped( axes, _grouped ) )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <availableMomentum> component element into fudge, including all forms in the component."""

    amc = component()
    for form in element :
        formClass = { tokensModule.groupedFormToken : grouped }.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown availableMomentum form: %s" % form.tag )
        try :
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e :
            raise Exception, "availableMomentum/%s: %s" % ( form.tag, e )
        amc.addForm( newForm )
    return amc

def calculateMomentumPoints( massInE, EMin, EMax, energyUnit, accuracy = 1e-4 ) :
    """This is a temporary function (hopefully) to calculate momentum vs E.
    What is really needed is a function like rriint in mcfgen."""

    import math
    momentum = []

    def insertPointIfNeeded( Ep1, Ep2 ) :

        E1, p1 = Ep1
        E2, p2 = Ep2
        s = ( p2 - p1 ) / ( E2 - E1 )
        Em = massInE / ( 2. * s * s )
        pc = math.sqrt( 2 * massInE * Em )
        pi = s * ( Em - E1 ) + p1
        if( abs( pi / pc  - 1 ) > accuracy ) :
            m = [ Em, pc ]
            momentum.append( m )
            insertPointIfNeeded( Ep1, m )
            insertPointIfNeeded( m, Ep2 )

    if( massInE == 0 ) :            # gammas
        momentum.append( [ EMin, EMin ] )
        momentum.append( [ EMax, EMax ] )
    else :
        momentum.append( [ EMin, math.sqrt( 2 * massInE * EMin ) ] )
        momentum.append( [ EMax, math.sqrt( 2 * massInE * EMax ) ] )
        insertPointIfNeeded( momentum[0], momentum[1] )
    momentum.sort( )
    axes = axesModule.axes( labelsUnits = { 0 : [ 'energy_in', energyUnit ], 1 : [ 'momentum', energyUnit + '/c' ] } )
    return( XYsModule.XYs( data = momentum, axes = axes, accuracy = accuracy ) )
