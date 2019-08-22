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

from fudge.gnd import baseClasses
import base
from fudge.gnd import miscellaneous, tokens
from fudge.core.math.xData import axes, XYs

__metaclass__ = type

availableMomentumForms = [ tokens.groupedFormToken ]

#
# availableMomentum genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.availableMomentumToken

    def __init__( self ) :

        baseClasses.componentBase.__init__( self, availableMomentumForms )

    def makeGrouped( self, processInfo, tempInfo ) :

        energyUnit = tempInfo['crossSection'].getDomainUnit( )
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        massInE = tempInfo['reactionSuite'].projectile.getMass( energyUnit + '/c**2' )
        xMin, xMax = tempInfo['crossSection'].getDomain( )
        momentum = calculateMomentumPoints( massInE, xMin, xMax, energyUnit )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, momentum )
        self.addForm( grouped( axes_, grouped_ ) )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedFormBase.__init__( self, axes_, groupData )

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

    momentum.append( [ EMin, math.sqrt( 2 * massInE * EMin ) ] )
    momentum.append( [ EMax, math.sqrt( 2 * massInE * EMax ) ] )
    insertPointIfNeeded( momentum[0], momentum[1] )
    momentum.sort( )
    axes_ = axes.defaultAxes( labelsUnits = { 0 : [ 'energy_in', energyUnit ], 1 : [ 'momentum', energyUnit + '/c' ] } )
    return( XYs.XYs( axes_, momentum, accuracy = accuracy ) )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <availableMomentum> component element into fudge, including all forms in the component."""
    amc = component()
    for form in element:
        formClass = {
                tokens.groupedFormToken: grouped,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown availableMomentum form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/availableMomentum/%s: %s" % (form.tag, e)
        amc.addForm( newForm )
    amc.nativeData = element.get('nativeData')
    return amc
