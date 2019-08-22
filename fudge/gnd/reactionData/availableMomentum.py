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

from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import gridded as griddedModule

from fudge.gnd import miscellaneous as miscellaneousModule
from fudge.gnd import tokens as tokensModule
from fudge.gnd import abstractClasses as abstractClassesModule

import base as baseModule

__metaclass__ = type

class baseAvailableMomentumForm( abstractClassesModule.form ) :

    pass
#
# availableMomentum forms
#
class multiGroup( baseAvailableMomentumForm, griddedModule.gridded ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded.__init__( self, **kwargs )
#
# availableMomentum forms
#
class component( abstractClassesModule.component ) :

    moniker = baseModule.availableMomentumToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( multiGroup, ) )

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
        formClass = { tokensModule.groupedFormToken : multiGroup }.get( form.tag )
        if( formClass is None ) : raise Exception( "encountered unknown availableMomentum form: %s" % form.tag )
        try :
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e :
            raise Exception, "availableMomentum/%s: %s" % ( form.tag, e )
        amc.addForm( newForm )
    return amc

def calculateMomentumPoints( massInE, EMin, EMax, energyUnit, accuracy = 1e-4 ) :
    """
    This is a temporary function (hopefully) to calculate momentum vs E.
    What is really needed is a function like rriint in mcfgen.
    """

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
    return( XYsModule.XYs1d( data = momentum, axes = axes, accuracy = accuracy ) )
