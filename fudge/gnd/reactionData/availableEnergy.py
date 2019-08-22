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

availableEnergyForms = [ tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken ]

#
# availableEnergy genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.availableEnergyToken

    def __init__( self ) :

        baseClasses.componentBase.__init__( self, availableEnergyForms )

    def makeGrouped( self, processInfo, tempInfo ) :

        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        energyUnit = tempInfo['crossSection'].getDomainUnit( )
        axes_ = axes.defaultAxes( labelsUnits = { 0 : [ 'energy_in', energyUnit ], 1 : [ 'energy', energyUnit ] } )
        xMin, xMax = tempInfo['crossSection'].getDomain( )
        energy = XYs.XYs( axes_, [ [ xMin, 0. ], [ xMax, xMax ] ], accuracy = 1e-8 )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, energy )
        self.addForm( grouped( axes_, grouped_ ) )
        axes_, grouped_ = miscellaneous.makeGrouped( self, processInfo, tempInfo, energy, normType = 'groupedFlux' )
        self.addForm( groupedWithCrossSection( axes_, grouped_ ) )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedFormBase.__init__( self, axes_, groupData )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes_, groupData )

def parseXMLNode( element, linkData={} ):
    """Reads an xml <availableEnergy> component element into fudge, including all forms in the component."""
    aec = component()
    for form in element:
        formClass = {
                tokens.groupedFormToken: grouped,
                tokens.groupedWithCrossSectionFormToken: groupedWithCrossSection,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown availableEnergy form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, linkData )
        except Exception as e:
            raise Exception, "/availableEnergy/%s: %s" % (form.tag, e)
        aec.addForm( newForm )
    aec.nativeData = element.get('nativeData')
    return aec
