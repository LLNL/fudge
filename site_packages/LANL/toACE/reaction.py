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
This module adds the method toACE to the reaction class.
"""

from fudge.gnd.reactions import reaction

def toACE( self, temperature, EMin, data ) :

    from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU

    MT = int( self.attributes['ENDF_MT'] )
    MTData = {}
    MTData['fission'] = self.outputChannel.isFission( )
    if( self.domainMin( ) > EMin ) :
        MTData['ESZ'] = self.crossSection.toPointwise_withLinearXYs( )
    else :
        EMinUnit = PQU( EMin, self.getDomainUnit( ) )
        MTData['ESZ'] = self.heatCrossSection( temperature, EMinUnit, heatBelowThreshold = False, heatAllEDomain = True, interpolationAccuracy = 0.002 )
    MTData['Q'] = self.getQ( 'MeV' )
    MTData['n'] = []
    MTData['gamma'] = []
    for product in self.outputChannel : product.toACE( MTData )
    data.append( ( MT, MTData ) )

reaction.reaction.toACE = toACE
