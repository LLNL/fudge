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

def _groupFunctionsAndFluxInit( projectile, processInfo, f1 ) :

    groupBoundaries = processInfo.getParticleGroups( projectile )
    flux = processInfo['flux']['data'][0]
    xMin, xMax = f1.getDomain( )
    return( groupBoundaries, flux.xSlice( xMin, xMax ) )

def _mutualifyGrouping3Data( f1, f2, f3 ) :

    xMin1, xMax1 = f1.getDomain( )
    xMin2, xMax2 = f2.getDomain( )
    xMin3, xMax3 = f2.getDomain( )
    if( ( xMin1 != xMin2 ) or ( xMin1 != xMin3 ) or ( xMax1 != xMax2 ) or ( xMax1 != xMax3 ) ) :
        xMin, xMax = max( xMin1, xMin2, xMin3 ), min( xMax1, xMax2, xMax3 )
        print "WARNING: making domains mutual for grouping,", xMin1, xMin2, xMin3, xMax1, xMax2, xMax3
        f1, f2, f3 = f1.xSlice( xMin = xMin, xMax = xMax ), f2.xSlice( xMin = xMin, xMax = xMax ), f3.xSlice( xMin = xMin, xMax = xMax )
    return( f1, f2, f3 )

def groupOneFunctionAndFlux( projectile, processInfo, f1, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux_, norm = norm ) )

def groupTwoFunctionsAndFlux( projectile, processInfo, f1, f2, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    f1, f2, flux_ = _mutualifyGrouping3Data( f1, f2, flux_ )
    return( f1.groupThreeFunctions( groupBoundaries, flux_, f2, norm = norm ) )
