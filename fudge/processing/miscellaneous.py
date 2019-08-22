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

def _groupFunctionsAndFluxInit( projectile, processInfo, f1 ) :

    groupBoundaries = processInfo.getParticleGroups( projectile )
    flux = processInfo['flux']['data'][0]
    domainMin, domainMax = f1.domain( )
    return( groupBoundaries, flux.domainSlice( domainMin, domainMax ) )

def _mutualifyGrouping3Data( f1, f2, f3 ) :

    domainMin1, domainMax1 = f1.domain( )
    domainMin2, domainMax2 = f2.domain( )
    domainMin3, domainMax3 = f2.domain( )
    if( ( domainMin1 != domainMin2 ) or ( domainMin1 != domainMin3 ) or ( domainMax1 != domainMax2 ) or ( domainMax1 != domainMax3 ) ) :
        domainMin, domainMax = max( domainMin1, domainMin2, domainMin3 ), min( domainMax1, domainMax2, domainMax3 )
        print "WARNING: making domains mutual for grouping,", domainMin1, domainMin2, domainMin3, domainMax1, domainMax2, domainMax3
        f1 = f1.domainSlice( domainMin = domainMin, domainMax = domainMax )
        f2 = f2.domainSlice( domainMin = domainMin, domainMax = domainMax )
        f3 = f3.domainSlice( domainMin = domainMin, domainMax = domainMax )
    return( f1, f2, f3 )

def groupOneFunctionAndFlux( projectile, processInfo, f1, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux_, norm = norm ) )

def groupTwoFunctionsAndFlux( projectile, processInfo, f1, f2, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    f1, f2, flux_ = _mutualifyGrouping3Data( f1, f2, flux_ )
    return( f1.groupThreeFunctions( groupBoundaries, flux_, f2, norm = norm ) )
