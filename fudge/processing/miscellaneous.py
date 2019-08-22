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

from xData import link as linkModule

def _groupFunctionsAndFluxInit( style, tempInfo, f1 ) :

    reactionSuite = tempInfo['reactionSuite']
    groupBoundaries = style.transportables[reactionSuite.projectile.name].group.boundaries
    domainMin, domainMax = f1.domain( )
    flux = style.flux[0]
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

def groupOneFunctionAndFlux( style, tempInfo, f1 ) :

    if( isinstance( f1, linkModule.link ) ) : f1 = style.findFormMatchingDerivedStyles( f1.link )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux, norm = tempInfo['groupedFlux'] ) )

def groupTwoFunctionsAndFlux( style, tempInfo, f1, f2, norm = None ) :

    groupBoundaries, flux = _groupFunctionsAndFluxInit( tempInfo, f1 )
    f1, f2, flux = _mutualifyGrouping3Data( f1, f2, flux )
    return( f1.groupThreeFunctions( groupBoundaries, flux, f2, norm = norm ) )
