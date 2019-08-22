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

from .. import endfFormats as endfFormatsModule
from .. import gndToENDF6 as gndToENDF6Module

from fudge.gnd.differentialCrossSection import CoulombElastic as CoulombElasticModule
from fudge.gnd.productData.distributions import angular as angularModule

#
# CoulombElastic.form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    targetInfo['LIDP'] = self.identicalParticles
    targetInfo['productFrame'] = self.productFrame
    self.data.toENDF6( MT, endfMFList, flags, targetInfo )
    del targetInfo['LIDP'], targetInfo['productFrame']

CoulombElasticModule.form.toENDF6 = toENDF6

#
# CoulombExpansion
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    def LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal, interferenceImaginary, lineData ) :

        if LIDP:
            NL = len( nuclear ) - 1
            NW = 3 * NL + 3
        else :
            NL = ( len( nuclear ) - 1 ) // 2
            NW = 4 * NL + 3
        lineData.append( endfFormatsModule.endfContLine( 0, nuclear.value, LTP, 0, NW, NL ) )
        legendreDat = nuclear.coefficients
        for j, r in enumerate( interferenceReal.coefficients ) :
            legendreDat.append( r )
            legendreDat.append( interferenceImaginary[j] )
        lineData += endfFormatsModule.endfDataList( legendreDat )

    counts, interpolationFlagsList, lineData = 0, [], []
    LTP = 1                     # indicates this is a nuclear + interference section

    reactionSuite = targetInfo['reactionSuite']
    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    if hasattr(projectile, 'nucleus'): projectile = projectile.nucleus
    projectileSpin = 0
    if( len( projectile.spin ) > 0 ) : projectileSpin = projectile.spin[0].value
    LIDP = targetInfo['LIDP']
    if( isinstance( self.nuclearTerm.data, angularModule.XYs2d ) ) :
        for ridx in xrange( len( self.nuclearTerm.data ) ) :
            counts += 1
            nuclear, interferenceReal, interferenceImaginary = (self.nuclearTerm.data[ridx],
                                                                self.realInterferenceTerm.data[ridx],
                                                                self.imaginaryInterferenceTerm.data[ridx] )
            LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal, interferenceImaginary, lineData )
        interpolationFlagsList += [ counts, gndToENDF6Module.gndToENDFInterpolationFlag( self.nuclearTerm.data.interpolation ) ]
    elif( isinstance( self.nuclearTerm.data, angularModule.regions2d ) ) :
        for regionIndex, region in enumerate( self.nuclearTerm.data ) :
            interferenceReal, interferenceImaginary = ( self.realInterferenceTerm.data[regionIndex],
                                                        self.imaginaryInterferenceTerm.data[regionIndex] )
            for energyIndex, nuclear in enumerate( region ) :
                if( ( regionIndex != 0 ) and ( energyIndex == 0 ) ) : continue
                counts += 1
                LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal[energyIndex], interferenceImaginary[energyIndex], lineData )
            interpolationFlagsList += [ counts, gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation ) ]
    else :
        raise NotImplementedError( "Unknown data storage inside CoulombExpansion: %s" % type(self.nuclearTerm.data) )
    interpolationFlags = endfFormatsModule.endfInterpolationList( interpolationFlagsList )
    ENDFDataList = [ endfFormatsModule.endfContLine( projectileSpin, 0, LIDP, 0,
            len( interpolationFlagsList ) / 2, counts ) ] + interpolationFlags + lineData
    LAW = 5
    gndToENDF6Module.toENDF6_MF6(MT, endfMFList, flags, targetInfo, LAW, targetInfo['productFrame'], ENDFDataList)

CoulombElasticModule.CoulombExpansion.toENDF6 = toENDF6

#
# NuclearPlusCoulombInterference
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    LI, LTT, MF6 = self.effectiveDistribution.data.toENDF6( flags, { 'doMF4AsMF6' : True } )
    reactionSuite = targetInfo['reactionSuite']
    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    if hasattr(projectile, 'nucleus'): projectile = projectile.nucleus
    projectileSpin = 0
    if( len( projectile.spin ) > 0 ) : projectileSpin = projectile.spin[0].value
    LIDP = targetInfo['LIDP']
    NR, NE = [ int(a) for a in MF6[0].split()[-2:] ]
    MF6[0] = endfFormatsModule.endfContLine( projectileSpin, 0, LIDP, 0, NR, NE )
    LAW = 5
    gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, targetInfo['productFrame'], MF6 )

CoulombElasticModule.NuclearPlusCoulombInterference.toENDF6 = toENDF6
