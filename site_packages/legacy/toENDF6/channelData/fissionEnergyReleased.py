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

from fudge.gnds.channelData import fissionEnergyReleased as fissionEnergyReleasedModule

from .. import endfFormats as endfFormatsModule
from .. import gndsToENDF6 as gndsToENDF6Module

#
# component
#
def toENDF6( self, endfMFList, flags, targetInfo ) :

    order = 0

    energyReleaseTerms = (
        'promptProductKE',      'promptNeutronKE',      'delayedNeutronKE', 'promptGammaEnergy',
        'delayedGammaEnergy',   'delayedBetaEnergy',    'neutrinoEnergy',   'nonNeutrinoEnergy',
        'totalEnergy')

    for term in energyReleaseTerms:
        form = getattr( self, term ).data
        if isinstance( form, fissionEnergyReleasedModule.polynomial1d ):
            order = len( form.coefficients )-1

    n, data = order+1, {}
    tabulatedTerms = []
    for i in range( n ) :
        data[i] = []
        for term in energyReleaseTerms:
            form = getattr( self, term ).data
            if isinstance( form, fissionEnergyReleasedModule.polynomial1d ):
                data[i] += [ form.coefficients[i], form.uncertainty.data.coefficients[i] ]
            elif isinstance( form, fissionEnergyReleasedModule.XYs1d ):
                data[i] += [ form.evaluate(1e-5), 0 ]
                tabulatedTerms.append( term )
            else:
                raise NotImplementedError("Unknown fission energy release form %s" % type(form))
    datalist = []
    for i in range( n ) : datalist += data[i]

    LFC = len(tabulatedTerms) > 0
    endfMFList[1][458] = [ endfFormatsModule.endfContLine( targetInfo['ZA'], targetInfo['mass'], 0, LFC, 0, len(tabulatedTerms) ) ]
    endfMFList[1][458].append( endfFormatsModule.endfContLine( 0, 0, 0, order, 18 * n, 9 * n ) )
    endfMFList[1][458] += endfFormatsModule.endfDataList( datalist )

    for term in tabulatedTerms:
        LDRV = 2 if term == 'promptProductKE' else 1    # FIXME, hard-coded
        IFC = energyReleaseTerms.index(term) + 1
        form = getattr( self, term ).data

        endfInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(form.interpolation)
        flatData = []
        for xy in form.copyDataToXYs(): flatData += xy
        interpolations = [len(flatData) / 2, endfInterpolation]
        endfMFList[1][458].append(endfFormatsModule.endfContLine(0, 0, LDRV, IFC, len(interpolations) / 2,
                                                                                len(flatData) / 2))
        endfMFList[1][458] += endfFormatsModule.endfInterpolationList(interpolations)
        endfMFList[1][458] += endfFormatsModule.endfDataList(flatData)

    endfMFList[1][458].append( endfFormatsModule.endfSENDLineNumber( ) )


fissionEnergyReleasedModule.fissionEnergyReleased.toENDF6 = toENDF6
