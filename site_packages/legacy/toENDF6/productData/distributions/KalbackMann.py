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

import xData.XYs as XYsModule
import pqu.PQU as PQUModule

from fudge.core.utilities import brb
import fudge.gnd.productData.distributions.KalbachMann as KalbachMannModule

import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module

#
# KalbachMann
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    fSubform = self.fSubform.data
    rSubform = self.rSubform.data
    aSubform = self.aSubform
    if( not( aSubform.isEmptyASubform( ) ) ) : raise 'hell - FIXME'

    outgoingInterpolation = set( [val.interpolation for val in fSubform] + [val.interpolation for val in rSubform] )
    if len(outgoingInterpolation) != 1:
        raise NotImplementedError("Only one outgoing interpolation supported when writing Kalbach-Mann to ENDF-6")

    LEP = { 'flat': 1, 'lin,lin': 2 }[outgoingInterpolation.pop()]
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 2, LEP, 1, len( fSubform ) ) ]
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( fSubform ), 2 ] )
    EInUnit = fSubform.axes[2].unit
    EpUnit = fSubform.axes[1].unit
    fUnit = fSubform.axes[0].unit
    EInFactor = PQUModule.PQU( 1, EInUnit ).getValueAs( 'eV' )
    EpFactor = PQUModule.PQU( 1, EpUnit ).getValueAs( 'eV' )
    fFactor = PQUModule.PQU( 1, fUnit ).getValueAs( '1/eV' )
    for i1, fAtEnergy in enumerate( fSubform ) :
        rAtEnergy = rSubform[i1]
        if( not( isinstance( fAtEnergy, XYsModule.XYs ) ) ) : raise 'hell - FIXME'
        if( not( isinstance( rAtEnergy, XYsModule.XYs ) ) ) : raise 'hell - FIXME'
        value = fAtEnergy.value
        if( value != rAtEnergy.value ) : raise 'hell - FIXME'
        value *= EInFactor
        coefficients = []
        for i1, ( Ep1, f1 ) in enumerate( fAtEnergy ) :
            Ep2, r1 = rAtEnergy[i1]
            if( Ep1 != Ep2 ) : raise 'hell - FIXME'
            coefficients += [ EpFactor * Ep1, fFactor * f1, r1 ]
        length = len( coefficients )
        ENDFDataList += [ endfFormatsModule.endfContLine( 0, value, 0, 1, length, length / 3 ) ]
        ENDFDataList += endfFormatsModule.endfDataList( coefficients )
    gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 1, self.productFrame, ENDFDataList )

KalbachMannModule.form.toENDF6 = toENDF6
