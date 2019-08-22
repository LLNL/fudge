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

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from xData import axes as axesModule
from xData import regions as regionsModule
from xData import standards as standardsModule

import fudge.gnd.productData.distributions.uncorrelated as uncorrelatedModule
import fudge.gnd.productData.distributions.energy as energyModule
import fudge.gnd.productData.distributions.angular as angularModule
import fudge.gnd.productData.distributions.energyAngular as energyAngularModule

from ... import endfFormats as endfFormatsModule
from ... import gndToENDF6 as gndToENDF6Module

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :
    """
    In ENDF MF=6, some distributions should really be treated as uncorrelated: NBodyPhaseSpace, and also 
    Legendre expansions when only L=0 is listed.
    For GND we split these into uncorrelated angular (isotropic) and energy distributions.
    Must put back in original format when writing back to ENDF.
    """

    frame = self.productFrame
    energySubform = self.energySubform.data
    angularSubform = self.angularSubform.data
    if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
        energySubform.toENDF6( MT, endfMFList, flags, targetInfo )
    elif( targetInfo['doMF4AsMF6'] ) :
        if( isinstance( energySubform, ( energyModule.discreteGamma, energyModule.primaryGamma ) ) ) :
# BRB - this needs to be checked.
            if( targetInfo['product'].id != IDsPoPsModule.photon ) : raise ValueError( 'This logic is only for discrete gammas' )
            energyForm = energyModule.form( self.label, frame, energySubform )
            angularForm = angularModule.form( self.label, frame, angularSubform )
            energyForm.toENDF6( MT, endfMFList, flags, targetInfo )
            angularForm.toENDF6( MT, endfMFList, flags, targetInfo )
        elif( isinstance( angularSubform, angularModule.isotropic ) ) :                # Change to energyAngular with Legendre
            if( not( isinstance( energySubform, energyModule.XYs2d ) ) ) : raise 'hell - fix me'
            axes = axesModule.axes( rank = 4 )
            axes[3] = axesModule.axis( 'energy_in', 3, 'eV' )
            axes[2] = axesModule.axis( 'energy_out', 2, 'eV' )
            axes[1] = axesModule.axis( 'l', 1, '' )
            axes[0] = axesModule.axis( 'C_l(energy_out|energy_in)', 0, '1/eV' )
            energyAngularSubform = energyAngularModule.XYs3d( axes = axes, interpolation = energySubform.interpolation,
                    interpolationQualifier = energySubform.interpolationQualifier )
            EInFactor = PQUModule.PQU( 1, energySubform.axes[2].unit ).getValueAs( 'eV' )

            for EIn in energySubform :
                if isinstance( EIn, regionsModule.regions1d ):
                    # writing to MF6 LAW=1, LANG=1 doesn't support multiple regions, must recombine.
                    if( len( set( [ ein.interpolation for ein in EIn ] ) ) != 1 ) :
                        raise NotImplemented, "ENDF MF6 LAW=1 LANG=1 doesn't support multiple E' interpolations!"
                    xyvals = EIn[0].copyDataToXYs()
                    for region in EIn[1:]:
                        xynew = region.copyDataToXYs()
                        xynew[0][0] *= 1.00000001
                        xyvals.extend( xynew )
                    EIn_copy = EIn[0].copy( )
                    EIn_copy.value = EIn.value
                    EIn_copy.axes = EIn.axes.copy( )
                    EIn_copy.setData( xyvals )
                    EIn = EIn_copy
                multiD_2d = energyAngularModule.XYs2d( value = EIn.value * EInFactor, interpolation = EIn.interpolation )
                EpCls = EIn.copyDataToXYs( )
                for e_out, Cls in EpCls :
                    multiD_2d.append( energyAngularModule.Legendre( [ Cls ], value = e_out ) )
                energyAngularSubform.append( multiD_2d )
            form = energyAngularModule.form( '', self.productFrame, energyAngularSubform )
            if( targetInfo.get( "gammaToENDF6" ) ) : return( form )
            form.toENDF6( MT, endfMFList, flags, targetInfo )
        elif( isinstance( energySubform, energyModule.XYs2d ) and isinstance( angularSubform, angularModule.XYs2d ) ) :

            LANG, LEP = 12, 2
            if( energySubform.interpolation == standardsModule.interpolation.flatToken ) :
                LEP = 1  # interpolation for E_out
                LANG = 11
            elif( energySubform.interpolation in (standardsModule.interpolation.loglinToken,
                standardsModule.interpolation.loglogToken ) ) :
                LANG = 14
            MF6 = [ endfFormatsModule.endfContLine( 0, 0, LANG, LEP, 1, len( energySubform ) ) ]
            EInInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( energySubform.interpolation )
            MF6 +=  endfFormatsModule.endfInterpolationList( [ len( energySubform ), EInInterpolation ] )
            for indexE, EEpP in enumerate( energySubform ) :
                EMuP = angularSubform[indexE]
                if( EEpP.value != EMuP.value ) : raise Exception( "EEpP.value = %s != EMuP.value = %s" % ( EEpP.value, EMuP.value ) )
                NA, NEP = 2 * len( EMuP ), len( EEpP )
                MF6.append( endfFormatsModule.endfContLine( 0, EEpP.value, 0, NA, NEP * ( NA + 2 ), NEP ) )
                data = []
                for EpP in EEpP :
                    data = [ EpP[0], EpP[1] ]
                    for muP in EMuP : data += muP
                    MF6 += endfFormatsModule.endfDataList( data )
            LAW = 1
            gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
        else :
            raise Exception( 'uncorrelated.toENDF6 not supported for energy subform = %s and angular subform = %s' %
                ( energySubform.moniker, angularSubform.moniker ) )
    else :                          # original data is in uncorrelated form
        if( MT not in [ 527, 528 ] ) :
            angularForm = angularModule.form( "", frame, self.angularSubform.data )
            angularForm.toENDF6( MT, endfMFList, flags, targetInfo )
        energyForm = energyModule.form( "", frame, self.energySubform.data )
        energyForm.toENDF6( MT, endfMFList, flags, targetInfo )
        if( MT == 527 ) : endfMFList[26][MT][0]

uncorrelatedModule.form.toENDF6 = toENDF6
