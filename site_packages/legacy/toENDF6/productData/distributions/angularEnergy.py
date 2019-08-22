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
from xData import standards as standardsModule

from fudge.gnds.productData.distributions import angularEnergy as angularEnergyModule

from ... import gndsToENDF6 as gndsToENDF6Module
from ... import endfFormats as endfFormatsModule

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    subform = self.angularEnergySubform
    if( hasattr( subform, 'toENDF6' ) ) :
        LAW, frame, MF6 = subform.toENDF6( flags, targetInfo )
        gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
    else :
        print 'WARNING: angularEnergy subform "%s" has no toENDF6 method' % subform.moniker

angularEnergyModule.form.toENDF6 = toENDF6

#
# XYs3d
#
def toENDF6( self, flags, targetInfo ) :

    MF6 = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
    EInInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    energyConversionFactor = PQUModule.PQU(1, self.axes[-1].unit ).getValueAs('eV')
    MF6 += endfFormatsModule.endfInterpolationList( [ len( self ), EInInterpolation ] )
    for oneEin in self :
        muInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
        Ein = oneEin.value * energyConversionFactor
        numMu = len( oneEin )
        MF6 += [ endfFormatsModule.endfContLine( 0, Ein, 0, 0, 1, numMu ) ]
        MF6 += endfFormatsModule.endfInterpolationList( [ numMu, muInterpolation ] )
        for entries in oneEin :
            pdf_of_EpInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( entries.interpolation )
            mu = entries.value
            numEout = len( entries )
            MF6 += [ endfFormatsModule.endfContLine( 0, mu, 0, 0, 1, numEout ) ]
            MF6 += endfFormatsModule.endfInterpolationList( [ numEout, pdf_of_EpInterpolation ] )
            xys = entries.copyDataToXYs( )
            MF6 += endfFormatsModule.endfNdDataList( xys )
    return( 7, standardsModule.frames.labToken, MF6 )

angularEnergyModule.XYs3d.toENDF6 = toENDF6

#
# LLNLAngularEnergyForm
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :    # FIXME appears to be broken code

    angularForm = self.angularForm
    angularEnergyForm = self.angularEnergyForm

    energy_inInterpolation, energy_inFunctionInterpolation, energy_inInterpolationQualifier = angularEnergyForm.axes[0].interpolation.getInterpolationTokens( )
    muInterpolation, muFunctionInterpolation, muQualifier = angularEnergyForm.axes[1].interpolation.getInterpolationTokens( )
    energy_outInterpolation, probabilityInterpolation, energy_outQualifier = angularEnergyForm.axes[2].interpolation.getInterpolationTokens( )
    frame = angularEnergyForm.getProductFrame( )
    axes = pointwise.defaultAxes( energyInterpolation = energy_inInterpolation, energyFunctionInterpolation = energy_inFunctionInterpolation, 
            energyInterpolationQualifier = energy_inInterpolationQualifier, muInterpolation = muInterpolation, 
            energy_outInterpolation = energy_outInterpolation, probabilityInterpolation = probabilityInterpolation )
    E_inRatio = PQUModule.PQU( 1, angularEnergyForm.axes[0].getUnit( ) ).getValueAs( 'eV' )
    E_outRatio = PQUModule.PQU( 1, angularEnergyForm.axes[2].getUnit( ) ).getValueAs( 'eV' )
    LAW7 = pointwise( axes, self.getProductFrame( ) )

    if( len( angularForm ) != len( angularEnergyForm ) ) :
        raise Exception( "len( angularForm ) = %s != len( angularEnergyForm ) = %s" % ( len( angularForm ), len( angularEnergyForm ) ) )
    for indexE, EMuP in enumerate( angularForm ) :
        EMuEpP = angularEnergyForm[indexE]
        if( EMuP.value != EMuEpP.value ) : raise Exception( "At indexE = %d, EMuP.value %s != EMuEpP.value = %s" % ( indexE, EMuP.value, EMuEpP.value ) )
        if( len( EMuP ) != len( EMuEpP ) ) :
            raise Exception( "At indexE = %d (E_in = %s), len( EMuP ) %s != len( EMuEpP ) = %s" % ( indexE, EMuP.value, len( EMuP ), len( EMuEpP ) ) )
        w_xys = W_XYs.W_XYs( axesW_XY, index = indexE, value = EMuP.value )
        for indexMu, muP in enumerate( EMuP ) :
            muEpP = EMuEpP[indexMu]
            if( muP[0] != muEpP.value ) : raise Exception( "At indexE = %d, mu = %s != muEpP.value = %s" % ( indexE, muP[0], muEpP.value ) )
            xys = [ [ E_outRatio * Ep, muP[1] * P / E_outRatio ] for Ep, P in muEpP ]
            xys = XYs.XYs( axesXY, xys, accuracy = muEpP.getAccuracy( ), value = muP[0], index = indexMu, parent = w_xys )
            w_xys.append( XYs.XYs( axesXY, muEpP * muP[1], accuracy = muEpP.getAccuracy( ), value = muP[0], index = indexMu, parent = w_xys ) )
        LAW7.append( w_xys )

    LAW, frame, MF6 = LAW7.toENDF6( flags, targetInfo )
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )

angularEnergyModule.LLNLAngularEnergyForm.toENDF6 = toENDF6
