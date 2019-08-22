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

from fudge.core.utilities import brb

import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module

import fudge.gnd.productData.distributions.angular as angularModule
import xData.standards as standardsModule
import xData.XYs as XYsModule
import xData.multiD_XYs as multiD_XYsModule
import xData.regions as regionsModule
import xData.series1d as series1dModule

import pqu.PQU as PQUModule

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    if( MT==455 ) and ( MT in endfMFList[4] ):
        return  # only write one 'isotropic' section for all delayed neutron groups
    angularSubform = self.angularSubform
    if( hasattr( angularSubform, 'toENDF6' ) ) :

        doMF4AsMF6 = targetInfo['doMF4AsMF6']
        MF = 4
        NM = 0
        frame = self.productFrame
        if( frame is None ) : frame = self.ancestor.productFrame  # Happens for uncorrelated distribution.
        if( isinstance( angularSubform, multiD_XYsModule.XYs2d ) ) :
            if( isinstance( angularSubform[0], XYsModule.XYs1d ) ) :
                LI, LTT, MF4 = toAngularPointwise( angularSubform, targetInfo, not( doMF4AsMF6 ) )
            elif( isinstance( angularSubform[0], series1dModule.LegendreSeries ) ) :
                interpolation, numberOfPoints, LI, LTT, NM, MF4Sub = toAngularLegendre( angularSubform, targetInfo, not( doMF4AsMF6 ) )
                MF4 = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, numberOfPoints ) ] + \
                        endfFormatsModule.endfInterpolationList( [ len( angularSubform ), interpolation ] )
                MF4 += MF4Sub
            else :
                raise 'hell - fix me'
        elif( isinstance( angularSubform, regionsModule.regions2d ) ) :
            LTT, MF4, numberOfPoints, LegendreInterpolations, LegendreData  = None, [], 0, [], []
            for ridx, region in enumerate(angularSubform) :
                targetInfo['skipFirstEnergy'] = False
                if ridx > 0 and type(angularSubform[ridx][0]) == type(angularSubform[ridx-1][-1]):
                    # FIXME: following should work once __eq__ fixed for xData. Shows up in ENDF-VII.1 Cu evaluations
                    #if angularSubform[ridx-1][-1] == angularSubform[ridx][0]: targetInfo['skipFirstEnergy'] = True
                    if isinstance(angularSubform[ridx][0], series1dModule.LegendreSeries):
                        if angularSubform[ridx][0].coefficients == angularSubform[ridx-1][-1].coefficients:
                            targetInfo['skipFirstEnergy'] = True
                    else:
                        raise NotImplementedError
                if( isinstance( region, angularModule.XYs2d ) ) :
                    if( isinstance( region[0], series1dModule.LegendreSeries ) ) :
                        interpolation, numberOfPointsSub, LI, LTTSub, NMtmp, MF4Sub = toAngularLegendre( region, targetInfo, not( doMF4AsMF6 ) )
                        numberOfPoints += numberOfPointsSub
                        NM = max(NM,NMtmp)
                        LegendreInterpolations += [ numberOfPoints, interpolation ]
                    elif( isinstance( region[0], XYsModule.XYs1d ) ) :
                        LI, LTTSub, MF4Sub = toAngularPointwise( region, targetInfo, not( doMF4AsMF6 ) )
                    else :
                        raise 'hell - fix me'
                    if( len( MF4 ) > 0 ) : MF4.pop( -1 )
                    MF4 += MF4Sub
                    if( LTT is None ) : LTT = LTTSub
                    if( LTT != LTTSub ) :
                        if( ( LTT == 1 ) and ( LTTSub == 2 ) ) :
                            LTT = 3
                        else :
                            raise 'hell - fix me'
                else :
                    raise 'hell - fix me'
                del targetInfo['skipFirstEnergy']
            if( len( LegendreInterpolations ) > 0 ) :
                MF4 = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, len( LegendreInterpolations ) / 2, LegendreInterpolations[-2] ) ] + \
                        endfFormatsModule.endfInterpolationList( LegendreInterpolations ) + MF4
            ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, len( angularSubform ) ) ]
            ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( angularSubform ), interpolation ] )
            LI = 0
        elif( isinstance( angularSubform, angularModule.recoil ) ) :
            if not targetInfo['doMF4AsMF6']:
                return      # recoil partners only get written to file 6
            LI, LTT, MF4 = angularSubform.toENDF6( flags, targetInfo )
            MF = 6
        elif( isinstance( angularSubform, angularModule.isotropic ) ) :
            LI, LTT, MF4 = angularSubform.toENDF6( flags, targetInfo )
            MF = 4
            if( doMF4AsMF6 ) : MF = 6
        else :
            brb.objectoutline( angularSubform )
            raise 'hell - fix me'
        if( doMF4AsMF6 ) :
            if( LTT in [ 0 ] ) :
                LAW = 3
            elif( LTT in [ 1, 2 ] ) :
                LAW = 2
            elif( LTT in [ 4 ] ) :
                LAW = 4
            else :
                raise Exception( 'LTT = %s needs a LAW' % LTT )
            gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF4 )
        else :
            LCT = { standardsModule.frames.labToken : 1, standardsModule.frames.centerOfMassToken : 2 }[frame]
            if( MT not in endfMFList[MF] ) : endfMFList[MF][MT] = []
            if LTT != 3: NM = 0
            endfMFList[MF][MT] += [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, LTT, 0, 0 ),
                                  endfFormatsModule.endfHeadLine(                0, targetInfo['mass'], LI, LCT, 0, NM ) ] + MF4
    else :
        print 'WARNING: subform %s does not have method toENDF6 for form %s' % ( targetInfo['style'], self.moniker )

angularModule.form.toENDF6 = toENDF6
angularModule.twoBodyForm.toENDF6 = toENDF6

def toAngularPointwise( angularSubform, targetInfo, insertSENDL ) :

    energyConversionFactor = PQUModule.PQU(1, angularSubform.axes[-1].unit ).getValueAs('eV')

    def angularPointwiseEnergy2ENDF6( self, targetInfo ) :

        interpolation = gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation )
        energy_in_eV = self.value * energyConversionFactor
        if( targetInfo['doMF4AsMF6'] ) :
            ENDFDataList = [ endfFormatsModule.endfContLine( 0, energy_in_eV, interpolation + 10, 0, 2 * len( self ), len( self ) ) ]
        else :
            ENDFDataList = [ endfFormatsModule.endfContLine( 0, energy_in_eV, 0, 0, 1, len( self ) ) ]
            ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( self ), interpolation ] )
        ENDFDataList += endfFormatsModule.endfNdDataList( self )
        return( ENDFDataList )

    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, len( angularSubform ) ) ]
    interpolation = ( gndToENDF6Module.gndToENDFInterpolationFlag( angularSubform.interpolation ) )
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( angularSubform ), interpolation ] )
    start = 0
    if targetInfo.get('skipFirstEnergy'): start = 1
    for energy_in in angularSubform[start:] : ENDFDataList += angularPointwiseEnergy2ENDF6( energy_in, targetInfo )
    if( insertSENDL ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( 0, 2, ENDFDataList )

def toAngularLegendre( angularSubform, targetInfo, insertSENDL ) :
    """This should only be called from this module."""

    NM = 0
    interpolation = gndToENDF6Module.gndToENDFInterpolationFlag( angularSubform.interpolation )
    energyConversionFactor = PQUModule.PQU(1, angularSubform.axes[-1].unit ).getValueAs('eV')
    ENDFDataList = []
    start = 0
    if targetInfo.get('skipFirstEnergy'): start = 1
    for energy in angularSubform[start:] :
        NW, NL = len( energy ) - 1, 0
        if( targetInfo['doMF4AsMF6'] ) : NL = NW
        ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.value * energyConversionFactor, 0, 0, NW, NL ) )
        ENDFDataList += endfFormatsModule.endfDataList( energy.coefficients[1:] )
        NM = max(NM, len(energy.coefficients[1:]))
    if( insertSENDL ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( interpolation, len( angularSubform[start:] ), 0, 1, NM, ENDFDataList )

#
# CoulombExpansionForm
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
    target, projectile = targetInfo['reactionSuite'].target, targetInfo['reactionSuite'].projectile
    LIDP = target==projectile
    if( isinstance( self.nuclear_term, angularModule.XYs2d ) ) :
        for ridx in xrange( len( self.nuclear_term ) ) :
            counts += 1
            nuclear, interferenceReal, interferenceImaginary = self.nuclear_term[ridx], self.interferenceReal_term[ridx], self.interferenceImaginary_term[ridx]
            LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal, interferenceImaginary, lineData )
        interpolationFlagsList += [ counts, gndToENDF6Module.gndToENDFInterpolationFlag( self.nuclear_term.interpolation ) ]
    elif( isinstance( self.nuclear_term, angularModule.regions2d ) ) :
        for regionIndex, region in enumerate( self.nuclear_term ) :
            interferenceReal, interferenceImaginary = self.interferenceReal_term[regionIndex], self.interferenceImaginary_term[regionIndex]
            for energyIndex, nuclear in enumerate( region ) :
                if( ( regionIndex != 0 ) and ( energyIndex == 0 ) ) : continue
                counts += 1
                LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal[energyIndex], interferenceImaginary[energyIndex], lineData )
            interpolationFlagsList += [ counts, gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation ) ]
    else :
        raise NotImplementedError( "Unknown data storage inside CoulombExpansion: %s" % type(self.nuclear_term) )
    interpolationFlags = endfFormatsModule.endfInterpolationList( interpolationFlagsList )
    ENDFDataList = [ endfFormatsModule.endfContLine( projectile.getSpin().value, 0, LIDP, 0, len( interpolationFlagsList ) / 2, counts ) ] + interpolationFlags + lineData
    if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    LAW = 5
    gndToENDF6Module.toENDF6_MF6(MT, endfMFList, flags, targetInfo, LAW, self.productFrame, ENDFDataList)

angularModule.CoulombExpansionForm.toENDF6 = toENDF6

#
# NuclearPlusCoulombInterference
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    LI, LTT, MF6 = self.angularSubform.toENDF6( flags, { 'doMF4AsMF6' : True } )
    target, projectile = targetInfo['reactionSuite'].target, targetInfo['reactionSuite'].projectile
    spin = projectile.getSpin().value
    LIDP = target==projectile
    NR, NE = [ int(a) for a in MF6[0].split()[-2:] ]
    MF6[0] = endfFormatsModule.endfContLine( spin, 0, LIDP, 0, NR, NE )
    LAW = 5
    gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, self.productFrame, MF6 )

angularModule.NuclearPlusCoulombInterferenceForm.toENDF6 = toENDF6

#
# isotropic
#
def toENDF6( self, flags, targetInfo ) :

    ENDFDataList = []
    if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( 1, 0, ENDFDataList )

angularModule.isotropic.toENDF6 = toENDF6

#
# forward
#
def toENDF6( self, flags, targetInfo ) :

    return( 1, 0, [] )

angularModule.forward.toENDF6 = toENDF6

#
# recoil
#
def toENDF6( self, flags, targetInfo ) :

    return( 1, 4, [] )

angularModule.recoil.toENDF6 = toENDF6

#
#XYs2d 
#
def toENDF6( self, flags, targetInfo ) :

    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
    interpolation = gndToENDF6Module.gndToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( self ), interpolation ] )
    for energy_in in self : 
        if( isinstance( energy_in, XYsModule.XYs1d ) ) :
            ENDFDataList += gndToENDF6Module.angularPointwiseEnergy2ENDF6( energy_in, targetInfo )
        elif( isinstance( energy_in, series1dModule.LegendreSeries ) ) :
            ENDFDataList += gndToENDF6Module.angularLegendreEnergy2ENDF6( energy_in, targetInfo )
        else :
            raise 'hell - fix me'
    if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( 0, 2, ENDFDataList )

angularModule.XYs2d.toENDF6 = toENDF6

#
# regions2d
#
def toENDF6( self, flags, targetInfo ) :

    raise 'hell - not implemented'

angularModule.regions2d.toENDF6 = toENDF6
