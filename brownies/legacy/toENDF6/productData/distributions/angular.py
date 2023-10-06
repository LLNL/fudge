# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule
from xData import series1d as series1dModule

from fudge.core.utilities import brb

from fudge.productData.distributions import angular as angularModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

#
# form
#
def toENDF6(self, MT, endfMFList, flags, targetInfo):

    if MT == 455: return                # No angular data for delayed neutrons in ENDF
    angularSubform = self.angularSubform
    if hasattr(angularSubform, 'toENDF6'):

        doMF4AsMF6 = targetInfo['doMF4AsMF6']
        MF = 4
        NM = 0
        frame = self.productFrame
        if frame is None: frame = self.ancestor.productFrame  # Happens for uncorrelated distribution.
        if isinstance(angularSubform, multiD_XYsModule.XYs2d):
            if isinstance(angularSubform[0], XYs1dModule.XYs1d):
                interpolation, numberOfPoints, LI, LTT, MF4Sub = toAngularPointwise(angularSubform, targetInfo)
                MF4 = [endfFormatsModule.endfContLine(0, 0, 0, 0, 1, numberOfPoints)] + \
                       endfFormatsModule.endfInterpolationList([len(angularSubform), interpolation])
                MF4 += MF4Sub
            elif isinstance(angularSubform[0], series1dModule.LegendreSeries):
                interpolation, numberOfPoints, LI, LTT, NM, MF4Sub = toAngularLegendre(angularSubform, targetInfo)
                MF4 = [endfFormatsModule.endfContLine(0, 0, 0, 0, 1, numberOfPoints)] + \
                       endfFormatsModule.endfInterpolationList([len(angularSubform), interpolation])
                MF4 += MF4Sub
            else:
                raise 'hell - fix me'
            if not doMF4AsMF6:
                MF4.append(endfFormatsModule.endfSENDLineNumber())
        elif isinstance(angularSubform, regionsModule.Regions2d):
            LTT = None
            NLegendre, LegendreInterpolations, LegendreData = 0, [], []
            NXY, XYInterpolations, XYData = 0, [], []
            for ridx, region in enumerate(angularSubform):
                targetInfo['skipFirstEnergy'] = False
                if ridx > 0 and type(angularSubform[ridx][0]) is type(angularSubform[ridx-1][-1]):
                    if angularSubform[ridx-1][-1] == angularSubform[ridx][0]: targetInfo['skipFirstEnergy'] = True
                if isinstance(region, angularModule.XYs2d):
                    if isinstance(region[0], series1dModule.LegendreSeries):
                        interpolation, numberOfPointsSub, LI, LTTSub, NMtmp, MF4Sub = toAngularLegendre(region, targetInfo)
                        NLegendre += numberOfPointsSub
                        NM = max(NM,NMtmp)
                        LegendreInterpolations += [NLegendre, interpolation]
                        LegendreData += MF4Sub
                    elif isinstance(region[0], XYs1dModule.XYs1d):
                        interpolation, numberOfPointsSub, LI, LTTSub, MF4Sub = toAngularPointwise(region, targetInfo)
                        NXY += numberOfPointsSub
                        XYInterpolations += [NXY, interpolation]
                        XYData += MF4Sub
                    else:
                        raise 'hell - fix me'
                    if LTT is None: LTT = LTTSub
                    if LTT != LTTSub:
                        if LTT in (1, 3) and LTTSub == 2:
                            LTT = 3
                        else:
                            raise 'hell - fix me'
                else:
                    raise 'hell - fix me'
                del targetInfo['skipFirstEnergy']
            MF4 = []
            if len(LegendreInterpolations) > 0:
                MF4 = ([endfFormatsModule.endfContLine(
                    0, 0, 0, 0, len(LegendreInterpolations) / 2, LegendreInterpolations[-2])] +
                       endfFormatsModule.endfInterpolationList(LegendreInterpolations) + LegendreData)
            if len(XYInterpolations) > 0:
                MF4 += ([endfFormatsModule.endfContLine(
                    0, 0, 0, 0, len(XYInterpolations) / 2, XYInterpolations[-2])] +
                       endfFormatsModule.endfInterpolationList(XYInterpolations) + XYData)
            LI = 0
            if not doMF4AsMF6:
                MF4.append(endfFormatsModule.endfSENDLineNumber())
        elif isinstance(angularSubform, angularModule.Recoil):
            if not targetInfo['doMF4AsMF6']:
                return      # recoil partners only get written to file 6
            LI, LTT, MF4 = angularSubform.toENDF6(flags, targetInfo)
            MF = 6
        elif isinstance(angularSubform, angularModule.Isotropic2d):
            LI, LTT, MF4 = angularSubform.toENDF6(flags, targetInfo)
            MF = 4
            if doMF4AsMF6: MF = 6
        else:
            brb.objectoutline(angularSubform)
            raise 'hell - fix me'
        if doMF4AsMF6:
            if LTT in [0]:
                LAW = 3
            elif LTT in [1, 2]:
                LAW = 2
            elif LTT in [4]:
                LAW = 4
            else:
                raise Exception('LTT = %s needs a LAW' % LTT)
            gndsToENDF6Module.toENDF6_MF6(MT, endfMFList, flags, targetInfo, LAW, frame, MF4)
        else:
            LCT = {xDataEnumsModule.Frame.lab: 1, xDataEnumsModule.Frame.centerOfMass: 2}[frame]
            if MT not in endfMFList[MF]: endfMFList[MF][MT] = []
            if LTT != 3: NM = 0
            endfMFList[MF][MT] += [endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], 0, LTT, 0, 0),
                                   endfFormatsModule.endfHeadLine(0, targetInfo['mass'], LI, LCT, 0, NM)] + MF4
    else:
        print('WARNING: subform %s does not have method toENDF6 for form %s' % (targetInfo['style'], self.moniker))

angularModule.Form.toENDF6 = toENDF6
angularModule.TwoBody.toENDF6 = toENDF6

def toAngularPointwise(angularSubform, targetInfo):

    energyConversionFactor = PQUModule.PQU(1, angularSubform.axes[-1].unit ).getValueAs('eV')

    def angularPointwiseEnergy2ENDF6(self, targetInfo):

        interpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(self.interpolation)
        energy_in_eV = self.outerDomainValue * energyConversionFactor
        if targetInfo['doMF4AsMF6']:
            ENDFDataList = [endfFormatsModule.endfContLine(
                0, energy_in_eV, interpolation + 10, 0, 2 * len(self), len(self))]
        else:
            ENDFDataList = [endfFormatsModule.endfContLine(0, energy_in_eV, 0, 0, 1, len(self))]
            ENDFDataList += endfFormatsModule.endfInterpolationList([len(self), interpolation])
        ENDFDataList += endfFormatsModule.endfNdDataList(self)
        return ENDFDataList

    ENDFDataList = []
    start = 0
    interpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(angularSubform.interpolation)
    if targetInfo.get('skipFirstEnergy'): start = 1
    for energy_in in angularSubform[start:]: ENDFDataList += angularPointwiseEnergy2ENDF6(energy_in, targetInfo)
    return interpolation, len(angularSubform[start:]), 0, 2, ENDFDataList

def toAngularLegendre(angularSubform, targetInfo):
    """This should only be called from this module."""

    NM = 0
    interpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(angularSubform.interpolation)
    energyConversionFactor = PQUModule.PQU(1, angularSubform.axes[-1].unit).getValueAs('eV')
    ENDFDataList = []
    start = 0
    if targetInfo.get('skipFirstEnergy'): start = 1
    for energy in angularSubform[start:]:
        NW, NL = len(energy) - 1, 0
        if targetInfo['doMF4AsMF6']: NL = NW
        ENDFDataList.append(endfFormatsModule.endfContLine(
            0, energy.outerDomainValue * energyConversionFactor, 0, 0, NW, NL))
        ENDFDataList += endfFormatsModule.endfDataList(energy.coefficients[1:])
        NM = max(NM, len(energy.coefficients[1:]))
    return interpolation, len(angularSubform[start:]), 0, 1, NM, ENDFDataList

#
# isotropic
#
def toENDF6( self, flags, targetInfo ) :

    ENDFDataList = []
    if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( 1, 0, ENDFDataList )

angularModule.Isotropic2d.toENDF6 = toENDF6

#
# Forward
#
def toENDF6( self, flags, targetInfo ) :

    return( 1, 0, [] )

angularModule.Forward.toENDF6 = toENDF6

#
# Recoil
#
def toENDF6( self, flags, targetInfo ) :

    return( 1, 4, [] )

angularModule.Recoil.toENDF6 = toENDF6

#
#XYs2d 
#
def toENDF6( self, flags, targetInfo ) :

    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
    interpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( self ), interpolation ] )
    EInFactor = PQUModule.PQU( 1, self.axes[-1].unit ).getValueAs( 'eV' )
    for energy_in in self : 
        if( isinstance( energy_in, XYs1dModule.XYs1d ) ) :
            ENDFDataList += gndsToENDF6Module.angularPointwiseEnergy2ENDF6( energy_in, targetInfo, EInFactor )
        elif( isinstance( energy_in, series1dModule.LegendreSeries ) ) :
            ENDFDataList += gndsToENDF6Module.angularLegendreEnergy2ENDF6( energy_in, targetInfo, EInFactor )
        else :
            raise 'hell - fix me'
    if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormatsModule.endfSENDLineNumber( ) )
    return( 0, 2, ENDFDataList )

angularModule.XYs2d.toENDF6 = toENDF6

#
# Regions2d
#
def toENDF6( self, flags, targetInfo ) :

    raise 'hell - not implemented'

angularModule.Regions2d.toENDF6 = toENDF6
