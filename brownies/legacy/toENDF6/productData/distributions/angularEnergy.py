# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule
from xData import standards as standardsModule

from fudge.productData.distributions import angularEnergy as angularEnergyModule

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
        print( 'WARNING: angularEnergy subform "%s" has no toENDF6 method' % subform.moniker )

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
        Ein = oneEin.outerDomainValue * energyConversionFactor
        numMu = len( oneEin )
        MF6 += [ endfFormatsModule.endfContLine( 0, Ein, 0, 0, 1, numMu ) ]
        MF6 += endfFormatsModule.endfInterpolationList( [ numMu, muInterpolation ] )
        for entries in oneEin :
            pdf_of_EpInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( entries.interpolation )
            mu = entries.outerDomainValue
            numEout = len( entries )
            MF6 += [ endfFormatsModule.endfContLine( 0, mu, 0, 0, 1, numEout ) ]
            MF6 += endfFormatsModule.endfInterpolationList( [ numEout, pdf_of_EpInterpolation ] )
            xys = entries.copyDataToXYs( )
            MF6 += endfFormatsModule.endfNdDataList( xys )
    return( 7, standardsModule.frames.labToken, MF6 )

angularEnergyModule.XYs3d.toENDF6 = toENDF6
