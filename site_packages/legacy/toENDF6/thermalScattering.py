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

import itertools

from fudge.gnd import thermalScattering as thermalScatteringModule

from . import endfFormats as endfFormatsModule
from . import gndToENDF6

#
# thermalScattering
#
def toENDF6( self, flags = {}, verbosityIndent = '', covarianceSuite = None ):
    endfMFList = { 1 : { 451 : [] }, 7 : {} }
    targetInfo = {'ZA':self.MAT + 100, 'mass':float(self.mass)}
    MAT = self.MAT
    NSUB, NVER = 12, 8  # 12: thermal scattering sub-library, 8: ENDF/B-VIII

    for subsection in (self.coherentElastic, self.incoherentElastic, self.incoherentInelastic):
        if subsection is not None:
            subsection.toENDF6( endfMFList, flags, targetInfo, verbosityIndent )

    endfDoc = self.documentation.getLines()
    docHeader = [
            endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], -1, 0, 0, 0 ),
            endfFormatsModule.endfHeadLine( 0, 0, 0, 0, 0, 6 ),    # ENDF-6
            endfFormatsModule.endfHeadLine( 1.0, self.cutoffEnergy.value, 1, 0, NSUB, NVER ),
            endfFormatsModule.endfHeadLine( 0.0, 0.0, 0, 0, len(endfDoc), 0 ) ]
    endfMFList[1][451] = docHeader + endfDoc

    return endfFormatsModule.endfMFListToFinalFile( endfMFList, MAT, lineNumbers=True )
thermalScatteringModule.thermalScattering.toENDF6 = toENDF6

#
# coherentElastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR = 1
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
    endf += self.S_table.toENDF6( flags, targetInfo, verbosityIndent = '' )
    endfMFList[7][2] = endf + [99999]
thermalScatteringModule.coherentElastic.toENDF6 = toENDF6

#
# S_table
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):

    gridded = self.gridded2d
    data = gridded.array.constructArray()
    Tlist = list( gridded.axes[-1].values )   #FIXME what if grid units aren't 'K'?
    Elist = list( gridded.axes[-2].values )   #FIXME      ""                  'eV'?

    LT = len(Tlist)-1
    # first temperature includes the energy list:
    endf = [endfFormatsModule.endfHeadLine( Tlist[0], 0, LT, 0, 1, len(data[0]) ) ]
    independentInterp = gndToENDF6.gndToENDFInterpolationFlag( gridded.axes[1].interpolation )
    dependentInterp = gndToENDF6.gndToENDFInterpolationFlag( gridded.axes[2].interpolation )
    endf += ['%11i%11i%44s' % (len(data[0]), independentInterp, '' )]    # no trailing zeros
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *zip(Elist, data[0]) ) ) )

    # remaining temperatures:
    for T, datList in zip( Tlist[1:], data[1:] ):
        endf += [endfFormatsModule.endfHeadLine( T, 0, dependentInterp, 0, len(datList), 0 )]
        endf += endfFormatsModule.endfDataList( datList )
    return endf
thermalScatteringModule.S_table.toENDF6 = toENDF6

#
# incoherentElastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR = 2
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
    targetInfo['characteristicCrossSection'] = self.characteristicCrossSection.getValueAs('b')
    endf += self.DebyeWaller.toENDF6( flags, targetInfo, verbosityIndent = '' )
    del targetInfo['characteristicCrossSection']
    endfMFList[7][2] = endf + [99999]
thermalScatteringModule.incoherentElastic.toENDF6 = toENDF6

#
# DebyeWaller
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
    NR = 1; NP = len(self)
    endf = [endfFormatsModule.endfHeadLine( targetInfo['characteristicCrossSection'], 0, 0, 0, NR, NP )]
    interp = gndToENDF6.gndToENDFInterpolationFlag( self.interpolation )
    #endf += ['%11i%11i%44s' % (len(self), interp, '')]
    endf += endfFormatsModule.endfInterpolationList( (len(self), interp) )
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *self.copyDataToXYs() ) ) )
    return endf
thermalScatteringModule.DebyeWaller.toENDF6 = toENDF6

#
# incoherentInelastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR, LAT, LASYM = 0, self.calculatedAtThermal, self.asymmetric
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, LAT, LASYM, 0 )]
    # describe scattering atoms:
    LLN, NS = 0, len( self.scatteringAtoms ) - 1
    NI = 6*(NS+1)
    endf += [endfFormatsModule.endfHeadLine( 0, 0, LLN, 0, NI, NS )]
    # principal scattering atom:
    atom = self.scatteringAtoms[0]
    endf += [endfFormatsModule.endfDataLine( [atom.freeAtomCrossSection.value,
        atom.e_critical.value, atom.mass.value, atom.e_max.value, 0, atom.numberPerMolecule] )]
    for atom in list(self.scatteringAtoms)[1:]:
        a1 = {'SCT':0.0, 'free_gas':1.0, 'diffusive_motion':2.0} [ atom.functionalForm ]
        endf += [endfFormatsModule.endfDataLine( [a1, atom.freeAtomCrossSection.value,
            atom.mass.value, 0, 0, atom.numberPerMolecule] )]

    # convert data form: sort first by beta, then E, then T
    gridded = self.S_alpha_beta.gridded3d
    array = gridded.array.constructArray()   # 3D numpy array

    Tlist = list( gridded.axes[3].values )
    betas = list( gridded.axes[2].values )
    alphas = list( gridded.axes[1].values )

    # switch array back to ENDF ordering:  1st beta, then T, then alpha:
    array = array.transpose( (1,0,2) )

    NR = 1; NB = len(betas)
    endf += [endfFormatsModule.endfHeadLine( 0, 0, 0, 0, NR, NB )]
    endf += endfFormatsModule.endfInterpolationList( (NB, 4) )
    #endf += ['%11i%11i%44s' % ( NB, 4, '' )]  # FIXME add 'suppressTrailingZeros' option to endfInterpolationList

    LT = len(Tlist)-1
    if LT:
        T_interp = gndToENDF6.gndToENDFInterpolationFlag( gridded.axes[3].interpolation )
    else:
        T_interp = None
    beta_interp = gndToENDF6.gndToENDFInterpolationFlag( gridded.axes[2].interpolation )
    alpha_interp = gndToENDF6.gndToENDFInterpolationFlag( gridded.axes[1].interpolation )

    for index, beta in enumerate( betas ):
        data = array[index,:,:] # 2D sub-array for this beta

        endf += [endfFormatsModule.endfHeadLine( Tlist[0], beta, LT, 0, 1, len(data[0]) )]
        endf += endfFormatsModule.endfInterpolationList( (len(alphas), alpha_interp) )
        #endf += ['%11i%11i%44s' % (len(alphas), alpha_interp, '' )]    # no trailing zeros
        # For each beta, the first temperature needs to include the energy list:
        endf += endfFormatsModule.endfDataList( list( itertools.chain( *zip(alphas, data[0]) ) ) )

        # remaining temperatures:
        for T, datList in zip( Tlist[1:], data[1:] ):
            endf += [endfFormatsModule.endfHeadLine( T, beta, T_interp, 0, len(datList), 0 )]
            endf += endfFormatsModule.endfDataList( datList )

    for atom in self.scatteringAtoms:
        if atom.T_effective is not None:
            endf += atom.T_effective.toENDF6( flags, targetInfo, verbosityIndent = '' )
    endfMFList[7][4] = endf + [99999]
thermalScatteringModule.incoherentInelastic.toENDF6 = toENDF6

#
# T_effective
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
    NR = 1; NP = len(self)
    endf = [endfFormatsModule.endfHeadLine( 0, 0, 0, 0, NR, NP )]
    interp = gndToENDF6.gndToENDFInterpolationFlag( self.interpolation )
    endf += endfFormatsModule.endfInterpolationList( (len(self), interp) )
    #endf += ['%11i%11i%44s' % (len(self), interp, '')]
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *self.copyDataToXYs() ) ) )
    return endf
thermalScatteringModule.T_effective.toENDF6 = toENDF6
