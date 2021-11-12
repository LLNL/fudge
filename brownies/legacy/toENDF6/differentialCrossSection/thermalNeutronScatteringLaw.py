# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import itertools

from fudge.reactionData.doubleDifferentialCrossSection import thermalNeutronScatteringLaw as thermalScatteringModule

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from brownies.legacy.toENDF6 import gndsToENDF6

#
# coherentElastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR = 1
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
    endf += self.S_table.toENDF6( flags, targetInfo, verbosityIndent = '' )
    endfMFList[7][2] = endf + [99999]
thermalScatteringModule.coherentElastic.form.toENDF6 = toENDF6

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
    independentInterp = gndsToENDF6.gndsToENDFInterpolationFlag( gridded.axes[1].interpolation )
    dependentInterp = gndsToENDF6.gndsToENDFInterpolationFlag( gridded.axes[2].interpolation )
    endf += ['%11i%11i%44s' % (len(data[0]), independentInterp, '' )]    # no trailing zeros
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *zip(Elist, data[0]) ) ) )

    # remaining temperatures:
    for T, datList in zip( Tlist[1:], data[1:] ):
        endf += [endfFormatsModule.endfHeadLine( T, 0, dependentInterp, 0, len(datList), 0 )]
        endf += endfFormatsModule.endfDataList( datList )
    return endf
thermalScatteringModule.coherentElastic.S_table.toENDF6 = toENDF6

#
# incoherentElastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR = 2
    if 7 in endfMFList and 2 in endfMFList[7]:
        # indicates that we already encountered coherentElastic data, so switch to LTHR=3:
        LTHR = 3
        endf = endfMFList[7][2][:-1]    # [:-1] drops SEND line from end of list
        endf[0] = endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )
    else:
        endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
    targetInfo['characteristicCrossSection'] = self.characteristicCrossSection.getValueAs('b')
    endf += self.DebyeWaller.toENDF6( flags, targetInfo, verbosityIndent = verbosityIndent )
    del targetInfo['characteristicCrossSection']
    endfMFList[7][2] = endf + [99999]
thermalScatteringModule.incoherentElastic.form.toENDF6 = toENDF6

#
# DebyeWaller
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
    data = self.function1d
    NR = 1; NP = len(data)
    endf = [endfFormatsModule.endfHeadLine( targetInfo['characteristicCrossSection'], 0, 0, 0, NR, NP )]
    interp = gndsToENDF6.gndsToENDFInterpolationFlag( data.interpolation )
    endf += endfFormatsModule.endfInterpolationList( (NP, interp) )
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *data.copyDataToXYs() ) ) )
    return endf
thermalScatteringModule.incoherentElastic.DebyeWaller.toENDF6 = toENDF6

#
# incoherentInelastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR, LAT, LASYM = 0, self.options.calculatedAtThermal, self.options.asymmetric
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, LAT, LASYM, 0 )]
    # describe scattering atoms:
    LLN, NS = 0, len( self.scatteringAtoms ) - 1
    NI = 6*(NS+1)
    endf += [endfFormatsModule.endfHeadLine( 0, 0, LLN, 0, NI, NS )]
    # principal scattering atom:
    atom = self.scatteringAtoms[0]
    neutronMass = self.findAttributeInAncestry('PoPs')['n'].mass.float('amu')
    endf += [endfFormatsModule.endfDataLine( [atom.freeAtomCrossSection.value * atom.numberPerMolecule,
        atom.e_critical.value, atom.mass.value / neutronMass, atom.e_max.value, 0, atom.numberPerMolecule] )]
    for atom in list(self.scatteringAtoms)[1:]:
        a1 = {'SCT':0.0, 'free_gas':1.0, 'diffusive_motion':2.0} [ atom.functionalForm ]
        endf += [endfFormatsModule.endfDataLine( [a1, atom.freeAtomCrossSection.value * atom.numberPerMolecule,
            atom.mass.value / neutronMass, 0, 0, atom.numberPerMolecule] )]

    # convert data form: sort first by beta, then E, then T
    gridded = self.S_alpha_beta.gridded3d
    array = gridded.array.constructArray()   # 3D numpy array

    Tlist = list( gridded.axes[3].values )
    betas = list( gridded.axes[2].values )
    alphas = list( gridded.axes[1].values )

    # switch array back to ENDF ordering:  1st beta, then T, then alpha:
    array = array.transpose( (1,0,2) )

    NR = 1; NB = len(betas)
    beta_interp = gndsToENDF6.gndsToENDFInterpolationFlag( gridded.axes[2].interpolation )
    endf += [endfFormatsModule.endfHeadLine( 0, 0, 0, 0, NR, NB )]
    endf += endfFormatsModule.endfInterpolationList( (NB, beta_interp) )
    #endf += ['%11i%11i%44s' % ( NB, 4, '' )]  # FIXME add 'suppressTrailingZeros' option to endfInterpolationList

    LT = len(Tlist)-1
    if LT:
        T_interp = gndsToENDF6.gndsToENDFInterpolationFlag( gridded.axes[3].interpolation )
    else:
        T_interp = None
    alpha_interp = gndsToENDF6.gndsToENDFInterpolationFlag( gridded.axes[1].interpolation )

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
thermalScatteringModule.incoherentInelastic.form.toENDF6 = toENDF6

#
# T_effective
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
    data = self.function1d
    NR = 1; NP = len(data)
    endf = [endfFormatsModule.endfHeadLine( 0, 0, 0, 0, NR, NP )]
    interp = gndsToENDF6.gndsToENDFInterpolationFlag( data.interpolation )
    endf += endfFormatsModule.endfInterpolationList( (NP, interp) )
    #endf += ['%11i%11i%44s' % (len(self), interp, '')]
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *data.copyDataToXYs() ) ) )
    return endf
thermalScatteringModule.incoherentInelastic.T_effective.toENDF6 = toENDF6
