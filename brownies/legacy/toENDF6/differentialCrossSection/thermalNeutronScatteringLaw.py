# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import itertools
import numpy

from PoPs import IDs as IDsPoPsModule

from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentInelastic

from xData import gridded as griddedModule

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
thermalScatteringModule.coherentElastic.Form.toENDF6 = toENDF6

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
    targetInfo['boundAtomCrossSection'] = self.boundAtomCrossSection.getValueAs('b')
    endf += self.DebyeWallerIntegral.toENDF6( flags, targetInfo, verbosityIndent = verbosityIndent )
    del targetInfo['boundAtomCrossSection']
    endfMFList[7][2] = endf + [99999]
thermalScatteringModule.incoherentElastic.Form.toENDF6 = toENDF6

#
# DebyeWallerIntegral
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
    data = self.function1d
    NR = 1; NP = len(data)
    endf = [endfFormatsModule.endfHeadLine( targetInfo['boundAtomCrossSection'], 0, 0, 0, NR, NP )]
    interp = gndsToENDF6.gndsToENDFInterpolationFlag( data.interpolation )
    endf += endfFormatsModule.endfInterpolationList( (NP, interp) )
    endf += endfFormatsModule.endfDataList( list( itertools.chain( *data.copyDataToXYs() ) ) )
    return endf
thermalScatteringModule.incoherentElastic.DebyeWallerIntegral.toENDF6 = toENDF6

#
# incoherentInelastic
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    LTHR, LAT = 0, self.calculatedAtThermal
    neutronMass = self.findAttributeInAncestry('PoPs')[IDsPoPsModule.neutron].mass.float('amu')

    scatteringKernel = self.getPrimaryScatterer().selfScatteringKernel.kernel
    if not isinstance(scatteringKernel, griddedModule.Gridded3d):
        raise TypeError("ENDF-6 does not support %s scattering kernel for primary scatterer" % type(scatteringKernel))

    Tlist = list( scatteringKernel.axes[3].values )
    betas = list( scatteringKernel.axes[2].values )
    alphas = list( scatteringKernel.axes[1].values )

    LASYM = min(betas) < 0  # S_ab not symmetric in beta
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, LTHR, LAT, LASYM, 0 )]

    # describe scattering atoms:
    LLN, NS = 0, len( self.scatteringAtoms ) - 1
    if self in targetInfo['ENDFconversionFlags']:
        flags = targetInfo['ENDFconversionFlags'][self]
        LLN, LLN_min = flags.split(',')
        LLN = int(LLN.split('=')[1])
        LLN_min = float(LLN_min.split('=')[1])
    NI = 6*(NS+1)
    endf += [endfFormatsModule.endfHeadLine( 0, 0, LLN, 0, NI, NS )]

    # principal scattering atom:
    atom1 = self.getPrimaryScatterer()
    freeAtomCrossSection = atom1.boundAtomCrossSection.value * atom1.mass.value**2 / (atom1.mass.value + neutronMass)
    endf += [endfFormatsModule.endfDataLine( [freeAtomCrossSection * atom1.numberPerMolecule,
        atom1.e_critical.value, atom1.mass.value / neutronMass, atom1.e_max.value, 0, atom1.numberPerMolecule] )]
    for atom in list(self.scatteringAtoms):
        if atom is atom1: continue
        a1 = {
            incoherentInelastic.SCTApproximation: 0.0,
            incoherentInelastic.FreeGasApproximation: 1.0
            # note: 'diffusive motion' (a1 = 2.0) appears to be unused and isn't supported by GNDS-2.0
        }[type(atom.selfScatteringKernel.kernel)]
        freeAtomCrossSection = atom.boundAtomCrossSection.value * atom.mass.value**2 / (atom.mass.value + neutronMass)
        endf += [endfFormatsModule.endfDataLine( [a1, freeAtomCrossSection * atom.numberPerMolecule,
            atom.mass.value / neutronMass, 0, 0, atom.numberPerMolecule] )]

    # convert data form: sort first by beta, then E, then T
    array = scatteringKernel.array.constructArray()   # 3D numpy array

    # switch array back to ENDF ordering:  1st beta, then T, then alpha:
    array = array.transpose( (1,0,2) )

    interpolations = {
        'T': scatteringKernel.axes[3].interpolation,
        'beta': scatteringKernel.axes[2].interpolation,
        'alpha': scatteringKernel.axes[1].interpolation
    }

    if LLN == 1:
        with numpy.errstate(invalid='ignore', divide='ignore'):
            array2 = numpy.log(array)
            array2[array==0] = LLN_min
            array = array2

        for key in interpolations:
            if interpolations[key] is None:
                continue
            if not interpolations[key].startswith('log'):
                raise ValueError("Expected log interpolation on S (linear on ln(S)) but got %s" % interpolations[key])
            interpolations[key] = interpolations[key].replace('log','lin',1)

    NR = 1; NB = len(betas)
    beta_interp = gndsToENDF6.gndsToENDFInterpolationFlag( interpolations['beta'] )
    endf += [endfFormatsModule.endfHeadLine( 0, 0, 0, 0, NR, NB )]
    endf += endfFormatsModule.endfInterpolationList( (NB, beta_interp) )
    #endf += ['%11i%11i%44s' % ( NB, 4, '' )]  # FIXME add 'suppressTrailingZeros' option to endfInterpolationList

    LT = len(Tlist)-1
    if LT:
        T_interp = gndsToENDF6.gndsToENDFInterpolationFlag( interpolations['T'] )
    else:
        T_interp = None
    alpha_interp = gndsToENDF6.gndsToENDFInterpolationFlag( interpolations['alpha'] )

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
thermalScatteringModule.incoherentInelastic.Form.toENDF6 = toENDF6

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
