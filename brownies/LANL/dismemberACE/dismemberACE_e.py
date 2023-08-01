# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module dismembers an continuous-energy electron ACE data file (identifier "e" is for continuous-energy electron tables).
'''

# Description of format per http://resource.npl.co.uk/docs/science_technology/ionising%20radiation/clubs_groups/monte_carlo/2004/workshop/electron_upgrade_paper.pdf.

from brownies.LANL.dismemberACE import dismemberACE_misc as dismemberACE_miscModule

def dismemberACE( args, NXS, JXS, XSS ) :

    NXS1  = NXS[1]                  # Length of XSS data.
    NXS2  = NXS[2]                  # Z
    NXS3  = NXS[3]                  # Number of radiation stopping power points.
    NXS4  = NXS[4]                  # Number of Mott scattering cross section corrections.
    NXS5  = NXS[5]                  # Number of electron energy points for bremsstrahlung.
    NXS6  = NXS[6]                  # Number of photon ratio points for bremsstrahlung.
# The following are only used for the "el3w" format for which NXS[16] is 3.
    NXS9  = NXS[9]                  # Number of points in the bremsstrahlung database for bremsstrahlung spectrum calculation.
    NXS10 = NXS[10]                 # Number of points in the bremsstrahlung database for angular/energy calculation.
    NXS11 = NXS[11]                 # Number of oscillator points for density effect calculation.
    NXS16 = NXS[16]                 # 3 for the "el3w" format and unused for the "el1" format.

    JXS1  = JXS[1]                  # Probably always 1 (Location of the line data).
    JXS2  = JXS[2]                  # Location of the radiation stopping power data.
    JXS3  = JXS[3]                  # Location of Mott scattering cross section.
    JXS4  = JXS[4]                  # Location of Riley cross section.
    JXS5  = JXS[5]                  # Location of ITS3.0 bremsstrahlung production data (ITS1.0 for the "el1" format and ITS3.0 for the "el3w" format).
    JXS7  = JXS[7]                  # Location of internally calcualed Riley cross section.
# The following are only used for the "el3w" format for which NXS[16] is 3.
    JXS8  = JXS[8]                  # Location of internally calcualed functions of Z.
    JXS9  = JXS[9]                  # Location of photon energy ratio.
    JXS10 = JXS[10]                 # Location of photon energy ratio for angular data.
    JXS11 = JXS[11]                 # Location of oscillator descriptions for density effect calcualion.

# Only used in the "el1" format.
    JXS6  = JXS[6]                  # Location of bremsstrahlung production data.

    numberInInfo = 2 if NXS16 == 3 else 4
    info = dismemberACE_miscModule.getData( 'Stopping power energies', XSS, JXS1, numberInInfo )
    EDG = info[0]                   # K edge below which no electron induced relaxation will occur.
    EEK = info[1]                   # el1 format: K x-ray of=r Auger electron emission energy.
                                    # el3w format: Auger electron emission energy (E_k - 2 * E_l).
    fOut = dismemberACE_miscModule.openFile( 0, 'Info' )
    if NXS16 == 3:
        fOut.write('# Format is %s.\n' % 'el3w')
        fOut.write(' %20.12e # K edge below which no electron induced relaxation will occur.\n' % EDG)
        fOut.write(' %20.12e # Auger electron emission energy (E_k - 2 * E_l).' % EEK)
    else :
        fOut.write('# Formant is %s.\n' % 'el1')
        fOut.write(' %20.12e # K edge below which no electron induced relaxation will occur.\n' % EDG)
        fOut.write(' %20.12e # K x-ray or Auger electron emission energy.\n' % EEK)
        fOut.write(' %20.12e\n' % info[2])
        fOut.write(' %20.12e\n' % info[3])
    fOut.close()

    fOut = dismemberACE_miscModule.openFile( 0, 'radiationStoppingPower' )
    fOut.write( '#        energy             stop power              nu (for the el1 format nu is all 0.0s).\n' )
    energies = dismemberACE_miscModule.getData( 'Stopping power energies', XSS, JXS2, NXS3 )
    stopPowers = dismemberACE_miscModule.getData( 'Stopping power', XSS, JXS2 + NXS3, NXS3 )
    nus = NXS3 * [ 0.0 ]
    if NXS16 == 3: nus = dismemberACE_miscModule.getData( 'Stopping power nus', XSS, JXS2 + 2 * NXS3, NXS3 )
    for index in range( NXS3 ) :
        fOut.write( ' %20.12e %20.12e %20.12e\n' % ( energies[index], stopPowers[index], nus[index] ) )
    
    fOut = dismemberACE_miscModule.openFile( 0, 'Mott' )
    Mott = dismemberACE_miscModule.getData( 'Stopping power energies', XSS, JXS3, 6 * NXS4 )
    fOut.write( '#        energy               h(0)                h(pi/4)              h(pi/2)            h(3pi/4)                h(pi)\n' )
    for index in range(NXS4):
        fOut.write( ' %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n' % ( Mott[index], Mott[index+NXS4], Mott[index+2*NXS4], Mott[index+3*NXS4], Mott[index+4*NXS4], Mott[index+5*NXS4] ) )

    fOut = dismemberACE_miscModule.openFile( 0, 'RileyScatteringCrossSection' )
    N9 = 9
    N14 = 14
    RileyScattering = dismemberACE_miscModule.getData( 'Riley scattering cross section parameters.', XSS, JXS4, N9 * N14 )
    fOut.write('# %d rows by %d columns\n' % (N9, N14))
    for index in range(N9):
        fOut.write( ' %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n' %
                (RileyScattering[index*N14], RileyScattering[index*N14+1], RileyScattering[index*N14+2], RileyScattering[index*N14+3], 
                RileyScattering[index*N14+4], RileyScattering[index*N14+5], RileyScattering[index*N14+6], RileyScattering[index*N14+7],
                RileyScattering[index*N14+8], RileyScattering[index*N14+9], RileyScattering[index*N14+10], RileyScattering[index*N14+11],
                RileyScattering[index*N14+12], RileyScattering[index*N14+13]))

    fOut = dismemberACE_miscModule.openFile( 0, 'BremsstrahlungCrossSection' )
    BremsstrahlungCrossSection = dismemberACE_miscModule.getData( 'Bremsstrahlung cross section', XSS, JXS5, NXS5 )
    for value in BremsstrahlungCrossSection: fOut.write( ' %20.12e\n' % value)

    fOut = dismemberACE_miscModule.openFile( 0, 'BremsstrahlungRatios' )
    if NXS16 == 3:
        BremsstrahlungRatios = dismemberACE_miscModule.getData( 'Bremsstrahlung cross section', XSS, JXS5 + NXS5, NXS6 )
    else:
        BremsstrahlungRatios = dismemberACE_miscModule.getData( 'Bremsstrahlung cross section', XSS, JXS5 + NXS5, NXS5 )
    for value in BremsstrahlungRatios: fOut.write( ' %20.12e\n' % value)

    if NXS16 == 3:
        fOut = dismemberACE_miscModule.openFile( 0, 'BremsstrahlungCrossSectionInterpolation' )
        BremsstrahlungCrossSectionInterpolation = dismemberACE_miscModule.getData( 'Bremsstrahlung cross section Interpolation', XSS, JXS5 + NXS5 + NXS6, NXS5 * NXS6 )
        for index1 in range(NXS5):
            for index2 in range(NXS6): fOut.write(' %20.12e' % BremsstrahlungCrossSectionInterpolation[index1*NXS6+index2])
            fOut.write('\n')
    else:
        fOut = dismemberACE_miscModule.openFile( 0, 'BremsstrahlungCrossSectionInterpolation' )
        BremsstrahlungCrossSectionInterpolation = dismemberACE_miscModule.getData('Bremsstrahlung cross section Interpolation', XSS, JXS6, 2 * NXS6)
        for index in range(NXS6):
            fOut.write(' %20.12e %20.12e' % (BremsstrahlungCrossSectionInterpolation[index], BremsstrahlungCrossSectionInterpolation[NXS6+index]))

    if JXS7 != 0 or JXS8 != 0:
        print('WARNING: JXS7 = %s and JXS8 = %s: currently only 0 is supported' % (JXS7, JXS8))

    fOut = dismemberACE_miscModule.openFile( 0, 'ratios1' )
    ratios1 = dismemberACE_miscModule.getData('ratios', XSS, JXS9, NXS9)
    for value in ratios1: fOut.write(' %20.12e\n' % value)

    fOut = dismemberACE_miscModule.openFile( 0, 'ratios2' )
    ratios2 = dismemberACE_miscModule.getData('ratios', XSS, JXS10, NXS10)
    for value in ratios2: fOut.write(' %20.12e\n' % value)

    fOut = dismemberACE_miscModule.openFile( 0, 'energyRatios' )
    energyRatios = dismemberACE_miscModule.getData('Energy ratio', XSS, JXS11, NXS11)
    for value in energyRatios : fOut.write(' %20.12e\n' % value)

    fOut = dismemberACE_miscModule.openFile( 0, 'BindingEnergies' )
    BindingEnergies = dismemberACE_miscModule.getData('Binding energies', XSS, JXS11 + NXS11, NXS11)
    for value in BindingEnergies: fOut.write(' %20.12e\n' % value)
