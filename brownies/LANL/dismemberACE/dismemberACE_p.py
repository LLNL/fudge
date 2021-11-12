# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.LANL.dismemberACE import dismemberACE_misc as dismemberACE_miscModule

def dismemberACE( args, NXS, JXS, XSS ) :

    NES = NXS[3]                    # Number of points of energy grid.

# ESZG Block.
    logEnergyGrid = dismemberACE_miscModule.getData( 'logEnergyGrid', XSS, JXS[1], NES )

    dismemberACE_miscModule.outputXYs1d( 504, 'incoherentCrossSection',     logEnergyGrid, dismemberACE_miscModule.getData( 'incoherentCrossSection',     XSS, JXS[1], NES,     NES ) )
    dismemberACE_miscModule.outputXYs1d( 502, 'coherentCrossSection',       logEnergyGrid, dismemberACE_miscModule.getData( 'coherentCrossSection',       XSS, JXS[1], NES, 2 * NES ) )
    dismemberACE_miscModule.outputXYs1d( 523, 'photoelectricCrossSection',  logEnergyGrid, dismemberACE_miscModule.getData( 'photoelectricCrossSection',  XSS, JXS[1], NES, 3 * NES ) )
    dismemberACE_miscModule.outputXYs1d( 516, 'pairProductionCrossSection', logEnergyGrid, dismemberACE_miscModule.getData( 'pairProductionCrossSection', XSS, JXS[1], NES, 4 * NES ) )

# JINC Block.
    NE = 21         # Fixed number of incoherent scattering function points.
    incoherentScatteringFunction = dismemberACE_miscModule.getData( 'incoherentScatteringFunction', XSS, JXS[2], NE )
    fOut = dismemberACE_miscModule.openFile( 504, 'incoherentScatteringFunction' )
    for value in incoherentScatteringFunction : fOut.write( '%20.12e\n' % value )
    fOut.close( )

    NE = 55         # Fixed number of coherent form factor points.
    integratedCoherentFormFactor = dismemberACE_miscModule.getData( 'integratedCoherentFormFactor', XSS, JXS[3], NE )
    fOut = dismemberACE_miscModule.openFile( 502, 'integratedCoherentFormFactor' )
    for value in integratedCoherentFormFactor : fOut.write( '%20.12e\n' % value )
    fOut.close( )

# JCOH Block.
    coherentFormFactor = dismemberACE_miscModule.getData( 'coherentFormFactor', XSS, JXS[3] + NE, NE )
    fOut = dismemberACE_miscModule.openFile( 502, 'coherentFormFactor' )
    for value in coherentFormFactor : fOut.write( '%20.12e\n' % value )
    fOut.close( )

# JFLO Block.
    if( NXS[4] != 0 ) :
        fOut = dismemberACE_miscModule.openFile( 0, 'fluorescence' )

        data = dismemberACE_miscModule.getData( 'Edge energies', XSS, JXS[4], NXS[4] )
        dismemberACE_miscModule.writeVector( fOut, 'Edge energies', '%20.12e', data )

        data = dismemberACE_miscModule.getData( 'Relative probabilities of ejection', XSS, JXS[4] + NXS[4] , NXS[4] )
        dismemberACE_miscModule.writeVector( fOut, 'Relative probabilities of ejection', '%20.12e', data )

        data = dismemberACE_miscModule.getData( 'Yields', XSS, JXS[4] + 2 * NXS[4], NXS[4] )
        dismemberACE_miscModule.writeVector( fOut, 'Yields', '%20.12e', data )

        data = dismemberACE_miscModule.getData( 'Fluorescent energies', XSS, JXS[4] + 3 * NXS[4], NXS[4] )
        dismemberACE_miscModule.writeVector( fOut, 'Fluorescent energies', '%20.12e', data )

        fOut.close( )

# LHNM Block.
    averageHeating = dismemberACE_miscModule.getData( 'averageHeating', XSS, JXS[5], NXS[3] )
    fOut = dismemberACE_miscModule.openFile( 0, 'averageHeating' )
    for value in averageHeating : fOut.write( '%20.12e\n' % value )
    fOut.close( )

    NXS5 = NXS[5]
    if( NXS5 != 0 ) :
        fOut = dismemberACE_miscModule.openFile( 0, 'shells' )

# LNEPS Block.
        data = dismemberACE_miscModule.getData( 'electrons', XSS, JXS[6], NXS5 )
        fOut.write( '# Electrons per shell\n' )
        for value in data : fOut.write( '    %20.12e\n' % value )

# LBEPS Block.
        data = dismemberACE_miscModule.getData( 'bindindEnergies', XSS, JXS[7], NXS5 )
        fOut.write( '# Binding energies per shell\n' )
        for value in data : fOut.write( '    %20.12e\n' % value )

# LPIPS Block.
        data = dismemberACE_miscModule.getData( 'probabilityOfInteraction', XSS, JXS[8], NXS5 )
        fOut.write( 'Probability of interaction per shell\n' )
        for value in data : fOut.write( '    %20.12e\n' % value )

# LSWD and SWD blocks.
        LSWD = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'electronsPerShell', XSS, JXS[9], NXS5 ) )
        for index, address in enumerate( LSWD ) :
            start = JXS[10] + address - 1
            fOut = dismemberACE_miscModule.openFile( 0, '%3.3d' % index, subDir = 'shell' )

            JJ = int( dismemberACE_miscModule.getData( 'JJ', XSS, start, 1 )[0] )
            NE = int( dismemberACE_miscModule.getData( 'NE', XSS, start + 1, 1 )[0] )
            start += 2

            data = dismemberACE_miscModule.getData( 'momentumComptonProfileShell', XSS, start, NE )
            fOut.write( '# Momentum parameter for compton profile\n' )
            for value in data : fOut.write( '    %20.12e\n' % value )

            data = dismemberACE_miscModule.getData( 'pdfComptonProfileShell', XSS, start, NE, offset = NE )
            fOut.write( '# pdf for compton profile\n' )
            for value in data : fOut.write( '    %20.12e\n' % value )

            data = dismemberACE_miscModule.getData( 'cdfComptonProfileShell', XSS, start, NE, offset = 2 * NE )
            fOut.write( '# cdf for compton profile\n' )
            for value in data : fOut.write( '    %20.12e\n' % value )

            fOut.close( )
