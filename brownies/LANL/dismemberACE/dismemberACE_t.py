# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module dismembers an ACE TNSL data file."""

from brownies.LANL.dismemberACE import dismemberACE_misc as dismemberACE_miscModule

def dismemberACE( args, NXS, JXS, XSS ) :

    IDPNI = NXS[2]                  # Inelastic scattering mode.
    NIL =   NXS[3]                  # Inelastic dimensioning parameter.
    NIEB =  NXS[4]                  # Number of inelastic exiting energies.
    IDPNC = NXS[5]                  # Elastic scattering mode (coherent = 4, incoherent != 4).
    NCL =   NXS[6]                  # Elastic dimensioning parameter.
    IFENG = NXS[7]                  # Secondary energy mode (discrete = 0, skewed = 1, continuous = 2).

    ITIE = JXS[1]                   # Location of inelastic energy table.
    ITIX = JXS[2]                   # Location of inelastic cross sections.
    ITXE = JXS[3]                   # Location of inelastic energy/angle distributions.
    ITCE = JXS[4]                   # Location of elastic energy table.
    ITCX = JXS[5]                   # Location of elastic cross sections.
    ITCA = JXS[6]                   # Location of elastic angular distributions.

    NE = int( dismemberACE_miscModule.getData( 'Number of inelastic energies', XSS, ITIE, 1 )[0] )
    energies = dismemberACE_miscModule.getData( 'Inelastic energies', XSS, ITIE + 1, NE )
    crossSections = dismemberACE_miscModule.getData( 'Inelastic cross section', XSS, ITIX, NE )

    MT = 4
    dismemberACE_miscModule.outputXYs1d( MT, 'crossSection', energies, crossSections )

    if( IFENG in [ 0, 1 ] ) :
        NILm1 = NIL - 1                  # Number of equal probable mus.
        offset = ITXE
        fOut = dismemberACE_miscModule.openFile( MT, 'energyAngular' )
        fOut.write( 'Secondary energy mode (discrete = 0, skewed = 1, continuous = 2) = %d\n' % NIEB )
        fOut.write( 'Number of equal probable outgoing energies  = %d\n' % NIEB )
        fOut.write( 'Number of equal probable mus = %d\n' % NILm1 )
        for index1 in range( NE ) :
            fOut.write( 'Incident energy = %20.12e\n' % energies[index1] )
            for index2 in range( NIEB ) :
                energyOut, pdf, cdf = dismemberACE_miscModule.getData( 'Energy out pdf and cdf', XSS, offset, 3 )
                fOut.write( '###     energy_out = %20.12e pdf = %20.12e cdf = %20.12e\n' % ( energyOut, pdf, cdf ) )
                offset += 3

                mus =  dismemberACE_miscModule.getData( 'Inelastic mus', XSS, offset, NILm1 )
                for mu in mus : fOut.write( '        %20.12e\n' % mu )
                offset += NILm1
    elif( IFENG == 2 ) :
        NILm1 = NIL - 1                  # Number of equal probable mus.
        fOut = dismemberACE_miscModule.openFile( MT, 'energyAngular' )
        locationOfDistributions = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'Location of distributions', XSS, ITXE, NE ) )
        numbersOfOutgoingEnergies = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'numbers of outgoing energies ', XSS, ITXE+NE, NE ) )
        for index1 in range( NE ) :
            locationOfDistribution = locationOfDistributions[index1]
            numberOfOutgoingEnergies = numbersOfOutgoingEnergies[index1]
            offset = locationOfDistribution + 1
            fOut.write( '\n\n### Incident energy = %20.12e, number of outgoing energies = %d, offset = %d\n' % ( energies[index1], numberOfOutgoingEnergies, offset ) )

            for index2 in range( numberOfOutgoingEnergies ) :
                energyOut, pdf, cdf = dismemberACE_miscModule.getData( 'Energy out pdf and cdf', XSS, offset, 3 )
                fOut.write( '###     energy_out = %20.12e pdf = %20.12e cdf = %20.12e\n' % ( energyOut, pdf, cdf ) )
                offset += 3

                mus = dismemberACE_miscModule.getData( 'Energy out distribution', XSS, offset, NILm1 )
                for index3 in range( NILm1 ) : fOut.write( '     %20.12e\n' % mus[index3] )
                offset += NILm1
    else :
        print( 'Unsupport IFENG (NXS[7]) flag = %s, skipping inelastic energy/angular data.' % IFENG )


    MT = 2
    if( ITCE != 0 ) :
        NE = int( dismemberACE_miscModule.getData( 'Number of elastic energies', XSS, ITCE, 1 )[0] )
        energies = dismemberACE_miscModule.getData( 'Elastic energies', XSS, ITCE + 1, NE )
        crossSections = dismemberACE_miscModule.getData( 'Elastic cross section', XSS, ITCX, NE )
        dismemberACE_miscModule.outputXYs1d( MT, 'crossSection', energies, crossSections )
        if( NXS[6] != -1 ) :
            NCLp1 = NCL + 1
            offset = ITCA
            fOut = dismemberACE_miscModule.openFile( MT, 'angular' )
            for index1 in range( NE ) :
                fOut.write( 'Incident energy = %20.12e\n' % energies[index1] )
                mus =  dismemberACE_miscModule.getData( 'Elastic mus', XSS, offset, NCLp1 )
                for mu in mus : fOut.write( '   %20.12e\n' % mu )
                offset += NCLp1
