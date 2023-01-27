# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.LANL.dismemberACE import dismemberACE_misc as dismemberACE_miscModule

def dismemberACE( args, NXS, JXS, XSS ) :

    NES = NXS[3]                    # Number of points of energy grid.
    NTR = NXS[4]                    # Number of reactions excluding elastic.
    NR =  NXS[5]                    # Number of reactions having secondary neutrons excluding elastic.

    energyGrid = dismemberACE_miscModule.getData( 'energyGrid', XSS, JXS[1], NES )
    energy1 = -1
    printHeader = 'Non-monotonic energy grid\n'
    for i1, energy2 in enumerate( energyGrid ) :                # Check that energy grid is monotonic.
        if( energy1 >= energy2 ) :
            print( '%s    index = %7d   energy1 = %.17e  energy2 = %.17e   (energy2 - energy1) = %e' % ( printHeader, i1, energy1, energy2, energy2 - energy1 ) )
            printHeader = ''
        energy1 = energy2

    dismemberACE_miscModule.outputXYs1d( 0, 'totalCrossSection',        energyGrid, dismemberACE_miscModule.getData( 'totalCrossSection',       XSS, JXS[1], NES,     NES ) )
    dismemberACE_miscModule.outputXYs1d( 0, 'absorptionCrossSection',   energyGrid, dismemberACE_miscModule.getData( 'absorptionCrossSection',  XSS, JXS[1], NES, 2 * NES ) )
    dismemberACE_miscModule.outputXYs1d( 2, 'crossSection',             energyGrid, dismemberACE_miscModule.getData( 'crossSection',            XSS, JXS[1], NES, 3 * NES ) )
    dismemberACE_miscModule.outputXYs1d( 0, 'averageHeatingNumbers',    energyGrid, dismemberACE_miscModule.getData( 'averageHeatingNumbers',   XSS, JXS[1], NES, 4 * NES ) )

    dismemberACE_miscModule.nu_bar( JXS[2], XSS )
    if( JXS[24] > 0 ) : dismemberACE_miscModule.nu_bar2( 'nubar_delayed', JXS[24], XSS )
    dismemberACE_miscModule.delayedNeutronData( NXS, JXS, XSS )

    MTR = JXS[3]                    # Location of MT array.
    MTs = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'MTs', XSS, MTR, NTR ) )
    LQR = JXS[4]                    # Location of Q-value array.
    Qs = dismemberACE_miscModule.getData( 'Qs', XSS, LQR ,NTR )
    TYR = JXS[5]                    # Location of reaction type array.
    Types = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'Types', XSS, TYR, NTR ) )

    LSIG = JXS[6]                   # Location of table of cross section locators.
    crossSectionLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LSIG, NTR ) )
    SIG =  JXS[7]                   # Location of cross sections.

    LAND = JXS[8]                   # Location of table of angular distribution locators.
    angularLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LAND, NR + 1 ) )
    AND =  JXS[9]                   # Location of angular distribution.

    LDLW = JXS[10]                   # Location of table of energy distribution locators.
    energyLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LDLW, NR ) )
    DLW =  JXS[11]                   # Location of energy distribution.

    MTPsNs = {}
    NTRP  = NXS[6]                  # Number of photon production reactions.
    if( NTRP > 0 ) :
        GPD   = JXS[12]

        if( GPD > 0 ) :
            offset = GPD
            dismemberACE_miscModule.outputXYs1d( 0, 'totalPhotonProductionCrossSection', energyGrid, dismemberACE_miscModule.getData( 'totalPhotonProductionCrossSection', XSS, offset, NES ) )
            offset += NES
            if( False ) :           # Does not seem to exists for data_endfb7.1/O_016_300K.ace
                fOut = dismemberACE_miscModule.openFile( 0, 'outgoingPhotonEnergies' )
                for i1 in range( 30 ) :
                    fOut.write( "\n\n### energy index %d\n" % i1 )
                    outgoingPhotonEnergies = dismemberACE_miscModule.getData( 'outgoingPhotonEnergies', XSS, offset, 20 )
                    for energy in outgoingPhotonEnergies : fOut.write( "    %20.12e\n" % energy )
                    offset += 20
                fOut.close( )

        MTRP  = JXS[13]
        MTPs  = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'MTPs', XSS, MTRP, NTRP ) )
        for MTP in MTPs :
            MT = MTP // 1000
            if( MT not in MTPsNs ) : MTPsNs[MT] = []
            MTPsNs[MT].append( MTP )

        LSIGP = JXS[14]
        photonProductionCrossSectionLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LSIGP, NTRP ) )
        SIGP  = JXS[15]

        LANDP = JXS[16]
        photonAngularLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LANDP, NTRP ) )
        ANDP  = JXS[17]

        LDLWP = JXS[18]
        photonEnergyLocators = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'crossSectionLocators', XSS, LDLWP, NTRP ) )
        DLWP  = JXS[19]

    if( args.addresses ) : print( "LAND = %d:  AND = %d:  LDLW = %d:  DLW = %d" % ( LAND, AND, LDLW, DLW ) )
#
# Elastic reaction.
#
    fOut = dismemberACE_miscModule.openFile( 2, 'info' )
    fOut.write( 'type = %d\n' % -1 )
    fOut.write( 'Q = %20.12e\n' % 0.0 )
    fOut.close( )

    dismemberACE_miscModule.dismemberAngular( 0, 2, -1, angularLocators, AND, XSS, 'neutron' )         # Special case for elastic.

#
# Non-elastic reaction.
#
    for i1, MT in enumerate( MTs ) :
        if( args.verbose > 1 ) : print( "MT = %3d" % MT )
        locator = SIG + crossSectionLocators[i1] - 1

        fOut = dismemberACE_miscModule.openFile( MT, 'info' )
        fOut.write( 'type = %d\n' % Types[i1] )
        fOut.write( 'Q = %20.12e\n' % Qs[i1] )
        if( args.addresses ) : fOut.write( 'locator = %d\n' % locator )
        line, column = dismemberACE_miscModule.lineColumn1BaseForXSS( locator )
        if( args.addresses ) : fOut.write( 'line, column = %d, %d\n' % ( line, column ) )
        fOut.close( )

        offset, length = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'cross section info for MT = %d' % MT, XSS, locator, 2 ) )
        crossSection = dismemberACE_miscModule.getData( 'cross section for MT = %d' % MT, XSS, locator + 2, length )
        dismemberACE_miscModule.outputXYs1d( MT, 'crossSection', energyGrid, crossSection, offset = offset, length = length )

        if( 4 < MT < 100 ) :
            dismemberACE_miscModule.dismemberAngular( i1 + 1, MT, Types[i1], angularLocators, AND, XSS, 'neutron' )
            dismemberACE_miscModule.dismemberEnergy( 0, MT, Types[i1], DLW, energyLocators[i1], XSS, 'neutron' )

        TY = abs( Types[i1] )
        if( TY > 100 ) :
            fOut = dismemberACE_miscModule.openFile( MT, 'multiplicityVsEnergy', subDir = 'neutron' )

            offset = DLW + TY - 101
            offset += dismemberACE_miscModule.interpolationData( fOut, offset, XSS, 'photonProductionCrossSection' )

            NE = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'NE (photonProductionCrossSection)', XSS, offset, 1 ) )[0]
            fOut.write( '### NE = %d\n' % NE )
            energies = dismemberACE_miscModule.getData( 'energies (photonProductionCrossSection)', XSS, offset + 1,      NE )
            yields   = dismemberACE_miscModule.getData( 'P(E) (photonProductionCrossSection)',     XSS, offset + 1 + NE, NE )
            for i1 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i1], yields[i1] ) )

            fOut.close( )

    for MT in sorted( MTPsNs ) :
        MTPNs = MTPsNs[MT]
        for MTPN in MTPNs :
            i1 = MTPs.index( MTPN )
            locator = SIGP + photonProductionCrossSectionLocators[i1] - 1
            dismemberACE_miscModule.photonProductionCrossSection( MT, MTPN, SIGP, locator, XSS, energyGrid )
            dismemberACE_miscModule.dismemberAngular( i1, MT, 0, photonAngularLocators, ANDP, XSS, 'photon', MTP = MTPN )
            dismemberACE_miscModule.dismemberEnergy( 0, MT, 1, DLWP, photonEnergyLocators[i1], XSS, 'photon', MTP = MTPN )

    YP = JXS[20]
    if( YP > 0 ) :
        NYP = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'YP', XSS, YP, 1 ) )[0]
        MTYs = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'MTY', XSS, YP + 1, NYP ) )
        fOut = dismemberACE_miscModule.openFile( 0, 'YP' )
        fOut.write( 'NYP = %d\n' % NYP )
        for MTY in MTYs :fOut.write( '%d\n' % MTY )
        fOut.close( )

    LUNR = JXS[23]              # Location of probability tables
    if( LUNR > 0 ) :
        fOut = dismemberACE_miscModule.openFile( 0, 'URR_probabilityTables' )
        offset = LUNR
        N, M, INT, ILF, IOA, IFF = dismemberACE_miscModule.toIntegers( dismemberACE_miscModule.getData( 'N (URR)', XSS, LUNR, 6 ) )
        offset += 6
        fOut.write( '### N = %d\n'   % N )
        fOut.write( '### M = %d\n'   % M )
        fOut.write( '### INT = %d\n' % INT )
        fOut.write( '### ILF = %d\n' % ILF )
        fOut.write( '### IOA = %d\n' % IOA )
        fOut.write( '### IFF = %d\n' % IFF )

        energies = dismemberACE_miscModule.getData( 'energies (URR)', XSS, offset, N )
        offset += N
        fOut.write( '### Energies\n' )
        for energy in energies : fOut.write( '%20.12e\n' % energy )

        Ps = dismemberACE_miscModule.getData( 'P(i,j,k) (URR)', XSS, offset, N * 6 * M )
        fOut.write( '\n\n### Probality tables\n' )
        URR_offset = 0
        for energyIndex in range(N):
            enetryStr = ' at energy %s' % energies[energyIndex]
            for type in ['cumulative probability', 'total', 'elastic', 'fission', 'capture', 'heating']:
                typeStr = '    # %s%s' % (type, enetryStr)
                for index in range(M):
                    fOut.write( '%20.12e%s\n' % (Ps[URR_offset], typeStr))
                    URR_offset += 1
                    typeStr = ''
                enetryStr = ''

        fOut.close( )
