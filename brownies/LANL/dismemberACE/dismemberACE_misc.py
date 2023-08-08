# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

cites = []
args = None

def initialize( _args, XSS ) :

    global args, cites

    args = _args
    cites = len( XSS ) * [ 0 ]

def openFile( MT, name, subDir = None, MTP = None ) :

    return( openFile2( args, MT, name, subDir = subDir, MTP = MTP ) )

def openFile2( args, MT, name, subDir = None, MTP = None ) :

    if( MT > 0 ) :
        path = os.path.join( args.output, "MTs", "%.3d" % MT )
    else :
        path = os.path.join( args.output )

    if( subDir is not None ) : path = os.path.join( path, subDir )
    if( MTP is not None ) : path = os.path.join( path, str( MTP ) )
    if( not( os.path.exists( path ) ) ) : os.makedirs( path )
    path = os.path.join( path, name )
    return( open( path, 'w' ) )        

def lineColumn1BaseForXSS( offset ) :

    offset -= 1
    line = offset // 4 + 13
    column = offset % 4 + 1
    return( line, column  )

def get8Integers( line ) :

    values = []
    for i1 in range( 8 ) :
        values.append( int( line[:9] ) )
        line = line[9:]

    return( values )

class ListBase1 :

        def __init__( self, values ) : self.values = values
        def __len__( self ) : return( len( self.values ) )
        def __getitem__( self, index ) : return( self.values[index-1] )

def get4Floats( line, number, lineNumber ) :

    initialLine = line
    count = min( 4, number )
    number -= count
    values = []
    for i1 in range( count ) :
        if( len( line ) < 2 ) : break
        try :
            values.append( float( line[:20] ) )
        except :
            print( lineNumber )
            print( initialLine[:-1] )
            print( "HI <%s>" % line )
            raise
        line = line[20:]

    return( number, values )

def getData( label, XSS, start, numberOfPoints, offset = 0, cite = True ) :

    global cites

    start += offset - 1
    if( cite ) :
        total = sum( cites[start:start+numberOfPoints] )
        if( total > 0 ) : print( '    citings already in range %9d to %9d of length %9d: %s' % ( start, start + numberOfPoints, numberOfPoints, label ) )
        for i1 in range( start, start + numberOfPoints ) : cites[i1] += 1
    return( XSS.values[start:start+numberOfPoints] )

def toIntegers( values ) :

    return( [ int( value ) for value in values ] )

def outputList(name, values):

    fOut = openFile(-1, name)
    for value in values:
        fOut.write('%s\n' % value)
    fOut.close( )

def outputXYs1d( MT, name, xs, ys, offset = 0, length = -1, subDir = None, MTP = None ) :

    if( length == -1 ) : length = len( xs )
    xs = xs[offset:offset+length]

    fOut = openFile( MT, name, subDir = subDir, MTP = MTP )
    for i1, x in enumerate( xs ) : fOut.write( "%20.12e %20.12e\n" % ( x, ys[i1] ) )
    fOut.close( )

def writeVector( fOut, label, format, data ) :

    fOut.write( '%s\n' % label )
    for datum in data : fOut.write( ( format % datum ) + '\n' )

def interpolationData( fOut, offset, XSS, distType ) :

    NR   = toIntegers( getData( 'NR (%s)' % distType,  XSS, offset, 1 ) )[0]
    if( NR != 0 ) :
        NBT = toIntegers( getData( 'NBT (%s)\n' % distType, XSS, offset + 1,      NR ) )
        INT = toIntegers( getData( 'INT (%s)\n' % distType, XSS, offset + 1 + NR, NR ) )
        for i1 in range( NR ) : fOut.write( '### NBT = %3d   INT = %2d\n' % ( NBT[i1], INT[i1] ) )

    return( 1 + 2 * NR )

def nu_bar( NU, XSS ) :

    if( NU == 0 ) : return

    KNU = toIntegers( getData( 'KNU', XSS, NU, 1 ) )[0]

    if( KNU > 0 ) :
        nu_bar2( 'nubar_prompt_or_total', NU, XSS, LNU = KNU )
    else :
        nu_bar2( 'nubar_prompt', NU + 1, XSS )
        nu_bar2( 'nubar_total', NU + abs( KNU ) + 1, XSS )

def nu_bar2( label, NU, XSS, LNU = 0 ) :

    fOut = openFile( 18, label )

    if( LNU == 0 ) : LNU = toIntegers( getData( 'LNU (%s)' % label,  XSS, NU, 1 ) )[0]
    fOut.write( '### LNU = %d\n' % LNU )
    offset = NU + 1

    if( LNU == 1 ) :
        NC = toIntegers( getData( 'NC (%s)' % label, XSS, offset, 1 ) )[0]
        fOut.write( '### NC = %d\n' % NC )
        Cs =  getData( 'Cs (%s)' % label, XSS, offset + 1,      NC )
        for C in Cs : fOut.write( '    %20.12e\n' % C )
    elif( LNU == 2 ) :
        offset += interpolationData( fOut, offset, XSS, 'energy' )
        NE = toIntegers( getData( 'NE (%s)' % label, XSS, offset, 1 ) )[0]
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (%s)' % label, XSS, offset + 1,      NE )
        nu_bar   = getData( 'nu_bar (%s)' % label,   XSS, offset + 1 + NE, NE )
        for i1 in range( NE ) : fOut.write( '    %20.12e %20.12e\n' % ( energies[i1], nu_bar[i1] ) )
    else :
        raise Exception( 'Invalid LNU = %d' % LNU )

    fOut.close( )

def delayedNeutronData( NXS, JXS, XSS ) :

    NPCR = NXS[8]
    if( NPCR == 0 ) : return

    subDir = 'delayedNeutron'
    BDD    = JXS[25]
    LDNEDL = JXS[26]
    DNEDL = JXS[27]
    energyLocators = toIntegers( getData( 'energyLocators (%s)' % subDir, XSS, LDNEDL, NPCR ) )
    offset = BDD
    for i1 in range( NPCR ) :
        label = 'DEC[%d]' % i1
        fOut = openFile( 18, 'probability', subDir = subDir, MTP = i1 )

        DEC = getData( label, XSS, offset, 1 )[0]
        fOut.write( '### fission delayed neutron %20.12e\n' % DEC )
        offset += 1

        offset += interpolationData( fOut, offset, XSS, 'energy' )

        NE = toIntegers( getData( 'NE (%s)' % label, XSS, offset, 1 ) )[0]
        offset += 1
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (%s)' % label, XSS, offset,      NE )
        Ps       = getData( 'nu_bar (%s)' % label,   XSS, offset + NE, NE )
        for i2 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i2], Ps[i2] ) )
        offset += 2 * NE

        dismemberEnergy( 0, 18, 1, DNEDL, energyLocators[i1], XSS, subDir = subDir, MTP = i1 )

    fOut.close( )

def table_12( MT, AND, locator, XSS, subDir ) :
    """MCNP angular data."""

    fOut = openFile( MT, 'angular', subDir )

    NE = toIntegers( getData( 'NE (angular)', XSS, locator, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    energies = getData( 'energies (angular)', XSS, locator + 1, NE )
    LCs = toIntegers( getData( 'energies (angular)', XSS, locator + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( 'energy = %20.12e\n' % energies[i1] )
        locator2 = AND + abs( LCs[i1] ) - 1
        if( args.addresses ) : fOut.write( '    locator = %d (%d)\n' % ( LCs[i1], locator2 ) )
        if( LCs[i1] < 0 ) :
            interpolation, NP = toIntegers( getData( 'interpolation, NP (angular)', XSS, locator2, 2 ) )
            locator2 += 2
            interpolationStr = { 1 : 'flat', 2 : 'lin-lin' }
            fOut.write( '    interpolation = %s (%s)\n' % ( interpolation, interpolationStr[interpolation] ) )
            fOut.write( '    NP = %s\n' % NP )

            mus = getData( 'mus', XSS, locator2, NP )
            locator2 += NP
            pdf = getData( 'mus', XSS, locator2, NP )
            locator2 += NP
            cdf = getData( 'mus', XSS, locator2, NP )
            for i2 in range( NP ) : fOut.write( '        %20.12e %20.12e %20.12e\n' % ( mus[i2], pdf[i2], cdf[i2] ) )

        elif( LCs[i1] == 0 ) :                  # Isotopic
            fOut.write( '        isotropic' )

        else :                                  # 32 equal probable bins.
            epbs = getData( 'epbs (angular)', XSS, AND + LCs[i1] - 1, 33 )
            for epb in epbs : fOut.write( '        %20.12e\n' % epb )

    fOut.close( )

def dismemberAngular( index, MT, Type, LAND, AND, XSS, subDir, MTP = None ) :

    frame = 'lab'
    if( Type < 0 ) : frame = 'center of mass'
    locator = 0
    if( LAND[index] > 0 ) : locator = AND + LAND[index] - 1

    fOut = openFile( MT, 'info', subDir = subDir, MTP = MTP )
    fOut.write( 'Type = %d\n' % Type )          # abs( Type ) is multiplicity. If 19 fission. If > 100 is non-fission energy dependent.
    fOut.write( '    multiplicity = %d\n' % abs( Type ) )
    fOut.write( '    frame = %s\n' % frame )
    if( args.addresses ) : fOut.write( 'LAND[%d] = %s\n' % ( index, LAND[index] ) )
    if( args.addresses ) : fOut.write( 'locator (1 based) = %s\n' % locator )
    fOut.close( )

    if( LAND[index] < 0 ) :                         # Kalback/Mann data given in the energy section.
        pass
    elif( LAND[index] == 0 ) :
        fOut = openFile( MT, 'angular', subDir = subDir, MTP = MTP )
        fOut.write( '### isotropic\n' )
        fOut.close( )
    else :
        table_12( MT, AND, locator, XSS, 'neutron' )

def dismemberEnergy( level, MT, Type, DLW, offset, XSS, subDir, MTP = None ) :

    fOut = openFile( MT, 'energy_%d' % level, subDir = subDir, MTP = MTP )

    if( args.verbose > 2 ) : print( '            DLW = %d  offset = %d  start = %d: %s' % ( DLW, offset, DLW + offset - 1, subDir ) )
    offset += DLW - 1

    fOut.write( 'Type = %d\n' % Type )
    fOut.write( 'DLW = %d\n' % DLW )

    LNW = toIntegers( getData( 'LNW (energy)',   XSS, offset,     1 ) )[0]
    fOut.write( 'LNW = %d\n' % LNW )
    LAW  = toIntegers( getData( 'LAW (energy)',  XSS, offset + 1, 1 ) )[0]
    fOut.write( 'LAW = %d\n' % LAW )
    IDAT = toIntegers( getData( 'IDAT (energy)', XSS, offset + 2, 1 ) )[0]
    fOut.write( 'IDAT = %d\n' % IDAT )
    offset += 3

    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( 'NE = %d\n' % NE )
    energies =  getData( 'energies (energy)', XSS, offset + 1,      NE )
    Ps       =  getData( 'P(E) (energy)',     XSS, offset + 1 + NE, NE )
    for i1 in range( NE ) : fOut.write( 'energy = %20.12e %7d\n' % ( energies[i1], Ps[i1] ) )

    if( args.verbose > 2 ) : print( '            %sLAW = %d' % ( level * '    ', LAW ) )
# Missing LAWs 1, 22 (UK), 24 (UK)
    if( LAW == 2 ) :
        LAW2( fOut, DLW, IDAT, XSS )
    elif( LAW == 3 ) :
        LAW3( fOut, DLW, IDAT, XSS )
    elif( LAW == 4 ) :
        LAW4( fOut, DLW, IDAT, XSS )
    elif( LAW == 5 ) :
        LAW5( fOut, DLW, IDAT, XSS )
    elif( LAW == 7 ) :
        LAW7( fOut, DLW, IDAT, XSS )
    elif( LAW == 9 ) :
        LAW9( fOut, DLW, IDAT, XSS )
    elif( LAW == 11 ) :
        LAW11( fOut, DLW, IDAT, XSS )
    elif( LAW == 44 ) :
        LAW44_KalbachMann( fOut, DLW, IDAT, XSS )
    elif( LAW == 61 ) :
        LAW61( fOut, DLW, IDAT, XSS )
    elif( LAW == 66 ) :
        LAW66( fOut, DLW, IDAT, XSS )
    elif( LAW == 67 ) :
        LAW67( fOut, DLW, IDAT, XSS )
    else :
        print( "Unsupported LAW = %s" % LAW )
    fOut.close( )

    if( LNW != 0 ) :
        dismemberEnergy( level + 1, MT, Type, DLW, LNW, XSS, subDir, MTP = MTP )

def LAW2( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW2\n' )

    offset += DLW - 1
    LDAT1, LDAT2 = getData( 'LAW2 (energy)', XSS, offset, 2 )
    fOut.write( 'LDAT[1] = %d\n' % int( LDAT1 ) )
    fOut.write( 'LDAT[2] = %20.12e\n' % LDAT2 )

def LAW3( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW3\n' )

    offset += DLW - 1
    LDAT1, LDAT2 = getData( 'LAW3 (energy)', XSS, offset, 2 )
    fOut.write( 'LDAT[1] = %20.12e\n' % LDAT1 )
    fOut.write( 'LDAT[2] = %20.12e\n' % LDAT2 )

def LAW4( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW4\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )
    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '### INTT = %s\n' % INTT )
        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '### NP = %s\n' % NP )
        EPs = getData( 'EPS (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e\n' % ( EPs[i2], pdf[i2], cdf[i2] ) )

def LAW5( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW5\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )
    offset += 2 * NE

    NET = toIntegers( getData( 'NET (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NET = %d\n' % NET )
    offset += 1

    energies = getData( 'energies (energy)', XSS, offset,     NET )
    pdf = getData( 'pdf (energy)', XSS, offset + NET,     NET )
    cdf = getData( 'cdf (energy)', XSS, offset + 2 * NET, NET )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e %20.12e\n' % ( energy, pdf[i1], cdf[i1] ) )

def LAW7( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW7\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NE, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )

def LAW9( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW9\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NE = %d\n' % NE )
    offset += 1

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NE, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    energies = getData( 'energies (energy)', XSS, offset,      NE )
    Ts =       getData( 'Ts (energy)',       XSS, offset + NE, NE )
    for i1, energy in enumerate( energies ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Ts[i1] ) )

def LAW11( fOut, DLW, offset, XSS ) :

    fOut.write( '### LAW11\n' )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NEa = toIntegers( getData( 'NEa (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NEa = %d\n' % NEa )
    offset += 1
    energies_a = getData( 'energies_a (energy)', XSS, offset,       NEa )
    As =         getData( 'As (energy)',         XSS, offset + NEa, NEa )
    offset += 2 * NEa

    offset += interpolationData( fOut, offset, XSS, 'energy' )
    NEb = toIntegers( getData( 'NEb (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NEb = %d\n' % NEb )
    offset += 1
    energies_b = getData( 'energies_b (energy)', XSS, offset,       NEb )
    Bs =         getData( 'Bs (energy)',         XSS, offset + NEb, NEb )

    U = toIntegers( getData( 'U (energy)', XSS, offset + 2 * NEb, 1 ) )[0]
    fOut.write( '### U = %20.12e\n' % U )

    fOut.write( '### as\n' )
    for i1, energy in enumerate( energies_a ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, As[i1] ) )

    fOut.write( '\n\n### bs\n' )
    for i1, energy in enumerate( energies_b ) : fOut.write( '    %20.12e %20.12e\n' % ( energy, Bs[i1] ) )

def LAW44_KalbachMann( fOut, DLW, offset, XSS ) :

    if( args.verbose > 2 ) : print( '            offset = %d  start = %d' % ( offset, DLW + offset - 1 ) )

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '###    INTT = %s\n' % INTT )

        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '###    NP = %s\n' % NP )
        EPs = getData( 'EPS (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        Rs  = getData( 'Rs (energy)',  XSS, offset + 2 + 3 * NP, NP )
        As  = getData( 'As (energy)',  XSS, offset + 2 + 4 * NP, NP )
        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e %20.12e %20.12e\n' % ( EPs[i2], pdf[i2], cdf[i2], Rs[i2], As[i2] ) )

def LAW61( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTT = toIntegers( getData( 'INTT (energy)', XSS, offset, 1 ) )[0]
        fOut.write( '###    INTT = %s\n' % INTT )

        NP = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
        fOut.write( '###    NP = %s\n' % NP )

        EPs = getData( 'EPs (energy)', XSS, offset + 2,          NP )
        pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP,     NP )
        cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP, NP )
        LC2s  = toIntegers( getData( 'LC2s (energy)',  XSS, offset + 2 + 3 * NP, NP ) )

        for i2 in range( NP ) :
            fOut.write( '    %20.12e %20.12e %20.12e %8d\n' % ( EPs[i2], pdf[i2], cdf[i2], LC2s[i2] ) )
        for i2 in range( NP ) :
            fOut.write( '\n\n' )
            fOut.write( '###     energy_out = %20.12e\n' % EPs[i2] )

            if( LC2s[i2] == 0 ) :
                fOut.write( '###         isotopic\n' )
            else :
                offset = DLW + LC2s[i2] - 1
                JJ = toIntegers( getData( 'JJ (energy)', XSS, offset, 1 ) )[0]
                fOut.write( '###        JJ = %s\n' % INTT )

                NP2 = toIntegers( getData( 'INTT (energy)', XSS, offset + 1, 1 ) )[0]
                mus = getData( 'mus (energy)', XSS, offset + 2,           NP2 )
                pdf = getData( 'pdf (energy)', XSS, offset + 2 + NP2,     NP2 )
                cdf = getData( 'cdf (energy)', XSS, offset + 2 + 2 * NP2, NP2 )
                for i3 in range( NP2 ) :
                    fOut.write( '        %20.12e %20.12e %20.12e\n' % ( mus[i3], pdf[i3], cdf[i3] ) )

def LAW66( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    NPSX = toIntegers( getData( 'NPSX (energy)', XSS, offset, 1 ) )[0]
    fOut.write( '### NPSX = %d\n' % NPSX )
    Ap = getData( 'Ap (energy)', XSS, offset + 1, 1 )[0]
    fOut.write( '### Ap = %20.12e\n' % Ap )

def LAW67( fOut, DLW, offset, XSS ) :

    offset += DLW - 1
    offset += interpolationData( fOut, offset, XSS, 'energy' )

    NE = toIntegers( getData( 'NE (energy)', XSS, offset, 1 ) )[0]
    energies =        getData( 'energies (energy)', XSS, offset + 1,      NE )
    LCs = toIntegers( getData( 'energies (energy)', XSS, offset + 1 + NE, NE ) )
    for i1 in range( NE ) :
        fOut.write( '\n\n' )
        fOut.write( '### energy = %20.12e\n' % energies[i1] )

        offset = DLW + LCs[i1] - 1
        INTMU, NMU = toIntegers( getData( 'INTMU (energy)', XSS, offset, 2 ) )
        offset += 2
        fOut.write( '###   INTMU = %d\n' % INTMU )
        fOut.write( '###   NMU = %d\n' % NMU )

        mus  = getData( 'mus (energy)', XSS, offset, NMU )
        LMU = toIntegers( getData( 'LMU (energy)', XSS, offset + NMU, NMU ) )

        for i2 in range( NMU ) :
            fOut.write( '\n###   mu = %20.12e\n' % mus[i2] )

            offset = DLW + LMU[i2] - 1
            INTEP, NPEP = toIntegers( getData( 'INTEP (energy)', XSS, offset, 2 ) )
            offset += 2
            fOut.write( '###     INTEP = %d\n' % INTEP )
            fOut.write( '###     NPEP = %d\n' % NPEP )

            eps = getData( 'eps (energy)', XSS, offset,            NPEP )
            pdf = getData( 'pdf (energy)', XSS, offset +     NPEP, NPEP )
            cdf = getData( 'cdf (energy)', XSS, offset + 2 * NPEP, NPEP )
            for i3 in range( NPEP ) :
                fOut.write( '        %20.12e %20.12e %20.12e\n' % ( eps[i3], pdf[i3], cdf[i3] ) )

def photonProductionCrossSection( MT, MTPN, SIGP, offset, XSS, energyGrid ) :

    fOut = openFile( MT, 'crossSection', subDir = 'photon', MTP = MTPN )
    MFTYPE = toIntegers( getData( 'MFTYPE', XSS, offset, 1 ) )[0]
    fOut.write( '### MFTYPE = %d\n' % MFTYPE )
    offset += 1

    if( MFTYPE in [ 12, 16 ] ) :
        MTMULT = toIntegers( getData( 'MFTYPE', XSS, offset, 1 ) )[0]
        offset += 1
        fOut.write( '### MTMULT = %d\n' % MTMULT )

        offset += interpolationData( fOut, offset, XSS, 'photonProductionCrossSection' )

        NE = toIntegers( getData( 'NE (photonProductionCrossSection)', XSS, offset, 1 ) )[0]
        fOut.write( '### NE = %d\n' % NE )
        energies = getData( 'energies (photonProductionCrossSection)', XSS, offset + 1,      NE )
        yields   = getData( 'P(E) (photonProductionCrossSection)',     XSS, offset + 1 + NE, NE )
        for i1 in range( NE ) : fOut.write( '%20.12e %20.12e\n' % ( energies[i1], yields[i1] ) )

    elif( MFTYPE == 13 ) :
        IE, NE = toIntegers( getData( 'cross section info for MT = %d' % MT, XSS, offset, 2 ) )
        energyGrid = energyGrid[IE:IE+NE]
        crossSection = getData( 'cross section for MT = %d' % MT, XSS, offset + 2, NE )
        for i1, energy in enumerate( energyGrid ) : fOut.write( "%20.12e %20.12e\n" % ( energy, crossSection[i1] ) )

    else :
        raise Exception( 'Invalid MFTYPE = %d at location %d' % ( MFTYPE, offset ) )
