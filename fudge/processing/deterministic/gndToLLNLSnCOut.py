# <<BEGIN-copyright>>
# <<END-copyright>>

import math
import datetime
from fudge.core.utilities import fudgeZA
from fudge.core.math import fudge2dGrouping
import fudge
from fudge.processing import processingInfo

__metaclass__ = type

"""
Known bugs:

cout    | data      | sections
--------+-----------+----------
id      | label     | comment
 1      | gp bds    | group boundaries
 2      | speed     | grouped speeds
 3      | gp flux   | grouped fluxes
 4      | gp x-sec  | grouped cross sections by reaction
 5      | en avail  | grouped total available energy
 6      | prod.     | grouped production energy sum over reaction of Q * cross section
 7      | edeptoyo  | grouped outgoing energy to yo
 8      | N/A
 9      | tr corr   | grouped LLNL transport corrections
10      | sigt tr   | grouped total cross section
11      | xfer mat  | grouped transfer matrices for yo
12      | nubrsigf  | grouped fission cross section times nuBar
13      | fiss mat  | grouped neutron fission transfer matrices for neutron
14      |
"""

mass_Energy_speedConstant = 1.38913860e+01 # c4 ((cm/sh)*sqrt(amu/mev)).

#               0        1        2        3        4           5                   6          7         8         9
SnCToLabel = [   "",      "",      "",      "",      "",     "prod",                 "",        "",   "lacs",    "n+i", #  0 -  9
       " elastic",  ",n'",    ",2n",   ",3n",   ",4n",       ",f",              "3np",     "n2p",     "2p",     "pd", # 10 - 19
           ",n'p",  ",pn'",  ",n'd", ",n'da", ",n't",    ",n'he3",             ",n'a",   ",n'2a",  ",n'ta",   ",2np", # 20 - 29
           ",gna", ",2npa",  ",2nd", ",2n'a",  ",npa",     ",dn'",             ",n3a",     ",2a", ",he3,a",    ",tp", # 30 - 39
             ",p",    ",d",    ",t",   ",ta",  ",he3",       ",a",           ",gamma",     ",da",    ",pa",   ",2pa", # 40 - 49
            ",xp",   ",xd",   ",xd", ",xhe3",   ",xa",  ",xgamma",              ",xn",     ",xe",       "",       "", # 50 - 59
               "",      "",      "",      "",      "", "activat.", "yield to (z2,a2)",        "",       "",       "", # 60 - 69
          "total",     "coherent scattering",                 "incoherent scattering",    "photo electric",
                           "pair production",                "pair prod.(electronic)",                  "",
                                          "",                                    "e-",                  "",           # 70 - 79
               "",              "ionization",                        "bremsstrahlung",                  "",
                                          "",                                      "",                  "",
                                          "",                                      "",                  "",           # 80 - 89
               "",  "atomic subshell params",              "subshell transition prob", "whole atom params",
                                          "",                                      "",                  "",
                                          "",                                      "",                  ""  ]         # 90 - 99
SnCToLabel[8] = 'large angle coulomb scat'
SnCToLabel[9] = 'scat+int'

class transferMatrices :

    def __init__( self, lMax, nGroupsProjectile, nGroupsProduct ) :

        self.dataPresent = False
        self.nGroupsProjectile = nGroupsProjectile
        self.nGroupsProduct = nGroupsProduct
        self.TMs = {}
        for l in xrange( lMax ) :
            self.TMs[l] = [ fudge2dGrouping.groupedData( nGroupsProduct * [ 0. ] ) for i in xrange( nGroupsProjectile ) ]

    def __len__( self ) :

        return( len( self.TMs ) )

    def __getitem__( self, l ) :

        return( self.TMs[l] )

    def __iadd__( self, other ) :

        self.dataPresent = True
        lMax = min( len( self.TMs ), len( other ) )
        for l in xrange( lMax ) :
            toTMs = self.TMs[l]
            for i in xrange( len( other[l] ) ) :
                toTMs[i] = toTMs[i] + fudge2dGrouping.groupedData( other[l][i] )
        return( self )

    def hasData( self ) :

        return( self.dataPresent )

def groupSpeeds( groupBoundaries, fluxes, mass, mass_Energy_speedConstant ) :

    fluxGrouped = fudge2dGrouping.groupOneFunction( groupBoundaries, fluxes )
    xs = fudge2dGrouping.commonXIntersection( groupBoundaries.data, fluxes.xArray( ).data )
    if( mass == 0. ) :
        speeds = ( len( groupBoundaries ) - 1 ) * [ 2.99792E+02 ]
    else :
        const = mass_Energy_speedConstant / ( 2. * math.sqrt( mass ) )
        speeds = []
        i = 1
        s = 0.
        x1 = None
        for x2 in xs :
            f2 = fluxes.getValue( x2 )
            sqx2 = math.sqrt( x2 )
            if( x1 != None ) :
                s += ( ( x2 * f1 - x1 * f2 ) + ( f2 - f1 ) * ( x2 + sqx1 * sqx2 + x1 ) / 3. ) / ( sqx2 + sqx1 )
                if( x2 == groupBoundaries[i] ) :
                    if( s != 0. ) : s = fluxGrouped.data[i-1] * const / s
                    speeds.append( s )
                    s = 0.
                    i += 1
            x1 = x2
            sqx1 = sqx2
            f1 = f2
    return( fudge2dGrouping.groupedData( speeds ) )

""" cmattoon, 3/2011: copied here from endl2.py, to cut down on legacy dependencies: """
def ZAToYo( ZA ) :
    """Returns yo designator for ZA.  This mapping is not one-to-one as example
    ZA = 1001 could be yo = 2 or 12.  The smaller yo values is returned."""

    ZA_ = ZA
    if( ZA_ > 10 ) and ( ZA_ < 20 ) : ZA_ -= 10
    if ZA_ in [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1001, 1002, 1003, 2003, 2004 ] :
        if ( ZA_ == 0 ) : return 7
        if ( ZA_ < 1001 ) : return ZA_
        if ( ZA_ == 1001 ) : return 2
        if ( ZA_ == 1002 ) : return 3
        if ( ZA_ == 1003 ) : return 4
        if ( ZA_ == 2003 ) : return 5
        if ( ZA_ == 2004 ) : return 6
    raise Exception( "\nError in ZAToYo: unsupported ZA value = %s" % `ZA` )

def toLLNLSnCOut( self, toOtherData, moreOtherData, first, last, verbosityIndent = '' ) :

    projectileName = self.projectile.getName( )
    targetName = self.target.getName( )
    particles = toOtherData['particles']

    moreOtherData.temperature = None
    moreOtherData.CNumbers = []
    moreOtherData.reactionWithProduct = {}

    EBars = {}
    moreOtherData.nGroups = {}
    for product in particles : 
        moreOtherData.reactionWithProduct[product] = 0
        moreOtherData.nGroups[product] = len( particles[product].groups ) - 1
        E1, EBar = None, []
        for E2 in particles[product].groups :
            if( not ( E1 is None ) ) : EBar.append( 0.5 * ( E1 + E2 ) )
            EBars[product] = EBar
            E1 = E2
    nGroupsProjectile = moreOtherData.nGroups[projectileName]

    productDepositionEnergy = {}
    productTMs = {}
    for product in particles :
        productDepositionEnergy[product] = fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] )
        TM = transferMatrices( particles[product].lMax + 1, nGroupsProjectile, moreOtherData.nGroups[product] )
        TME = transferMatrices( particles[product].lMax + 1, nGroupsProjectile, moreOtherData.nGroups[product] )
        productTMs[product] = { processingInfo.conserveParticle : TM, processingInfo.conserveEnergy : TME }
    moreOtherData.coutData = { 4 : [], 5 : fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] ), 
        6 : fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] ), 7 : productDepositionEnergy, 11 : productTMs,
        12 : fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] ), 
        13 : transferMatrices( 1, nGroupsProjectile, nGroupsProjectile ), 14 : {} }     # 14 data is stored as [ decayTime, { data } ]
                            # where data must have key 'multiplicity' and may have keys 'depositionEnergy' and 'chi'.
    for channel in self : channel.toLLNLSnCOut( toOtherData, moreOtherData, verbosityIndent + '    ' )

    TMs = {}
    iecflg = {}
    energyDepositionCorrection = {}
    for product in particles :
        if( particles[product].conservationFlag == processingInfo.conserveParticle ) :
            TMs[product] = moreOtherData.coutData[11][product][processingInfo.conserveParticle]
            iecflg[product] = 0
        elif( particles[product].conservationFlag == processingInfo.conserveEnergy ) :
            TMs[product] = moreOtherData.coutData[11][product][processingInfo.conserveEnergy]
            EBar = EBars[product]
            for l in xrange( len( TMs[product] ) ) :
                for i, TMRow in enumerate( TMs[product][l] ) :
                    for j, TMCell in enumerate( TMRow ) : TMRow[j] = TMCell / EBar[j]
            iecflg[product] = 1
        elif( particles[product].conservationFlag == processingInfo.conserveParticleAndEnergy ) :
            TMs[product], energyDepositionCorrection[product] = calcConserveParticleAndEnergy( moreOtherData.coutData[11][product], nGroupsProjectile, 
                particles[product].groups )
            iecflg[product] = 3
        else :
            raise Exception( 'Conservation flag for %s = %s not supported' % ( product, particles[product].conservationFlag ) )

    corrections = {}
    for l_ in xrange( len( TMs[projectileName] ) ) :
        l = l_ + 1
        try :
            TM_l = TMs[projectileName][l]
            correction = []
            for TM_lRow in TM_l :
                s = 0
                for d in TM_lRow : s += d
                correction.append( s )
            corrections[l] = fudge2dGrouping.groupedData( correction )
        except :
            corrections[l] = fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] )

    flux = toOtherData['flux']
    fid = 0
    if( flux['name'][:9] == "LLNL_fid_" ) : fid = int( flux['name'][9:] )
    nflx = fid
    
    groups = particles[projectileName].groups
    gidProjectile = 0
    if( groups.name[:9] == "LLNL_gid_" ) : gidProjectile = int( groups.name[9:] )
    productAscStrs = {}
    productDepositionEnergies = moreOtherData.coutData[7]
    productList = { 'n' : 1, 'H1' : 1001, 'H2' : 1002, 'H3' : 1003, 'He3' : 2003, 'He4' : 2004, 'gamma' : 7 }
    nofiss, nout, nppc = 0, 1, 0

    excitationLevel = 0.                    # ???? This line needs work
    targetMass = self.target.getMass( 'amu' )     # Unit 'amu' should not be hard-wired.
    try :
        date = self.attributes['date']
        if( date[:8] == 'YYYYMMDD' ) : date = int( date[11:] )
    except :
        date = 991231                       # bogus date.

    for product in particles :
        ascStr = []
        if( first ) :
            today = datetime.date.today( )
            today = ( 10000 * ( today.year - 2000 ) + 100 * today.month + today.day )
            ascStr += [ '       0    date of ndf file generation:  %.6d' % today, '            date of last physics change:  %.6d' % date ]

        groups = particles[product].groups
        gid = 0
        if( groups.name[:9] == "LLNL_gid_" ) : gid = int( groups.name[9:] )

        ntcflgi = 7
        if( product == projectileName ) : ntcflgi = 5
        ascStr.append( '%6d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d' %
            ( fudgeZA.gndNameToZ_A_Suffix( targetName )[3], particles[product].lMax, nofiss, nout, moreOtherData.nGroups[product], 
                nGroupsProjectile, ntcflgi, moreOtherData.reactionWithProduct[product], nppc, gidProjectile, nflx, iecflg[product], gid, 0, 0 ) )
        ascStr.append( '%2d%2d%11.4E%11.4E%11.4E %.6d unset   ' % ( ZAToYo( productList[projectileName] ), ZAToYo( productList[product] ), 
            excitationLevel, moreOtherData.temperature.getValueAs( 'MeV' ), targetMass, date ) )
        if( product == projectileName ) :
            Cs = [ " %2d" % C for C in moreOtherData.CNumbers ]
            Cs = ''.join( Cs )
            while( len( Cs ) > 0 ) :
                ascStr += [ Cs[:72] ]
                Cs = Cs[72:]
            if( first ) :
                ascStr += groupedDataToCOutStr( groups, '   1 %3d gp bds  ' % gid )
                mass = self.projectile.getMass( 'amu' )
                ascStr += groupedDataToCOutStr( groupSpeeds( groups, flux['data'][0], mass, mass_Energy_speedConstant ), '   2   0 speed   ' )
                groupedFlux = fudge2dGrouping.groupOneFunction( groups, flux['data'][0] )
                ascStr += groupedDataToCOutStr( groupedFlux, '   3 %3d gp flux ' % fid )
            ascStr += [ '   4   0 gp x-sec' ]
            totalCrossSection = fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] )
            for reactionCrossSection in moreOtherData.coutData[4] :
                ascStr.append( '%4d     %-24s' % ( reactionCrossSection.C, SnCToLabel[reactionCrossSection.C] ) )
                ascStr.append( '%3d%11.4E%11.4E%11.4E%11.4E%11.4E' % ( reactionCrossSection.S, reactionCrossSection.Qp, reactionCrossSection.X1p, 
                reactionCrossSection.X2, reactionCrossSection.X3, reactionCrossSection.QEff ) )
                ascStr += groupedDataToCOutStr( reactionCrossSection.data )
                n = len( reactionCrossSection.products )
                if( 'gamma' in reactionCrossSection.products ) : n -= 1
                ascStr.append( '%2d' % n )
                s = ''
                for product2 in [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4' ] :
                    if( product2 in reactionCrossSection.products ) : s += '%6d%6d' % ( productList[product2], reactionCrossSection.products[product2] )
                for product2 in reactionCrossSection.products : 
                    if( product2 not in productList ) :
                        s += '%6d%6d' % ( fudgeZA.gndNameToZ_A_Suffix( product2 )[3], reactionCrossSection.products[product2] )
                        break
                ascStr.append( s )
                totalCrossSection = totalCrossSection + reactionCrossSection.data
            ascStr += groupedDataToCOutStr( moreOtherData.coutData[5] + moreOtherData.coutData[6], '   5   0 en avail' )
            ascStr += groupedDataToCOutStr( moreOtherData.coutData[6], '   6   0  prod.  ' )
        if( product in productDepositionEnergies ) :
            productDepositionEnergy = productDepositionEnergies[product]
            if( productDepositionEnergy.start != productDepositionEnergy.end ) :
                ascStr += groupedDataToCOutStr( productDepositionEnergy, '   7 %3d edeptoyo' % ZAToYo( productList[product] ) )
        if( product == projectileName ) :
            for l in corrections : ascStr += groupedDataToCOutStr( corrections[l], '   9  %2d tr corr ' % l )
            ascStr += groupedDataToCOutStr( totalCrossSection, '  10   0 sigt tr ' )
        if( TMs[product].hasData( ) ) :
            for l in xrange( len( TMs[product] ) ) :
                TM = TMs[product][l]
                ascStr.append( '  11  %2d xfer mat' % l )
                for i in xrange( len( TM ) - 1, -1, -1 ) : ascStr += groupedDataToCOutStr( TM[i] )
            fissionMatrix = moreOtherData.coutData[13]
            if( ( product == 'n' ) and ( fissionMatrix.hasData( ) ) ) :
                ascStr += groupedDataToCOutStr( moreOtherData.coutData[12], label = '  12   0 nubrsigf' )
                ascStr.append( '  13   0 fiss mat' )
                fissionMatrix_0 = fissionMatrix[0]
                for i in xrange( len( fissionMatrix_0 ) - 1, -1, -1 ) : ascStr += groupedDataToCOutStr( fissionMatrix_0[i] )
        delayedFissionData = moreOtherData.coutData[14]
        delayedFissions = delayedFissionData.keys( )
        if( product == 'n' ) :
            if( len( delayedFissions ) > 0 ) :
                delayedFissions.sort( )
                ascStr.append( '  14  %2d delayed fission data' % len( delayedFissions ) )
                for delayedFission in delayedFissions :
                    chiPresent = 0
                    depositionEnergyPresent = 0
                    if( 'depositionEnergy' in delayedFissionData[delayedFission] ) :
                        depositionEnergyPresent = 1
                        depositionEnergy = delayedFissionData[delayedFission]['depositionEnergy']
                    if( 'chi' in delayedFissionData[delayedFission] ) :
                        chi, depEn = delayedFissionData[delayedFission]['chi']
                        chiPresent = 1
                        if( depositionEnergyPresent == 0 ) :
                            depositionEnergy = depEn
                            depositionEnergyPresent = -1
                    ascStr.append( ' %2d %2d %12.5e' % ( chiPresent, depositionEnergyPresent, delayedFission ) )
                    ascStr += groupedDataToCOutStr( delayedFissionData[delayedFission]['multiplicity'] )
                    if( chiPresent != 0 ) : ascStr += groupedDataToCOutStr( chi )
                    if( depositionEnergyPresent != 0 ) : ascStr += groupedDataToCOutStr( depositionEnergy )
        ascStr.append( '  99' )
        if( last ) : ascStr.append( "999999  0  0  0  0  0  0  0  0  0  0                                   *" )
        ascStr.append( '' )
        productAscStrs[product] = '\n'.join( ascStr )
    return( productAscStrs )

def groupedDataToCOutStr( grouped, label = None ) :

    ascStr = []
    if( label is not None ) : ascStr.append( label )
    j, s = 0, ''
    for i in xrange( len( grouped.data ) - 1, -1, -1 ) :
        s += '%12.5E' % grouped.data[i]
        j += 1
        if( j == 6 ) :
            ascStr.append( s )
            j, s = 0, ''
    if( len( s ) ) : ascStr.append( s )
    return( ascStr )

def calcConserveParticleAndEnergy( productTMs, nGroupsProjectile, productGroupBoundaries ) :

    TM1 = productTMs[processingInfo.conserveParticle]
    TME = productTMs[processingInfo.conserveEnergy]

    nGroupsProduct = len( productGroupBoundaries ) - 1
    energyDepositionCorrection = fudge2dGrouping.groupedData( nGroupsProjectile * [ 0. ] )
    energyDepositionCorrection.end = nGroupsProjectile
    EpBar = []
    E1 = None
    for E2 in productGroupBoundaries :      # productGroupBoundaries is endl1dmath which has a good __getitem__.
        if( not( E1 is None ) ) : EpBar.append( 0.5 * ( E1 + E2 ) )
        E1 = E2

    PAndE = []
    if( TM1.hasData( ) ) :                          # Do l = 0 only as per ndfgen
        data = TM1[0]
        dataE = TME[0]
        for iEp, d1 in enumerate( data ) :
            de = dataE[iEp]
            dOut = nGroupsProduct * [ 0. ]
            for i, d in enumerate( d1 ) :
                if( d != 0 ) :
                    energy = de[i] / d
                    f = 1.
                    if( energy < EpBar[i] ) :
                        if( i > 0 ) :
                            f = ( energy - EpBar[i-1] ) / ( EpBar[i] - EpBar[i-1] )
                            dOut[i-1] = dOut[i-1] + ( 1. - f ) * d
                    else :
                        if( i < ( nGroupsProduct - 1 ) ) :
                            f = ( EpBar[i+1] - energy ) / ( EpBar[i+1] - EpBar[i] )
                            dOut[i+1] = dOut[i+1] + ( 1. - f ) * d
                    dOut[i] = dOut[i] + f * d
            PAndE.append( fudge2dGrouping.groupedData( dOut ) )
        for k, TM_1g in enumerate( PAndE ) :                                         # Loop over projectile's group
            dE = TM_1g[-1] * EpBar[-1] - dataE[k][-1]
            if( dE > 0. ) : dE = 0.
            energyDepositionCorrection[k] = dE
    PAndE = { 0 : PAndE }
    for l in xrange( len( TM1 ) ) :                                        # Copy other l order as is.
        if( l == 0 ) : continue
        PAndE[l] = TM1[l]
    TMs = transferMatrices( len( TM1 ), nGroupsProjectile, nGroupsProduct )
    TMs += PAndE

    return( TMs, energyDepositionCorrection )

class LLNLSnCOutCrossSection :

    def __init__( self, reaction ) :

        CS_MT = reaction.getENDL_CS_ENDF_MT( )
        self.C, self.S = CS_MT['C'], CS_MT['S']
        self.Qp = reaction.getQ( 'MeV' )
        self.X1p = 0.
        self.X2 = 0.
        self.X3 = 0.
        self.QEff = self.Qp

        if( self.S == 1 ) : self.QEff = reaction.getQ( 'MeV', final = False )
        if( len( reaction.outputChannel ) == 2 ) : self.X1p = reaction.outputChannel[1].getLevelAsFloat( 'MeV' )
        if( self.S == 8 ) : self.Qp = self.QEff = -self.X1p
        self.data = reaction.crossSection.getFormByToken( fudge.gnd.token.groupedFormToken ).data
        self.products = {}
