# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.utilities import fudgeZA

endfMATBases = { 0 :   1,
                 1 :   1,  2 :   3,  3 :   6,  4 :   9,  5 :  10,  6 :  12,  7 :  14,  8 :  16,  9 :  19, 10 :  -1,
                11 :  23, 12 :  24, 13 :  27, 14 :  28, 15 :  31, 16 :  32, 17 :  35, 18 :  36, 19 :  39, 20 :  40,
                21 :  45, 22 :  46, 23 :  50, 24 :  50, 25 :  55, 26 :  54, 27 :  59, 28 :  58, 29 :  63, 30 :  64,
                31 :  69, 32 :  70, 33 :  75, 34 :  74, 35 :  79, 36 :  78, 37 :  85, 38 :  84, 39 :  89, 40 :  90,
                41 :  93, 42 :  92, 43 :  99, 44 :  96, 45 : 103, 46 : 102, 47 : 107, 48 : 106, 49 : 113, 50 : 112,
                51 : 121, 52 : 120, 53 : 127, 54 : 124, 55 : 133, 56 : 130, 57 : 138, 58 : 136, 59 : 141, 60 : 142,
                61 : 139, 62 : 144, 63 : 151, 64 : 152, 65 : 159, 66 : 156, 67 : 165, 68 : 162, 69 : 169, 70 : 168,
                71 : 175, 72 : 174, 73 : 180, 74 : 180, 75 : 185, 76 : 184, 77 : 191, 78 : 190, 79 : 197, 80 : 196,
                81 : 203, 82 : 204, 83 : 209, 84 : 206, 85 :  -1, 86 :  -1, 87 :  -1, 88 : 223, 89 : 225, 90 : 227,
                91 : 229, 92 : 234, 93 : 230, 94 : 235, 95 : 235, 96 : 240, 97 : 240, 98 : 240, 99 : 265, 100 : 244 }

def endfZAPFromMT( MT ) :
    """This function identifies the outgoing particle (i.e., ZAP) from the MT number. The outgoing 
    particle is almost always a neutron, except in a few cases."""

    if  ( ( ( MT >= 600 ) and ( MT <= 649 ) ) or ( MT == 103 ) ) :   # proton
        return( "H1" )
    elif( ( ( MT >= 650 ) and ( MT <= 699 ) ) or ( MT == 104 ) ) :   # deuteron
        return( "H2" )
    elif( ( ( MT >= 700 ) and ( MT <= 749 ) ) or ( MT == 105 ) ) :   # triton
        return( "H3" )
    elif( ( ( MT >= 750 ) and ( MT <= 799 ) ) or ( MT == 106 ) ) :   # helium-3
        return( "He3" )
    elif( ( ( MT >= 800 ) and ( MT <= 849 ) ) or ( MT == 107 ) ) :   # helium-4
        return( "He4" )
    return( "n" )                                                    # neutron

def ZAAndMATFromParticleName( particleName ) :

    Z, A, suffix, ZA = fudgeZA.gndNameToZ_A_Suffix( particleName )
    m = 0
    if( len( suffix ) ) :
        if( suffix[0] == 'm' ) : m = int( suffix[1:] )
    if( m > 2 ) : raise Exception( 'Unsupport ENDF MAT for particle = %s' % particleName )
    if( A == 0 ) :
        MAT = 100 * Z
    else :
        Zp, AJumps = Z, 3
        if( Z >= 99 ) : Zp, AJumps  = 99, 1
        MAT = 100 * Zp + 25 + AJumps * ( A - endfMATBases[Z] ) + m
        # Kludge for Es254_m1 (MAT logic doesn't allow for isomers above Z=98, so Es254_m1 takes what should
        # be the Es255 MAT):
        if Z==99 and A>=255: MAT += 1
    return( ZA, MAT )

def getParticleNameFromMAT( MAT ):

    Z, MATstuff = divmod( MAT, 100 )
    nStable, nIsomer = divmod( (MATstuff-25), 3 )
    A = endfMATBases[Z] + nStable
    name = fudgeZA.ZAToGNDName( 1000*Z + A )
    if nIsomer: name += '_m%i' % nIsomer
    return name

class endfMTtoC_ProductList :

    def __init__( self, C, reactionLabel, isFission = 0, ns = 0, H1s = 0, H2s = 0, H3s = 0, He3s = 0, He4s = 0, gammas = 0, residualLevel = None ) :

        self.C = C
        self.residualLevel = residualLevel
        self.reactionLabel = reactionLabel
        self.isFission = isFission
        self.productCounts = { 'n' : ns, 'H1' : H1s, 'H2' : H2s, 'H3' : H3s, 'He3' : He3s, 'He4' : He4s, 'gamma' : gammas }

    def __getitem__( self, product ) :

        return( self.productCounts[product] )

    def __repr__( self ) :

        s = ''
        for p in [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4', 'gamma' ] :
            if( self.productCounts[p] != 0 ) : s += " %5s = %d:" % ( p, self.productCounts[p] )
        s = "C = %s: isFission = %5s:%s --- %s" % ( self.C, self.isFission != 0, s, self.reactionLabel )
        return( s )

def endfMTtoC_ProductList_excitedStateInitializer( list, MTGround, MTContinuum, C, label, ns = 0, H1s = 0, H2s = 0, H3s = 0, He3s = 0, He4s = 0, gammas = 0 ) :

    levelSuffixes = [ "", "st", "nd", "rd" ]
    list[MTGround] = endfMTtoC_ProductList( C, "(z,%s[0]) -- to ground state" % label, 0, ns, H1s, H2s, H3s, He3s, He4s, gammas, 0 )
    for i in xrange( MTGround + 1, MTContinuum ) :
        level = i - MTGround 
        try :
            levelSuffix = levelSuffixes[level]
        except :
            levelSuffix = "th"
        list[i] = endfMTtoC_ProductList( C, "(z,%s[%d]) -- to %d%s excited state" % ( label, level, level, levelSuffix ), \
                0, ns, H1s, H2s, H3s, He3s, He4s, gammas, level )
    list[MTContinuum] = endfMTtoC_ProductList( C, "(z,%s[c]) -- excitation to continuum" % label, 0, ns, H1s, H2s, H3s, He3s, He4s, gammas, 'c' )

endfMTtoC_ProductLists = {}
endfMTtoC_ProductLists[1]   = endfMTtoC_ProductList(  1, "(n,total)" )
endfMTtoC_ProductLists[2]   = endfMTtoC_ProductList( 10, "(z,elastic)" )
endfMTtoC_ProductLists[3]   = endfMTtoC_ProductList( 55, "(z,non-elastic)" )
endfMTtoC_ProductLists[4]   = endfMTtoC_ProductList( 55, "(z,n), redundant", 0, 1, 0, 0, 0, 0, 0, -1 )
endfMTtoC_ProductLists[5]   = endfMTtoC_ProductList(  5, "(z,anything)" )
endfMTtoC_ProductLists[10]  = endfMTtoC_ProductList( -1, "(z,continuum)" )
endfMTtoC_ProductLists[11]  = endfMTtoC_ProductList( 32, "(n,2nd)", 0, 2, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[16]  = endfMTtoC_ProductList( 12, "(z,2n)", 0, 2, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[17]  = endfMTtoC_ProductList( 13, "(z,3n)", 0, 3, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[18]  = endfMTtoC_ProductList( 15, "(z,f)", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[19]  = endfMTtoC_ProductList( -1, "(n,f) -- 1st chance fission.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[20]  = endfMTtoC_ProductList( -1, "(n,nf) -- 2nd chance fission.", -1, 1, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[21]  = endfMTtoC_ProductList( -1, "(n,2nf) -- 3rd chance fission.", -1, 2, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[22]  = endfMTtoC_ProductList( 26, "(z,na)", 0, 1, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[23]  = endfMTtoC_ProductList( 36, "(z,n3a)", 0, 1, 0, 0, 0, 0, 3, 0 )
endfMTtoC_ProductLists[24]  = endfMTtoC_ProductList( 33, "(z,2na)", 0, 2, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[25]  = endfMTtoC_ProductList( 16, "(z,3na)", 0, 3, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[27]  = endfMTtoC_ProductList( -1, "(n,abs)" )
endfMTtoC_ProductLists[28]  = endfMTtoC_ProductList( 20, "(z,np)", 0, 1, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[29]  = endfMTtoC_ProductList( 27, "(z,n2a)", 0, 1, 0, 0, 0, 0, 2, 0 )
endfMTtoC_ProductLists[30]  = endfMTtoC_ProductList( -1, "(z,2n2a)", 0, 2, 0, 0, 0, 0, 2, 0 )
endfMTtoC_ProductLists[32]  = endfMTtoC_ProductList( 22, "(z,nd)", 0, 1, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[33]  = endfMTtoC_ProductList( 24, "(z,nt)", 0, 1, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[34]  = endfMTtoC_ProductList( 25, "(z,nH)", 0, 1, 0, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[35]  = endfMTtoC_ProductList( -1, "(z,nd2a)", 0, 1, 0, 1, 0, 0, 2, 0 )
endfMTtoC_ProductLists[36]  = endfMTtoC_ProductList( -1, "(z,nt2a)", 0, 1, 0, 0, 1, 0, 2, 0 )
endfMTtoC_ProductLists[37]  = endfMTtoC_ProductList( 14, "(z,4n)", 0, 4, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[38]  = endfMTtoC_ProductList( -1, "(n,3nf) -- 4th chance fission.", -1, 3, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[41]  = endfMTtoC_ProductList( 29, "(z,2np)", 0, 2, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[42]  = endfMTtoC_ProductList( 16, "(z,3np)", 0, 3, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[44]  = endfMTtoC_ProductList( 17, "(n,n2p)", 0, 1, 2, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[45]  = endfMTtoC_ProductList( 34, "(n,npa)", 0, 1, 1, 0, 0, 0, 1, 0 )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 50, 91, 11, "n", ns = 1 )
endfMTtoC_ProductLists[101] = endfMTtoC_ProductList( -1, "(n,disappearance)" )
endfMTtoC_ProductLists[102] = endfMTtoC_ProductList( 46, "(z,g)", 0, 0, 0, 0, 0, 0, 0, 1 )
endfMTtoC_ProductLists[103] = endfMTtoC_ProductList( 40, "(z,p)", 0, 0, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[104] = endfMTtoC_ProductList( 41, "(z,d)", 0, 0, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[105] = endfMTtoC_ProductList( 42, "(z,t)", 0, 0, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[106] = endfMTtoC_ProductList( 44, "(z,H)", 0, 0, 0, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[107] = endfMTtoC_ProductList( 45, "(z,a)", 0, 0, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[108] = endfMTtoC_ProductList( 37, "(z,2a)", 0, 0, 0, 0, 0, 0, 2, 0 )
endfMTtoC_ProductLists[109] = endfMTtoC_ProductList( -1, "(z,3a)", 0, 0, 0, 0, 0, 0, 3, 0 )
endfMTtoC_ProductLists[111] = endfMTtoC_ProductList( 18, "(z,2p)", 0, 0, 2, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[112] = endfMTtoC_ProductList( 48, "(z,pa)", 0, 0, 1, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[113] = endfMTtoC_ProductList( 42, "(z,t2a)", 0, 0, 0, 0, 1, 0, 2, 0 )
endfMTtoC_ProductLists[114] = endfMTtoC_ProductList( -1, "(z,d2a)", 0, 0, 0, 1, 0, 0, 2, 0 )
endfMTtoC_ProductLists[115] = endfMTtoC_ProductList( 19, "(z,pd)", 0, 0, 1, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[116] = endfMTtoC_ProductList( 39, "(z,pt)", 0, 0, 1, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[117] = endfMTtoC_ProductList( 47, "(z,da)", 0, 0, 0, 1, 0, 0, 1, 0 )
endfMTtoC_ProductLists[151] = endfMTtoC_ProductList( -1, "(n,resonance)" )
# cmattoon Septmeber 2011, additional MT #s defined by CSEWG in 2010:
endfMTtoC_ProductLists[152] = endfMTtoC_ProductList( -1, "(z,5n)", 0, 5, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[153] = endfMTtoC_ProductList( -1, "(z,6n)", 0, 6, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[154] = endfMTtoC_ProductList( -1, "(z,2nt)", 0, 2, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[155] = endfMTtoC_ProductList( 43, "(z,ta)", 0, 0, 0, 0, 1, 0, 1, 0 )
endfMTtoC_ProductLists[156] = endfMTtoC_ProductList( -1, "(z,4np)", 0, 4, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[157] = endfMTtoC_ProductList( -1, "(z,3nd)", 0, 3, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[158] = endfMTtoC_ProductList( 23, "(z,nda)", 0, 1, 0, 1, 0, 0, 1, 0 )
endfMTtoC_ProductLists[159] = endfMTtoC_ProductList( 31, "(z,2npa)", 0, 2, 1, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[160] = endfMTtoC_ProductList( -1, "(z,7n)", 0, 7, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[161] = endfMTtoC_ProductList( -1, "(z,8n)", 0, 8, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[162] = endfMTtoC_ProductList( -1, "(z,5np)", 0, 5, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[163] = endfMTtoC_ProductList( -1, "(z,6np)", 0, 6, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[164] = endfMTtoC_ProductList( -1, "(z,7np)", 0, 7, 1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[165] = endfMTtoC_ProductList( -1, "(z,4na)", 0, 4, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[166] = endfMTtoC_ProductList( -1, "(z,5na)", 0, 5, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[167] = endfMTtoC_ProductList( -1, "(z,6na)", 0, 6, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[168] = endfMTtoC_ProductList( -1, "(z,7na)", 0, 7, 0, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[169] = endfMTtoC_ProductList( -1, "(z,4nd)", 0, 4, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[170] = endfMTtoC_ProductList( -1, "(z,5nd)", 0, 5, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[171] = endfMTtoC_ProductList( -1, "(z,6nd)", 0, 6, 0, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[172] = endfMTtoC_ProductList( -1, "(z,3nt)", 0, 3, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[173] = endfMTtoC_ProductList( -1, "(z,4nt)", 0, 4, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[174] = endfMTtoC_ProductList( -1, "(z,5nt)", 0, 5, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[175] = endfMTtoC_ProductList( -1, "(z,6nt)", 0, 6, 0, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[176] = endfMTtoC_ProductList( -1, "(z,2nH)", 0, 2, 0, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[177] = endfMTtoC_ProductList( -1, "(z,3nH)", 0, 3, 0, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[178] = endfMTtoC_ProductList( -1, "(z,4nH)", 0, 4, 0, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[179] = endfMTtoC_ProductList( -1, "(z,3n2p)", 0, 3, 2, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[180] = endfMTtoC_ProductList( -1, "(z,3n2a)", 0, 3, 0, 0, 0, 0, 2, 0 )
endfMTtoC_ProductLists[181] = endfMTtoC_ProductList( -1, "(z,3npa)", 0, 3, 1, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[182] = endfMTtoC_ProductList( -1, "(z,dt)", 0, 0, 0, 1, 1, 0, 0, 0 )
endfMTtoC_ProductLists[183] = endfMTtoC_ProductList( -1, "(z,npd)", 0, 1, 1, 1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[184] = endfMTtoC_ProductList( -1, "(z,npt)", 0, 1, 1, 0, 1, 0, 0, 0 )
endfMTtoC_ProductLists[185] = endfMTtoC_ProductList( -1, "(z,ndt)", 0, 1, 0, 1, 1, 0, 0, 0 )
endfMTtoC_ProductLists[186] = endfMTtoC_ProductList( -1, "(z,npH)", 0, 1, 1, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[187] = endfMTtoC_ProductList( -1, "(z,ndH)", 0, 1, 0, 1, 0, 1, 0, 0 )
endfMTtoC_ProductLists[188] = endfMTtoC_ProductList( -1, "(z,ntH)", 0, 1, 0, 0, 1, 1, 0, 0 )
endfMTtoC_ProductLists[189] = endfMTtoC_ProductList( 28, "(z,nta)", 0, 1, 0, 0, 1, 0, 1, 0 )
endfMTtoC_ProductLists[190] = endfMTtoC_ProductList( -1, "(z,2n2p)", 0, 2, 2, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[191] = endfMTtoC_ProductList( -1, "(z,pH)", 0, 0, 1, 0, 0, 1, 0, 0 )
endfMTtoC_ProductLists[192] = endfMTtoC_ProductList( -1, "(z,dH)", 0, 0, 0, 1, 0, 1, 0, 0 )
endfMTtoC_ProductLists[193] = endfMTtoC_ProductList( 38, "(z,Ha)", 0, 0, 0, 0, 0, 1, 1, 0 )
endfMTtoC_ProductLists[194] = endfMTtoC_ProductList( -1, "(z,4n2p)", 0, 4, 2, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[195] = endfMTtoC_ProductList( -1, "(z,4n2a)", 0, 4, 0, 0, 0, 0, 2, 0 )
endfMTtoC_ProductLists[196] = endfMTtoC_ProductList( -1, "(z,4npa)", 0, 4, 1, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[197] = endfMTtoC_ProductList( -1, "(z,3p)", 0, 0, 3, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[198] = endfMTtoC_ProductList( -1, "(z,n3p)", 0, 1, 3, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[199] = endfMTtoC_ProductList( -1, "(z,3n2pa)", 0, 3, 2, 0, 0, 0, 1, 0 )
endfMTtoC_ProductLists[200] = endfMTtoC_ProductList( -1, "(z,5n2p)", 0, 5, 2, 0, 0, 0, 0, 0 )
# end of new MT #s
endfMTtoC_ProductLists[201] = endfMTtoC_ProductList( -1, "(z,Xn)", 0, -1, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[202] = endfMTtoC_ProductList( -1, "(z,Xg)", 0, 0, 0, 0, 0, 0, 0, -1 )
endfMTtoC_ProductLists[203] = endfMTtoC_ProductList( 50, "(z,Xp)", 0, 0, -1, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[204] = endfMTtoC_ProductList( 51, "(z,Xd)", 0, 0, 0, -1, 0, 0, 0, 0 )
endfMTtoC_ProductLists[205] = endfMTtoC_ProductList( 52, "(z,Xt)", 0, 0, 0, 0, -1, 0, 0, 0 )
endfMTtoC_ProductLists[206] = endfMTtoC_ProductList( 53, "(z,XH)", 0, 0, 0, 0, 0, -1, 0, 0 )
endfMTtoC_ProductLists[207] = endfMTtoC_ProductList( 54, "(z,Xa)", 0, 0, 0, 0, 0, 0, -1, 0 )
endfMTtoC_ProductLists[208] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[209] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[210] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[211] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[212] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[213] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[214] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[215] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[216] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[217] = endfMTtoC_ProductList( -1, "Various meson and antiparticle production $sigma$'s" )
endfMTtoC_ProductLists[251] = endfMTtoC_ProductList( -1, "Various elastic neutrons scattering parameters." )
endfMTtoC_ProductLists[252] = endfMTtoC_ProductList( -1, "Various elastic neutrons scattering parameters." )
endfMTtoC_ProductLists[253] = endfMTtoC_ProductList( -1, "Various elastic neutrons scattering parameters." )
endfMTtoC_ProductLists[301] = endfMTtoC_ProductList( -1, "Energy release for total and partial $sigma$'s." )
endfMTtoC_ProductLists[451] = endfMTtoC_ProductList( -1, "Heading or title information, MF=1 only." )
endfMTtoC_ProductLists[452] = endfMTtoC_ProductList( 15, "(z,f) $bar{nu}$ total, i.e. prompt plus delayed, fission.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[454] = endfMTtoC_ProductList( 15, "(z,f) Independent fission product yields.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[455] = endfMTtoC_ProductList( 15, "(z,f) $bar{nu}$ for delayed fission.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[456] = endfMTtoC_ProductList( 15, "(z,f) $bar{nu}$ for prompt fission.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[457] = endfMTtoC_ProductList( -1, "(z,f) Radioactive decay data.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[458] = endfMTtoC_ProductList( 15, "(z,f) Energy release in fission for incident $n$'s.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[459] = endfMTtoC_ProductList( 15, "(z,f) Cumulative fission product yields.", -1, 0, 0, 0, 0, 0, 0, 0 )
endfMTtoC_ProductLists[500] = endfMTtoC_ProductList( -1, "Total charged particle stopping power." )
endfMTtoC_ProductLists[501] = endfMTtoC_ProductList( 70, "Total photon interaction $sigma$." )
endfMTtoC_ProductLists[502] = endfMTtoC_ProductList( 71, "Photon coherent scattering." )
endfMTtoC_ProductLists[504] = endfMTtoC_ProductList( 72, "Photon incoherent scattering." )
endfMTtoC_ProductLists[505] = endfMTtoC_ProductList( -1, "Imaginary scattering factor." )
endfMTtoC_ProductLists[506] = endfMTtoC_ProductList( -1, "Real scattering factor." )
endfMTtoC_ProductLists[515] = endfMTtoC_ProductList( -1, "Pair production, electron field." )
endfMTtoC_ProductLists[516] = endfMTtoC_ProductList( 74, "Pair production." )
endfMTtoC_ProductLists[517] = endfMTtoC_ProductList( -1, "Pair production, nuclear field." )
endfMTtoC_ProductLists[522] = endfMTtoC_ProductList( 73, "Photoelectric absorption." )
endfMTtoC_ProductLists[534] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[535] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[536] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[537] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[538] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[539] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[540] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[541] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[542] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[543] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[544] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[545] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[546] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[547] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[548] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[549] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[550] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[551] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[552] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[553] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[554] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[555] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[556] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[557] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[558] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[559] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[560] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[561] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[562] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[563] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[564] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[565] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[566] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[567] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[568] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[569] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[570] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductLists[571] = endfMTtoC_ProductList( -1, "Various subshell photoelectric $sigma$'s." )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 600, 649, 40, "p", H1s = 1 )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 650, 699, 41, "d", H2s = 1 )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 700, 749, 42, "t", H3s = 1 )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 750, 799, 44, "He3", He3s = 1 )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 800, 849, 45, "a", He4s = 1 )
endfMTtoC_ProductLists[851] = endfMTtoC_ProductList( -1, "Lumped reaction covariances." )
endfMTtoC_ProductList_excitedStateInitializer( endfMTtoC_ProductLists, 875, 891, 12, "2n", ns = 2 )

def getCSFromMT( MT ) :

    if( MT ==   1 ) : return(   1, 0 )     # (z,total)
    if( MT ==   2 ) : return(  10, 0 )     # (z,elas)
    if( MT in [ 3, 4, 5, 10, 25, 27, 30, 35, 36, 101, 109, 113, 114,
        152, 153, 154, 156, 157, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
        172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188,
        190, 191, 192, 194, 195, 196, 197, 198, 199, 200] ) : return(   -MT, 0 )
    if( MT ==  11 ) : return(  32, 0 )
    if( MT ==  16 ) : return(  12, 0 )
    if( MT ==  17 ) : return(  13, 0 )
    if( MT ==  18 ) : return(  15, 0 )
    if( MT ==  19 ) : return(  15, 0 )
    if( MT ==  20 ) : return(  15, 0 )
    if( MT ==  21 ) : return(  15, 0 )
    if( MT ==  22 ) : return(  26, 0 )
    if( MT ==  23 ) : return(  36, 0 )
    if( MT ==  24 ) : return(  33, 0 )
    if( MT ==  28 ) : return(  20, 0 )
    if( MT ==  29 ) : return(  27, 0 )
    if( MT ==  32 ) : return(  22, 0 )
    if( MT ==  33 ) : return(  24, 0 )
    if( MT ==  34 ) : return(  25, 0 )
    if( MT ==  37 ) : return(  14, 0 )
    if( MT ==  38 ) : return(  15, 0 )
    if( MT ==  41 ) : return(  29, 0 )
    if( MT ==  42 ) : return(  16, 0 )
    if( MT ==  44 ) : return(  17, 0 )
    if( MT ==  45 ) : return(  34, 0 )
    if( 50 <= MT < 91 ) : return( 11, 1 ) 
    if( MT ==  91 ) : return(  11, 0 )
    if( MT == 102 ) : return(  46, 0 )
    if( MT == 103 ) : return(  40, 0 )
    if( MT == 104 ) : return(  41, 0 )
    if( MT == 105 ) : return(  42, 0 )
    if( MT == 106 ) : return(  44, 0 )
    if( MT == 107 ) : return(  45, 0 )
    if( MT == 108 ) : return(  37, 0 )
    if( MT == 111 ) : return(  18, 0 )
    if( MT == 112 ) : return(  48, 0 )
    if( MT == 115 ) : return(  19, 0 )
    if( MT == 116 ) : return(  39, 0 )
    if( MT == 117 ) : return(  47, 0 )
    if( MT == 155 ) : return(  43, 0 )
    if( MT == 158 ) : return(  23, 0 )
    if( MT == 159 ) : return(  31, 0 )
    if( MT == 189 ) : return(  28, 0 )
    if( MT == 193 ) : return(  38, 0 )
    if( MT == 452 ) : return( 15, 0 )    # prompt plus delayed fission neutrons
    if( MT == 455 ) : return( 15, 7 )    # delayed fission neutrons
    if( MT == 456 ) : return( 15, 0 )    # prompt fission neutrons
    if( MT == 458 ) : return( 15, 0 )    # prompt fission neutron energy
    if( 600 <= MT < 649 ) : return( 40, 1 ) 
    if( MT == 649 ) : return( 40, 0 ) 
    if( 650 <= MT < 699 ) : return( 41, 1 ) 
    if( MT == 699 ) : return( 41, 0 ) 
    if( 700 <= MT < 749 ) : return( 42, 1 ) 
    if( MT == 749 ) : return( 42, 0 ) 
    if( 750 <= MT < 799 ) : return( 44, 1 ) 
    if( MT == 799 ) : return( 44, 0 ) 
    if( 800 <= MT < 849 ) : return( 45, 1 ) 
    if( MT == 849 ) : return( 45, 0 ) 
    if( 875 <= MT < 891 ) : return( 12, 1 ) 
    if( MT == 891 ) : return( 12, 0 ) 
    raise Exception( 'MT = %d is not supported for conversion to C, S' % MT )

def getMTFromC( C ) :

    if( C ==  1 ) : return(   1 )   # (z,total)
    if( C ==  5 ) : return(  -5 )   # (z,prod)
    if( C ==  8 ) : return(  -8 )   # (z,lacs)
    if( C ==  9 ) : return(  -9 )   # (z,n+i)
    if( C == 10 ) : return(   2 )   # (z,elas)
    if( C == 11 ) : return(  50 )   # (z,n)
    if( C == 12 ) : return(  16 )   # (z,2n)
    if( C == 13 ) : return(  17 )   # (z,3n)
    if( C == 14 ) : return(  37 )   # (z,4n)
    if( C == 15 ) : return(  18 )   # (z,f)

    if( C == 16 ) : return(  42 )   # (z,3np)
    if( C == 17 ) : return(  44 )   # (z,n2p)
    if( C == 18 ) : return( 111 )   # (z,2p)
    if( C == 19 ) : return( 115 )   # (z,pd)

    if( C == 20 ) : return(  28 )   # (z,np)
    if( C == 21 ) : return( -20 )   # (z,pn)
    if( C == 22 ) : return(  32 )   # (z,nd)
    if( C == 23 ) : return( 158 )   # (z,nda)
    if( C == 24 ) : return(  33 )   # (z,nt)
    if( C == 25 ) : return(  34 )   # (z,nHe3)
    if( C == 26 ) : return(  22 )   # (z,na)
    if( C == 27 ) : return(  29 )   # (z,n2a)
    if( C == 28 ) : return( 189 )   # (z,nta)
    if( C == 29 ) : return(  41 )   # (z,2np)

    if( C == 30 ) : return( -30 )   # (z,gna)
    if( C == 31 ) : return( 159 )   # (z,2npa)
    if( C == 32 ) : return(  11 )   # (z,2nd)
    if( C == 33 ) : return(  24 )   # (z,2na)
    if( C == 34 ) : return(  45 )   # (z,npa)
    if( C == 35 ) : return(  32 )   # (z,dn), ENDF does not have an (z,dn) reaction only an (z,nd) reaction and ENDL's only (z,dn) reaction is not two-body so order does not matter.
    if( C == 36 ) : return(  23 )   # (z,n3a)
    if( C == 37 ) : return( 108 )   # (z,2a)
    if( C == 38 ) : return( 193 )   # (z,He3 a)
    if( C == 39 ) : return( 116 )   # (z,pt)

    if( C == 40 ) : return( 600 )   # (z,p)
    if( C == 41 ) : return( 650 )   # (z,d)
    if( C == 42 ) : return( 700 )   # (z,t)
    if( C == 43 ) : return( 155 )   # (z,ta)
    if( C == 44 ) : return( 750 )   # (z,He3)
    if( C == 45 ) : return( 800 )   # (z,a)
    if( C == 46 ) : return( 102 )   # (z,g)
    if( C == 47 ) : return( 117 )   # (z,da)
    if( C == 48 ) : return( 112 )   # (z,pa)
    if( C == 49 ) : return( -49 )   # (z,2pa)

    if( C == 50 ) : return( -50 )   # (z,Xp)
    if( C == 51 ) : return( -51 )   # (z,Xd)
    if( C == 52 ) : return( -52 )   # (z,Xt)
    if( C == 53 ) : return( -53 )   # (z,XHe3)
    if( C == 54 ) : return( -54 )   # (z,Xa)
    if( C == 55 ) : return( -55 )   # (z,Xg)
    if( C == 56 ) : return( -56 )   # (z,Xn)
    if( C == 57 ) : return( -57 )   # (z,Xe)

    if( C == 70 ) : return( 501 )   # (z,totp)
    if( C == 71 ) : return( 502 )   # (z,coh)
    if( C == 72 ) : return( 504 )   # (z,incoh)
    if( C == 73 ) : return( 522 )   # (z,photo)
    if( C == 74 ) : return( 516 )   # (z,pair)
    if( C == 75 ) : return( -75 )   # (z,triplet)
    if( C == 78 ) : return( -78 )   # (z,ic)

    if( C == 81 ) : return( -81 )   # (z,ion)
    if( C == 82 ) : return( -82 )   # (z,brem)
    if( C == 83 ) : return( -83 )   # (z,excit)
    if( C == 84 ) : return( -84 )   # (z,coll)

    if( C == 91 ) : return( -91 )   # (z,shell)
    if( C == 92 ) : return( -92 )   # (z,trans)
    if( C == 93 ) : return( -93 )   # (z,whole)
    raise Exception( 'C = %d is not supported for conversion to MT' % C )

class ENDLCS_To_ENDFMT :

    def __init__( self, projectile ) :

        self.projectile = projectile
        self.MTPrimes = {}

    def getMTFromCS( self, C, S, CCounts = 0 ) :

        def MTPrimes( self, MT, S, projectile ) :

            if( S == 0 ) :
                if( MT == 50 ) :
                    MT = 91
                else :
                    MT += 49
            else :
                if( MT not in self.MTPrimes ) :
                    self.MTPrimes[MT] = -1
                    if( self.projectile == projectile ) : self.MTPrimes[MT] += 1
                self.MTPrimes[MT] += 1
                MT += self.MTPrimes[MT]
            return( MT )

        MT = getMTFromC( C )
        if(   MT ==  50 ) :
            MT = MTPrimes( self, MT, S, 'n' )
        elif( MT == 600 ) :
            MT = MTPrimes( self, MT, S, 'H1' )
            if( ( S == 0 ) and ( CCounts == 1 ) ) : MT = 103
        elif( MT == 650 ) :
            MT = MTPrimes( self, MT, S, 'H2' )
            if( ( S == 0 ) and ( CCounts == 1 ) ) : MT = 104
        elif( MT == 700 ) :
            MT = MTPrimes( self, MT, S, 'H3' )
            if( ( S == 0 ) and ( CCounts == 1 ) ) : MT = 105
        elif( MT == 750 ) :
            MT = MTPrimes( self, MT, S, 'He3' )
            if( ( S == 0 ) and ( CCounts == 1 ) ) : MT = 106
        elif( MT == 800 ) :
            MT = MTPrimes( self, MT, S, 'He4' )
            if( ( S == 0 ) and ( CCounts == 1 ) ) : MT = 107
        return( MT )

def ENDF_MTZAEquation( projectileZA, targetZA, MT ) :
    """This function returns a python list of length 2. The first element is a list of all outgoing particle ZA's
    (include the residual) for the reaction of projectileZA + targetZA with ENDF's reaction identifier MT. 
    The second element is a reaction equation for this projectileZA, targetZA and MT. For example
    ENDF_MTZAEquation( 1, 95242, 22 ) returns
        ([1, 2004, 93238], 'n + Am242 -> n + He4 + Np238')
    That is, for a neutron ( projectileZA = 1 ) hitting Am242 ( targetZA = 95242 ) with MT = 22 - ENDF (z,na) reaction - 
    the outgoing particle ZA's are [1, 2004, 93238] and the reaction equation is 'n + Am242 -> n + He4 + Np238'."""

    if( ( MT < 0 ) or ( MT > 999 ) or ( MT in [ 1, 3, 5, 10, 18, 19, 20, 21, 27, 38, 101, 151 ] or ( 200 < MT < 600 ) or ( 850 < MT < 875 ) ) ) :
        raise Exception( 'MT = %s is no supported' % MT )
    elif( MT == 2 ) :
        productCounts = { fudgeZA.ZAToGNDName( projectileZA ) : 1 }
        level = None
    elif( MT == 4 ) :
        productCounts = { fudgeZA.ZAToGNDName( projectileZA ) : 1 }
        level = None
    else :
        productCounts = endfMTtoC_ProductLists[MT].productCounts
        level = endfMTtoC_ProductLists[MT].residualLevel
    compoundZA = projectileZA + targetZA
    residualZA = compoundZA
    productCountList = []
    adder, equationZA, equation = '', [], '%s + %s ->' % ( fudgeZA.ZAToGNDName( projectileZA ), fudgeZA.ZAToGNDName( targetZA ) )
    for product in productCounts :
        if( product == 'gamma' ) :
            productZA = 0
        else :
            productZA = fudgeZA.gndNameToZ_A_Suffix( product )[-1]
        if( productCounts[product] > 0 ) : productCountList.append( [ productZA, product, productCounts[product] ] )
    productCountList.sort( )
    for productZA, token, count in productCountList :
        residualZA -= count * productZA
        for i in xrange( count ) :
            equation += ' %s%s' % ( adder, token )
            equationZA.append( productZA )
            adder = '+ '
    levelStr = ''
    if( not( level is None ) ) :
        if( type( level ) == type( 0 ) ) :
            if( level > 0 ) : levelStr = "_e%s" % level
        elif( type( level ) == type( '' ) ) :
            if( level == 'c' ) :
                levelStr = '_c'
            else :
                raise Exception( 'Unknown level specifier = %s' % level )
        else :
            raise Exception( 'Unknown level specifier = %s' % level )
    equation += ' %s%s%s' % ( adder, fudgeZA.ZAToGNDName( residualZA ), levelStr )
    equationZA.append( residualZA )
    return( equationZA, equation )
