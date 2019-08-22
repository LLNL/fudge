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

"""
This module contains miscellaneous functions.
"""

import math
from fudge.core.utilities import brb
from fudge.particles import nuclear
import endl_C
import endl_I
import endl_Z
import endl_y

def CInfo( w = 4 ) :
    """Prints a list of C and brief description for supported C-values."""

    n = endl_C.endl_nCMax( )
    l = endl_C.endl_CLabelMaxWidth( )
    if( l < 4 ) : l = 4
    Fmt = "%%4d %%%ds" % l
    Cs = []
    i = 1
    while( i <= n ) :
        s = endl_C.endl_CLabel( i )
        if( s == None ) : s = "----"
        Cs.append( Fmt % ( i, s ) )
        i += 1
    brb.tylist( Cs, w = w, sep = "  |" )

def IInfo( ) :
    """Prints a list of I, number of columns and brief description for supported I-values."""

    I = 0
    print "   I cols  comment"
    print "-----------------"
    while( I < 100 ) :
        s = endl_I.endl_ILabel( I )
        if( s != None ) :
            print "%4d %5d  %s" % ( I, endl_I.endl_IColumns( I ), s )
        I += 1

def yInfo( ) :
    """Prints a list of y-id, symbol and name for supported projectiles."""

    ls = endl_y.endl_yLabelMaxWidth( )
    ll = endl_y.endl_yLongLabelMaxWidth( )
    Fmt = "   y %%%ds %%%ds" % ( ls, ll )
    print Fmt % ( "Symbol", "Name" )
    dashes = ( 6 + ls + ll ) * "-"
    print dashes
    Fmt = "%%4d %%%ds %%%ds" % ( ls, ll )
    y = 1
    while( 1 ) :
        ls = endl_y.endl_yLabel( y )
        if( ls == None ) : break
        ll = endl_y.endl_yLongLabel( y )
        print Fmt % ( y, ls, ll )
        y += 1

def ZInfo( w = 2 ) :
    """Prints a list of Z, symbol and name for most elements."""

    n = endl_Z.endl_nZs( )
    Zs = []
    Z = 0
    while( Z < n ) :
        Zs.append( "%4d %4s %15s" % ( Z, endl_Z.endl_ZSymbol( Z ), endl_Z.endl_ZLabel( Z ) ) )
        Z += 1
    brb.tylist( Zs, w = w, sep = "  |" )

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

def yoToZA( yo_i ) :
    """Returns ZA designator, or special code, from yo id.  That is, ::
      yo_i = 0, 1,    2,    3,    4,    5,    6, 7,  8,  9,  10, 11,   12,   13,   14,   15,   16
      is converted to
             0, 1, 1001, 1002, 1003, 2003, 2004, 0, -8, -9, -10,  1, 1001, 1002, 1003, 2003, 2004
    with the following meanings ::
          none, n,    p,    d,    t,  He3,   He, g, b+, b-,  EC,  n,    p,    d,    t,  He3,   He.
    yo_i must be a python integer."""

    if ( yo_i < 0 ) or ( yo_i > 16 ) : raise Exception( "\nError in yoToZA: unsupported yo value = %s" % `yo_i` )
    yo = yo_i
    if ( yo > 10 ) : yo -= 10
    return ( 0, 1, 1001, 1002, 1003, 2003, 2004, 0, -8, -9, -10 )[yo]

def yoToZA_forNuclei( yo ) :
    """Same as yoToZA except returns 0 if yo is in [ 7, 8, 9, 17, 18, 19 ] (i.e., is a gamma, positron or electron)."""

    if( yo in [ 7, 8, 9, 17, 18, 19 ] ) : return( 0 )
    return( yoToZA( yo ) )

def ZandAFromZA( ZA ) :
    """Returns the tuple ( Z, A ) from a ZA designator.  ZA must be a python integer."""

    import endlmisc
    ZA_ = endlmisc.intZA( ZA )
    Z = ZA_ / 1000
    A = ZA_ - 1000 * Z
    return ( Z, A )

def AFromZA( ZA ) :
    """Returns A from ZA designator ( ZA = 1000 * Z + A ).  ZA must be a python integer."""

    return ZandAFromZA( ZA )[1]

def ZFromZA( ZA ) :
    """ZFromZA( ZA )
    Returns Z from ZA designator ( ZA = 1000 * Z + A ).  ZA must be a python integer."""

    return ZandAFromZA( ZA )[0]

def residualZA_yos_Q( yi_, ZA_, C, targetELevel = 0., X4 = 0, specialCases = 1, printQWarning = True, bdflsFile = None ) :
    """Returns a tuple containing the residual ZA', yos and Q for reaction of type 
    C involving incident particle yi and target ZA.  yos is a list.  For example,
    the reation, ZA + yi -> ZA' + yo_1 + yo_2 + ... + yo_n + Q,
    would return ( ZA', [ yo_1, yo_2, ..., yo_n ], Q )."""

    if( bdflsFile == None ) :
        import bdfls
        bdflsFile = bdfls.getDefaultBdfls( )
    import endlmisc
    ZA, Suffix = endlmisc.intZASuffix( ZA_ )
    if( ( ( ZA % 1000 ) / 100 ) > 7 ) : ZA = ZA % 100 + 1000 * ( ZA / 1000 )
    yi = endlmisc.incidentParticleTags( yi_ )[0]
    if( ( ( yi == 3 ) and ( ZA == 1003 ) and ( C == 30 ) ) or       # Added kludge to handle t(d,g n)He reaction.
        ( ( yi == 4 ) and ( ZA == 1002 ) and ( C == 30 ) ) ) :
        yos = [ 7, 1 ]
    else :
        yos = endl_C.endl_C_yoInfo( C )
    if ( yos == None ) or ( yos == () ) : return (None, ( ), None) # No such C-value or C-value does is not a reaction (e.g., C=1).
    if ( yos[0] == -1 ) : return ( ZA, ( yi, ), None )            # yo unknown.
    if ( yos[0] == -2 ) : return ( ZA, ( yi, ), 0. )              # yo = yi.
    if ( yos[0] == -3 ) : return ( None, ( -3, ), None )          # Fission
    if ( yos[0] == -4 ) : return ( None, ( -4, yos[1] ), None )  # Xyo
    if ( ZA >= 99000 ) and ( ZA <= 99200 ) : return( None, yos, None )
    Zr, Ar = ZandAFromZA( ZA )
    Z_A = ZandAFromZA( yoToZA_forNuclei( ZAToYo( yi ) ) )             # yi's Z and A.
    Zr += Z_A[0]
    if( Ar != 0 ) : Ar += Z_A[1]
    ZAr = ZA + yoToZA_forNuclei( ZAToYo( yi ) )                   # Add yi's ZA to ZAr.
    for yo in yos :                                         # Subtract yos' Z, A and ZA from ZAr.
        if( yo in [ 7, 8, 9 ] ) : continue
        yo_ = yo
        if ( yo < 1001 ) : yo_ = yoToZA_forNuclei( yo )
        Z_A = ZandAFromZA( yo_ )
        Zr -= Z_A[0]
        if( Ar != 0 ) :
            Ar -= Z_A[1]
            ZAr -= yo_
    if( Zr < 0 ) :                                          # Check all ok.
        raise Exception( "\nError in residualZA_yos_Q: bad C-value = %d for ZA = %d and yi = %d" % ( C, ZA, yi ) )
    if( specialCases and ( ZAr == 4008 ) ) :        # Special case for Be8 ->  2 He (e.g., n + Be9  -> 2 n + (Be8 ->  2 He)).
        yos = yos + ( 2004, )
        Zr, Ar, ZAr = 2, 4, 2004
    if( Ar < 0 ) : Ar = 0
    if( Ar == 0 ) : return( 1000 * Zr, yos, None )
    Remove_yi = True
    yos_ = []
    for yo in yos :                                         # To make math better, do not add and then subtract when yo = yi.
        try :
            yo_ = ZAToYo( yo )
        except :
            yo_ = yo
        if( Remove_yi and ( yo_ == yi ) ) :
            Remove_yi = False
        else :
            if ( yo != 0 ) : yos_.append( yo )
    if ( ZAr == ZA ) :
        ZAyos = 0
        for yo in yos :
            if( yo < 1000 ) :
                ZAyos += yoToZA_forNuclei( yo )
            else :
                ZAyos += yo
        q = bdflsFile.mass( yi )
        for yo in yos : q -= bdflsFile.mass( yo )
        if ( len( yos_ ) == 0 ) or ( ZAyos == yoToZA_forNuclei( yi ) ) :
            q = ( bdflsFile.constant( 4 ) * q )
            q += targetELevel - X4
            return ( ZAr, yos, q )
        endlmisc.printWarning( "Printing yo list")
        for yo in yos_ : endlmisc.printWarning( `yo` )
        endlmisc.printWarning( "\n" )
        raise Exception( "Error in residualZA_yos_Q: ZAr == ZA, but yo list not empty" )
    q = 0.
    if( Remove_yi ) : q += bdflsFile.mass( yi )
    for yo in yos_ : q -= bdflsFile.mass( yo )
    ZAm = bdflsFile.mass( ZA )
    if ( ZAm == None ) :
        if( printQWarning ) : endlmisc.printWarning( "Warning in residualZA_yos_Q: target ZA = %d not in bdfls file (yi = %d, ZA = %d, C = %d)" % \
            ( ZA, yi, ZA, C ) )
        q = None
    else :
        m = bdflsFile.mass( ZAr )                              # Mass of residual.
        if ( m == None ) : 
            if( printQWarning ) : endlmisc.printWarning( "Warning in residualZA_yos_Q: could not calculate Q as residual ZA = %d" % ZAr + 
                " not in bdfls file (yi = %d, ZA = %d, C = %d)" % ( yi, ZA, C ) )
            q = None
        else :
            q += bdflsFile.mass( ZA ) - m
            q *= bdflsFile.constant( 4 )                       # Convert AMU to MeVs.
    if q is not None : q += targetELevel - X4
    return ( ZAr, yos, q )

def reactionEquations( yi_, ZA_, C, level = None, specialCases = 1, printQWarning = True, NameASeperator = "", bdflsFile = None ) :
    """Given the incident particle yi, target ZA and reaction type C,
    this routine returns a list of three string for the reaction's equation.  The first string
    is of the form 'yi + ZA -> products + residual + Q' and the second is of the form
    ZA(yi,products)residual.  The third string is of the form (yi,products).  For example, 
    reactionEquations( 1, 94239, 13 ) returns the list 
    [ 'n + Pu239n -> 3n + Pu237 + -12.64503 MeV', 'Pu239(n,3n)Pu237', '(n,3n)' ].
    If level is not None, """

    import endlmisc
    yi = endlmisc.incidentParticleTags( yi_ )[0]
    ZA, Suffix = endlmisc.intZASuffix( ZA_ )
    Z, A = ZandAFromZA( ZA )
    s1 = endl_y.endl_yLabel( ZAToYo( yi ) ) + " + " + symbolForYoOrZA( ZA_, ZAOnly = 1, AddNatural = 0, NameASeperator = NameASeperator ) + " ->"
    s2 = symbolForYoOrZA( ZA_, ZAOnly = 1, AddNatural = 0, NameASeperator = NameASeperator ) + "(" + endl_y.endl_yLabel( ZAToYo( yi ) ) + ","
    ZAr, yos, q = residualZA_yos_Q( yi, ZA, C, specialCases = specialCases, printQWarning = printQWarning, bdflsFile = bdflsFile )
    if ( ( len( yos ) > 0 ) and ( yos[0] < 0 ) ) :
        if   ( yos[0] == -1 ) :     # yo is unknown
            s1 += " ?"
            s2 += "?)?"
        elif ( yos[0] == -2 ) :     # yo = yi
            s1 += " %s + %s" % ( endl_y.endl_yLabel( ZAToYo( yi ) ), symbolForYoOrZA( ZA_, ZAOnly = 1, AddNatural = 0, NameASeperator = NameASeperator ) )
            s2 += "%s)%s" % ( endl_y.endl_yLabel( ZAToYo( yi ) ), symbolForYoOrZA( ZA_, ZAOnly = 1, AddNatural = 0, NameASeperator = NameASeperator ) )
        elif ( yos[0] == -3 ) :     # fission
            s1 += " fission + ?"
            s2 += "fission)?"
        elif ( yos[0] == -4 ) :     # Xyo
            s1 += " X%s + ?" % endl_y.endl_yLabel( ZAToYo( yos[1] ) )
            s2 += "X%s)?" % endl_y.endl_yLabel( ZAToYo( yos[1] ) )
        else :
            raise Exception( "\nError in reactionName: invalid yos[0] = %s" % `yos[0]` )
    else :
        n = 1
        yo_ = None
        syo = ""
        for yo in yos :
            if ( yo_ == None ) :
                yo_ = yo
            elif ( yo != yo_ ) :
                syo += " + "
                if ( n > 1 ) :
                    syo += `n`
                    s2 += `n`
                syo += symbolForYoOrZA( yo_, NameASeperator = NameASeperator )
                s2 += symbolForYoOrZA( yo_, NameASeperator = NameASeperator ) + " "
                n = 1
                yo_ = yo
            else :
                n += 1
        if ( yo_ != None ) :
            syo += " + "
            if ( n > 1 ) :
                syo += `n`
                s2 += `n`
            syo += endl_y.endl_yLabel( ZAToYo( yo_ ) )
            s2 += endl_y.endl_yLabel( ZAToYo( yo_ ) ) + " "
            if ( ( yi == 1 ) and ( C == 11 ) ) or ( ( yi == 2 ) and ( C == 40 ) ) or ( ( yi == 3 ) and ( C == 41 ) ) or \
               ( ( yi == 4 ) and ( C == 42 ) ) or ( ( yi == 5 ) and ( C == 44 ) ) or ( ( yi == 6 ) and ( C == 45 ) ) : 
                syo += "'"
                s2 = s2[:-1] + "'"
        syo = syo[2:]                                       # remove first " +"
        if ( q == None ) :
            qStr = "?"
        else :
            qStr = "%.8g MeV" % q
        if( ZA > 99000 ) and ( ZA <= 99200 ) or ( ZAr == None ) :
            ZArStr = "?"
        else :
            ZArStr = symbolForYoOrZA( ZAr, ZAOnly = 1, AddNatural = 0, NameASeperator = NameASeperator )
        if( level != None ) : ZArStr += '[' + str( level ) + ']'
        s1 += syo + " + " + ZArStr + " + " + qStr
        if ( s2[-1] == " " ) : s2 = s2[:-1]
        s2 += ")" + ZArStr
    s3 = ""
    if ( s2.find( "(" ) != -1 ) : s3 = "(" + s2.split( "(" )[1].split( ")" )[0] + ")"
    return [ s1, s2, s3 ]

def reactionQ( yi, ZA, C, targetELevel = 0., X4 = 0, specialCases = 1, printQWarning = True, bdflsFile = None ) :
    """Returns the Q for reaction of type C involving incident particle yi and target ZA. Note,
    the target is taken to be in its ground state."""

    q = residualZA_yos_Q( yi, ZA, C, targetELevel = targetELevel, X4 = X4, specialCases = specialCases, printQWarning = printQWarning, bdflsFile = bdflsFile )
    if ( q[0] == None ) : return None
    return q[2]

def reactionThreshold( yi, ZA, C, targetELevel = 0., residualELevel = 0., Q = None, specialCases = 0, S = None, bdflsFile = None ) :
    """Returns the threshold for this reaction. targetELevel is the excited state of the target
    and residualELevel is the excited state of the residual."""

    if( bdflsFile == None ) :
        import bdfls
        bdflsFile = bdfls.getDefaultBdfls( )
    QRequired = False
    if( AFromZA( ZA ) == 0 ) :        # User must enter Q
        QRequired = True
    elif( 99120 <= ZA <= 99200 ) :          # User must enter Q
        QRequired = True
    elif( C == 1 ) :
        raise Exception( 'Error from reactionThreshold: cannot calculate threshold for C = 1 reaction' )
    elif( C == 5 ) :
        raise Exception( 'Error from reactionThreshold: cannot calculate threshold for C = 5 reaction' )
    elif( C == 15 ) :                       # User must enter Q
        QRequired = True
    elif( C >= 50 ) :                  # User must enter Q
        QRequired = True
    else :
        C_ = C
        if( S == 1 ) :
            yos = endl_C.endl_C_yoInfo( C )
            if( len( yos ) > 1 ) :              # Attempt to handle breakup of residual properly (e.g., n + C12 --> n + (C12' -> 3a )).
                yo = yos[:1]
                CList = endl_C.endl_C_CList( )
                for Cp in CList :
                    if( yo == endl_C.endl_C_yoInfo( Cp ) ) :
                        C_ = Cp
                        break
        Q = reactionQ( yi, ZA, C_, specialCases = specialCases, bdflsFile = bdflsFile )
    if( QRequired and ( Q == None ) ) : raise Exception( 'Error from reactionThreshold: user must supply Q for requested ZA = %s or C = %d' % ( `ZA`, C ) )
    Threshold = -Q  - targetELevel + residualELevel
    if( Threshold < 0. ) : Threshold = 0.
    r = 1. + bdflsFile.mass( yi ) / bdflsFile.mass( ZA )
    Threshold *= r
    return( Threshold )

def symbolForYoOrZA( yoOrZA, NameASeperator = "", ZAOnly = 0, AddNatural = 1, m_to_m1 = False, suffixSeperator = None ) :
    """If ZAOnly == 0 and yoOrZA is a yo type particle then the symbol for its yo is returned;
    otherwise, the symbol for ZA is returned.  If AddNatural is true, then natural ZAs (i.e.,
    those with A = 000 like 12000) have the string 'natural' append to the symbol for the Z, 
    otherwise, natural targets only contain the symbol. If m_to_m1 is True then when Suffix = 'm'
        it is changed to 'm1'."""

    return symOrNameForYoOrZA( yoOrZA, 0, ASep = NameASeperator, ZAOnly = ZAOnly, AddNatural = AddNatural, m_to_m1 = m_to_m1, suffixSeperator = suffixSeperator )

def nameForYoOrZA( yoOrZA, NameASeperator = "", ZAOnly = 0, AddNatural = 1, m_to_m1 = False, suffixSeperator = None ) :
    """If ZAOnly == 0 and yoOrZA is a yo type particle then the name for its yo is returned;
    otherwise, the name for ZA is returned. If AddNatural is true, then natural ZAs (i.e.,
    those with A = 000 like 12000) have the string 'natural' append to the name for the Z, 
    otherwise, natural targets only contain the name. If m_to_m1 is True then when Suffix = 'm'
    it is changed to 'm1'."""

    return symOrNameForYoOrZA( yoOrZA, 1, ASep = NameASeperator, ZAOnly = ZAOnly, AddNatural = AddNatural, m_to_m1 = m_to_m1, suffixSeperator = suffixSeperator )

def symOrNameForYoOrZA( yoOrZA, wantName, ASep = "", ZAOnly = 0, AddNatural = 1, m_to_m1 = False, suffixSeperator = None ) :
    """For internal use only (see symbolForYoOrZA or nameForYoOrZA)."""

    import endlmisc
    if( yoOrZA == 1 ) : return 'n'
    if( yoOrZA == 7 ) : return 'gamma'
    if ( ZAOnly == 0 ) :
        try :
            ZAyo = ZAToYo( yoOrZA )                         # Is ZA a valid yo.
        except :
            ZAOnly = 1
    if ( ZAOnly == 0 ) :
        if ( wantName == 0 ) : return endl_y.endl_yLabel( ZAyo )
        return endl_y.endl_yLongLabel( ZAyo )
    ZA, Suffix = endlmisc.intZASuffix( yoOrZA )
    if( m_to_m1 and ( Suffix == 'm' ) ) : Suffix = 'm1'
    if( ZA > 99000 ) and ( ZA <= 99200 ) : return "FissionProductENDL%d" % ZA
    if( suffixSeperator is None ) : suffixSeperator = ASep
    if( Suffix != "" ) : Suffix = suffixSeperator + Suffix
    Z, A = ZandAFromZA( ZA )
    if ( wantName == 0 ) :
        ZStr = endl_Z.endl_ZSymbol( Z )
        if ( ZStr == None ) : ZStr = "Unk"
    else : 
        ZStr = endl_Z.endl_ZLabel( Z )
        if ( ZStr == None ) : ZStr = "Unknown"
    if( A == 0 ) :
        if ( AddNatural ) : return ZStr + suffixSeperator + 'natural' + Suffix
        return ZStr
    return ZStr + ASep + `A` + Suffix

def endlToGNDName( yoOrZA ) :

    return( symbolForYoOrZA( yoOrZA, NameASeperator = '', ZAOnly = True, AddNatural = True, m_to_m1 = True, suffixSeperator = '_' ) )

def gndNameToEndlZ_A_Suffix( name ) :
    """Returns the tuple (Z, A, suffix, ZA) for a gnd isotope name (e.g., gnd name = 'Am242_m1' returns ( 95, 242, 'm1', 95242 )."""

    return( nuclear.getZ_A_suffix_andZAFromName( name ) )

def fileInfo( yi, ZA, fileName ) :
    """Obtains yo, C, I and S from fileName and calls reactionInfo."""

    import endlmisc
    ( yo, C, I, S ) = endlmisc.intYoCIS( fileName )
    return reactionInfo( yi, ZA, yo, C, I, S )

def reactionInfo( yi, ZA, yo, C, I, S, bdflsFile = None ) :
    """Returns a string of 5 lines. The first 3 lines are the item return by reactionEquations.
    The fourth is "Elastic" if C = 10 else it is a blank line, and the last is "total cross section"
    if I == 0 and C == 0, else it is the string return by endl_I.endl_ILabel."""

    r = reactionEquations( yi, ZA, C, bdflsFile = bdflsFile )
    s  = r[0] + "\n"                            # Line 1 done
    s += r[1] + "\n"                            # Line 2 done
    s += r[2] + "\n"                            # Line 3 done
    if   ( C == 10 ) : s += "Elastic"
    s += "\n"                                   # Line 4 done
    if ( I == 0 ) and ( C == 0 ) :
        s += "total cross section"
    else :
        s += endl_I.endl_ILabel( I )
    s += "\n"                                   # Line 5 done
    return s

def reactionQByZAs( incomingZAs, outgoingZAs, bdflsFile = None ) :
    """Returns the Q in MeV for the reaction, that is, the (sum of incoming masses) - (sum of outgoing mass) in MeV.
    Returns None if at least one of the masses is unknown."""

    if( bdflsFile == None ) :
        import bdfls
        bdflsFile = bdfls.getDefaultBdfls( )
    ZAMultiplicities = {}
    for ZA in incomingZAs :
        if( ZA < 1000 ) : ZA = yoToZA_forNuclei( ZA )
        if( ZA not in ZAMultiplicities ) : ZAMultiplicities[ZA] = 0
        ZAMultiplicities[ZA] += 1
    for ZA in outgoingZAs :
        if( ZA < 1000 ) : ZA = yoToZA_forNuclei( ZA )
        if( ZA not in ZAMultiplicities ) : ZAMultiplicities[ZA] = 0
        ZAMultiplicities[ZA] -= 1
    mass = 0.
    ZAs = ZAMultiplicities.keys( )
    ZAs.sort( )
    for ZA in ZAs :
        m = bdflsFile.mass( ZA )
        if( m == None ) : return( None )
        mass += m * ZAMultiplicities[ZA]
    return( mass * bdflsFile.constant( 4 ) )

def reactionQByIsotopeNames( incomingNames, outgoingNames, bdflsFile = None ) :
    """Returns the Q in MeV for the reaction, that is, the (sum of incoming masses) - (sum of outgoing mass) in MeV.
    Returns None if at least one of the masses is unknown. This routine only works for isotopes in the ground state as it
    converts each name to it ZA and calls reactionQByZAs (see reactionQByZAs for more details). Example,
    incomingNames = [ 'n', 'W186' ]
    outgoingNames = [ 'H3', 'He4', 'Tm176' ]
    Q = reactionQByIsotopeNames( incomingNames, outgoingNames )"""

    incomingZAs = [ gndNameToEndlZ_A_Suffix( name )[-1] for name in incomingNames ]
    outgoingZAs = [ gndNameToEndlZ_A_Suffix( name )[-1] for name in outgoingNames ]
    return( reactionQByZAs( incomingZAs, outgoingZAs, bdflsFile ) )

def possibleReactions( projectile, target, energy_MeV, maxProducts = 4, bdflsFile = None ) :
    """This function returns a list of all reactions for the projectile with kinetic energy energy_MeV 
    hitting the stationary target that are allowed by the amount of energy_MeV available in the 
    center-of-mass and a reaction's Q.  Reactions considered will only contain products of n_1, 
    H_1, H_2, H_3, He_3, He_4 and a residual. The reactions considered are also limited to the 
    number of outgoing products, not including the residual, being no greater than maxProducts."""

    if( bdflsFile == None ) :
        import bdfls
        bdflsFile = bdfls.getDefaultBdfls( )
    import endlmisc

    projectile_ = endlmisc.incidentParticleTags( projectile )
    projectiles = [ [ 1, 0 ], [ 1001, 0 ], [ 1002, 0 ], [ 1003, 0 ], [ 2003, 0 ], [ 2004, 0 ] ]
    nProjectiles = len( projectiles )
    if( not ( 0 < projectile_[0] <= nProjectiles ) and projectile != 7 ) : raise Exception( '%s is an unsupported projectile' % projectile )
    projectileZA = yoToZA( projectile_[0] )
    incidentZA = projectileZA + target

    if projectile != 7:
        for p in projectiles :
            p[1] = bdflsFile.mass( p[0] )
            if( p[0] == projectileZA ) : projectileMass = p[1]
    else: projectileMass = 0.0
    targetMass = bdflsFile.mass( target )
    incidentMass = projectileMass + targetMass
    kineticMass = energy_MeV / bdflsFile.constant( 4 )
    if( kineticMass < 1e-6 * projectileMass ) :                         # Use non-relativity
        incidentMass += targetMass * kineticMass / incidentMass
    else :                                                              # Need relativity
        incidentMass = math.sqrt( incidentMass * incidentMass + 2. * kineticMass * targetMass )

    reactions = []

    def possibleReactions2( currentZA, products, productMass, level ) :

        if( level > maxProducts ) : return
        for ZA, mass in projectiles :
            residualZA = currentZA - ZA
            if( residualZA < ZA ) : break
            Z, A = ZandAFromZA( residualZA )
            if( A < Z ) : continue
            newProducts = products + [ ZA ]
            newProductMass = productMass + mass
            if( not( 1 < residualZA < 1001 ) ) : 
                residualMass = bdflsFile.mass( residualZA )
                if( residualMass != None ) : 
                    if( incidentMass >= ( residualMass + newProductMass ) ) : reactions.append(  newProducts + [ residualZA ] )
            possibleReactions2( residualZA, newProducts, newProductMass, level + 1 )

    possibleReactions2( incidentZA, [], 0., 1 )
    return( reactions )
