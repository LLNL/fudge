/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <limits.h>
#include <math.h>

#include "coulfg2.h"

#define ACCUR 1.0e-16

static int JWKB( double XX, double ETA, double XL, double *FJWKB, double *GJWKB, int IEXP );
static int IDINT( double d );
static double DSIGN( double a, double b );
static int MAX0( int i1, int i2 );
static double DMIN1( double a, double b );
static double DMAX1( double a, double b );
/*
************************************************************
*/
int Coulomb_FG( double XX, double ETA, double XLMIN, double XLMAX, double *FC, double *GC, int *M1 ) {
/*
C  REVISED IJT WITH L-T ALGORITHMN FOR CONTINUED FRACTIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL      "      "     "                C
C            = 2 CYLINDRICAL BESSEL      "      "     "                C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR DOUBLEPRECISION ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/
    int ETANE0, XLTURN;
    int IEXP = 1, L, L1, LXTRA, LP;
    double ACC, ACCH, FCL, E2MM1, F;
/*  double ACC4; */
    double GAM, W, BETA, FCM, GCL, GPL, ALPHA, A, B, XL, EL, RL, SL, GCL1, GJWKB, FJWKB = 0, FPL, FCL1;
    double DELL, XLL, XI, PK, PX, D, C, PK1, EK, RK2, TK, DF, TA, WI, P, Q, AR, AI, BR, BI, DR, DI, DP, DQ;
    double PACCQ = 1.0, ABORT = 2.0e4;

    GJWKB = 0.0;
    ETANE0  = ETA != 0.0;
    ACC   = ACCUR * 10.0;
/*  ACC4  = ACC * 1e4; */
    ACCH  = sqrt( ACC );

    RL = 0; EL = 0; XL = 0;

    if( XX <= ACCH ) return( -1 );      /* TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE */
    if( ( XLMIN <= -1.0 ) || ( XLMAX < XLMIN) ) return( -2 );
    E2MM1 = ETA * ETA + XLMIN * XLMIN + XLMIN;
    XLTURN = ( XX * ( XX - 2.0 * ETA ) ) < ( XLMIN * XLMIN + XLMIN );
    DELL  = XLMAX - XLMIN + ACC;
    F = fmod( DELL, 1.0 );
/*    if( fabs( F ) > ACC4 ) WRITE(6,2040) XLMAX, XLMIN, DELL, F ????????? */
    LXTRA = IDINT( DELL );
    XLL   = XLMIN + (double) LXTRA;
/* LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED XLL IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN */
    *M1 = MAX0( IDINT( XLMIN + ACC ), 0 );
    L1 = *M1 + LXTRA;
/*
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,XX)/F(XL,ETA,XX)
*/
    XI  = 1.0 / XX;
    FCL = 1.0;
    PK  = XLL + 1.0;
    PX  = PK  + ABORT;
    F   =  ETA / PK + PK * XI;
    if( fabs( F ) < 1e-30 ) F = 1e-30;
    D = 0.0;
    C = F;
/*
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
*/
    do {
        PK1 = PK + 1.0;
        EK  = ETA / PK;
        RK2 = 1.0 + EK * EK;
        TK  = ( PK + PK1 ) * ( XI + EK / PK1 );
        D   =  TK - RK2 * D;
        C   =  TK - RK2 / C;
        if( fabs( C ) < 1e-30 ) C = 1e-30;
        if( fabs( D ) < 1e-30 ) D = 1e-30;
        D = 1.0 / D;
        DF = D * C;
        F  = F * DF;
        if( D < 0.0 ) FCL = -FCL;
        PK = PK1;
        if( PK > PX ) return( -10 );            /* Failed to converge. */
    } while( fabs( DF - 1.0 ) >= ACC );
    if( LXTRA != 0 ) {
/*
C *** DOWNWARD RECURRENCE TO LAMBDA = XLMIN. ARRAY GC, IF PRESENT, STORES RL
*/
        FCL = FCL * 1e-30;
        FPL = FCL * F;
        FC[L1] = FCL;
        XL  = XLL;
        RL  = 1.0;
        EL  = 0.0;
        for( LP = 0; LP < LXTRA; LP++ ) {
            if( ETANE0 ) {
                EL = ETA / XL;
                RL = sqrt( 1.0 + EL * EL );
            }
            SL    =  EL  + XL * XI;
            L     =  L1  - LP;
            FCL1  = ( FCL * SL + FPL ) / RL;
            FPL   =  FCL1 * SL - FCL * RL;
            FCL   =  FCL1;
            FC[L-1] =  FCL;
            if( ETANE0 ) GC[L] = RL;
            XL = XL - 1.0;
        }
        if( FCL == 0.0 ) FCL = ACC;
        F  = FPL / FCL;
    }
/****    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM */
    if( XLTURN ) JWKB( XX, ETA, DMAX1( XLMIN, 0.0 ), &FJWKB, &GJWKB, IEXP );
    if( !( ( IEXP > 1 ) || ( GJWKB > ( 1.0 / ( ACCH * 100.0 ) ) ) ) ) {
        XLTURN = 0;
        TA = 2.0 * ABORT;
        PK = 0.0;
        WI = ETA + ETA;
        P  = 0.0;
        Q  = 1.0 - ETA * XI;
        AR = -E2MM1;
        AI =  ETA;
        BR =  2.0 * ( XX - ETA );
        BI =  2.0;
        DR =  BR / ( BR * BR + BI * BI );
        DI = -BI / ( BR * BR + BI * BI );
        DP = -XI * ( AR * DI + AI * DR );
        DQ =  XI * ( AR * DR - AI * DI );
        do {
            P     = P  + DP;
            Q  = Q  + DQ;
            PK = PK + 2.0;
            AR = AR + PK;
            AI = AI + WI;
            BI = BI + 2.0;
            D  = AR * DR - AI * DI + BR;
            DI = AI * DR + AR * DI + BI;
            C  = 1.0 / ( D * D + DI * DI );
            DR =  C * D;
            DI = -C * DI;
            A  = BR * DR - BI * DI - 1.0;
            B  = BI * DR + BR * DI;
            C  = DP * A  - DQ * B;
            DQ = DP * B  + DQ * A;
            DP = C;
            if( PK > TA ) return( -20 );
        } while( ( fabs( DP ) + fabs( DQ ) ) >= ( ( fabs( P ) + fabs( Q ) ) * ACC ) );
        PACCQ = 0.5 * ACC / DMIN1( fabs( Q ), 1.0 );
        if ( fabs( P ) > fabs( Q ) ) PACCQ = PACCQ * fabs( P );
/*
C *** SOLVE FOR FCM = F AT LAMBDA = XLMIN, THEN FIND NORM FACTOR W=W/FCM
*/
        GAM = (F - P) / Q;
/*
        if( Q <= ( ACC4 * fabs( P ) ) ) write(4,*) 'COULFG Error:', ' Q,P,XX,XLMIN =', real(q), real(p), real(XX), xlmin
*/
        W   = 1.0 / sqrt( ( F - P ) * GAM + Q ); }
    else {
/* *** ARRIVE HERE IF G(XLMIN) .GT. 10**6 OR IEXP .GT. 70 & XLTURN = .TRUE. */
        W   = FJWKB;
        GAM = GJWKB * W;
        P   = F;
        Q   = 1.0;
    }
/*
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
*/
    ALPHA = 0.0;
    BETA  = 1.0;
    FCM  = DSIGN( W, FCL ) * BETA;
    FC[*M1]  = FCM;
    if( XLTURN ) {
        GCL = GJWKB * BETA; }
    else {
        GCL = FCM * GAM;
    }
    GC[*M1] = GCL;
    GPL = GCL * ( P - Q / GAM ) - ALPHA * GCL;
    if( LXTRA == 0 ) return( 0 );
/**** UPWARD RECURRENCE FROM GC(M1), GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC, FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLMIN HERE  AND RL = 1.0, EL = 0.0 FOR BESSELS */
    W    = BETA * W / fabs( FCL );
    for( L = *M1; L < L1; L++ ) {
        XL = XL + 1.0;
        if( ETANE0 ) EL = ETA / XL;
        if( ETANE0 ) RL = GC[L+1];
        SL = EL + XL * XI;
        GCL1    = ( ( SL - ALPHA ) * GCL - GPL ) / RL;
        GPL     = RL * GCL - ( SL + ALPHA ) * GCL1;
        GCL     = GCL1;
        GC[L+1] = GCL1;
        FC[L+1] = W * FC[L+1];
    }
    return( 0 );
}
/*
************************************************************
*/
static int JWKB( double XX, double ETA, double XL, double *FJWKB, double *GJWKB, int IEXP ) {
/*
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C *** CALLS DMAX1,SQRT,DLOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981
*/
    double GH2, XLL1, HLL, HL, SL, RL2, GH, PHI, PHI10;

    GH2  = XX * ( ETA + ETA - XX );
    XLL1 = DMAX1( XL * XL + XL, 0.0 );
    if( ( GH2 + XLL1 ) <= 0.0 ) return( IEXP );
    HLL  = XLL1 + 6.0 / 35.0;
    HL   = sqrt( HLL );
    SL   = ETA / HL + HL / XX;
    RL2  = 1.0 + ETA * ETA / HLL;
    GH   = sqrt( GH2 + HLL ) / XX;
    PHI  = XX * GH - 0.5 * ( HL * log( ( GH + SL ) * ( GH + SL ) / RL2 ) - log( GH ) );
    if( ETA != 0.0 ) PHI = PHI - ETA * atan2( XX * GH, XX - ETA );
    PHI10 = -PHI * 0.43429448190325182;           /* 0.43429448190325182 = log10( e ). */
    IEXP = (int) PHI10;
    if( IEXP > 70 ) {
        *GJWKB = pow( 10.0, PHI10 - (double) IEXP ); }
    else {
      *GJWKB = exp( -PHI );
      IEXP  = 0;
    }
    *FJWKB = 0.5 / ( GH * *GJWKB );
    return( IEXP );
}
/*
************************************************************
*/
static int IDINT( double d ) {

    if( d < INT_MIN ) return( INT_MIN );
    if( d > INT_MAX ) return( INT_MAX );
    return( (int) d );
}
/*
************************************************************
*/
static double DSIGN( double a, double b ) {

    if( b < 0.0 ) return( -fabs( a ) );
    return( fabs( a ) );
}
/*
************************************************************
*/
static int MAX0( int i1, int i2 ) {

    if( i1 > i2 ) return( i1 );
    return( i2 );
}
/*
************************************************************
*/
static double DMIN1( double a, double b ) {

    if( a > b ) return( a );
    return( b );
}
/*
************************************************************
*/
static double DMAX1( double a, double b ) {

    if( a > b ) return( a );
    return( b );
}
