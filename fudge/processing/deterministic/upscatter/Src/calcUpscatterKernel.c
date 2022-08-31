/*
  # <<BEGIN-copyright>>
  Copyright 2022, Lawrence Livermore National Security, LLC.
  See the top-level COPYRIGHT file for details.
  
  SPDX-License-Identifier: BSD-3-Clause
  # <<END-copyright>>
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
   This code generates a table of E, mu, E' and P(T,E,mu,E') for elastic scattering of a neutron off a target of mass M at
   temperature T, where E and E' are the incident and outgoing neutron's energies, mu is the cosine of the angle between 
   the incident and outgoing neutron's velocities and P is the unnormalized probability for a neutron to be produced with 
   mu, E' given E and T. The formula for P(T,E,mu,E') is from "The Slowing Down and Thermalization of Neutrons", M.M.R. 
   Williams, 1966. In the reference, P(T,E,mu,E') is defined as

       P(T,E,mu,E') = Exp( -( a^2 + b^2 ) / ( 4 a ) ) Exp( -b/2 ) sqrt( E'/ E ) / ( T sqrt( 4 pi a ) )

   where 
       b = -( E - E' ) / T, 
       a = k_vec^2 / ( 2 M T ), 
       k_vec = K_vec - K_vec', and K_vec (K_vec') in the momentum of the incident (outgoing) neutron.
       a = k_vec^2 / ( 2 M T ) = ( K_vec - K_vec' ) * ( K_vec - K_vec' ) / ( 2 M T ) = ( E + E' - 2 sqrt( E E' ) mu ) * ( m / M ) / T.
   For numerical reasons, it is more convenient to write 

       Exp( -( a^2 + b^2 ) / ( 4 a ) ) Exp( -b/2 ) = Exp( -( a + b )^2 / ( 4 a ) )

   xsec

*/

/*

TODO:
    Cleanups:
        - remove older file read write after everything is stable. (i.e.  int getLine( FILE *f, char *fileName )  )
        -  
    
    Additions:
        - pass groups as an input instead of a file
        - 
 
*/

#define Pi 3.1415926535897932385
static double xsecMin = 0., T_MeV;

#define bufferSize 512
static int lMaxPlus1 = 4;

#define myPrintf myPrintfNothing

typedef
    struct EpX {
        double Ep, X;
    } EpX;

typedef
    struct EEpMuP {
        long Index;
        double E, Ep, Mu, P;
    } EEpMuP;

typedef
    struct DummyList {
        int n, nAllocated;
        double Dummy;
        void *v;
    } DummyList;

typedef
    struct EpList {
        int n, nAllocated;
        double Mu;
        EpX *EpX;
    } EpList;

typedef
    struct MuList {
        int n, nAllocated;
        double E;
        EpList *EpL;
    } MuList;

typedef
    struct XEplList {
        double E;
        double *XEpl;
    } XEplList;
    
typedef
    struct XEEplList {
        int n, nAllocated;
        double Dummy;
        XEplList *XEplL;
    }  XEEplList;
    
void printTransferMatrixAndData( int ng, double *gb, int lMax, double *XEEpl , int gMax );
void interpolateXE( int n, int lMax, double E, XEplList *priorXEplL, XEplList *nextXEplL, XEplList *XEplL );
void convertEMuEp_To_EEpl( XEEplList *XEEplL, MuList *MuL, int lMax, double E, int ng, double *gb, FILE *fp, int gMax );
double CalculateEpAve( MuList *MuL, double E );
static void P_ls( int l, double mu, double Pls[5] );
static double muPl( int l, double mu1, double X1, double mu2, double X2 );
void addToMuList( MuList *l, EpList *EpL, double Mu );
void addToEpList( EpList *l, double d1, double d2, double d3, double d4 );
int mergeXXYLists( int nx, double *x,  int nxy, double *xy, int nmerge, double *mergexxy );
void checkListSizeToAdd( DummyList *l, size_t s, int nToAdd, int increment );
void *myMalloc( size_t s, char *m );
void *myRealloc( void *p, size_t s, char *m );
void myFree( void *p );
int getLine( FILE *f, char *fileName );
void myPrintfNothing( char *s, ... );
void printMsg( char *s );

static int fileLine = 0;
static double E, Mu, Ep, prob ;

void calcKernel( double A, double T_MeV, double MinEValue,  double MaxEValue, double* gb, int ng, int* gMax ); 
int checkPoint( int indexMin, int indexMax, double y, double Mu, double A );
int insertPoint( int indexPrior, double x, double P );
void insertEEMPpoint( double E, double Ep, double Mu, double P );
int getNextEEMPline(  );
void clearPoints( void );
double xsec( double y, double x, double Mu, double A );

typedef
    struct xPData_s {
        int Prior;
        int Next;
        int Index;
        double x, P;
    } xPData;

static double EmaxScalingFactor = 200. ;
static double fraction = 0.02;
static int nxPData = 0, currentxPData = 0;
static xPData *xPDatas = NULL;

static long nEEMP = 0, currentEEMP;
static EEpMuP *EEMPdata = NULL;

/*
============================================================
*/
int main( int argc, char **argv ) {

    int i, j, k, l, gid, fid, ng, ngb, nf, iMax, nTmp, gMax;
    double EPrior = -1., muPrior = -1., *gb, *fe, *ff, *d1, *d2, *XEEpl;
    double MinEValue = 1e-11, MaxEValue = 20., massRatio;
    int lMax;
    char *errStrg;
    EpList EpL = { 0, 0, 0., NULL };
    MuList MuL = { 0, 0, 0., NULL };
    XEEplList XEEplL = { 0, 0, 0., NULL };
    XEplList tmpXEplL1 = { 0., NULL }, tmpXEplL2 = { 0., NULL }, *priorXEplL, *nextXEplL;
    FILE *f,*g,*fp;

    // get inputs
    massRatio = strtod( argv[2], &errStrg );
    if( *errStrg != 0 ) printMsg( "could not convert mass ratio to double" );
    T_MeV = strtod( argv[3], &errStrg );
    if( *errStrg != 0 ) printMsg( "could not convert T to double" );
    MinEValue = strtod( argv[4], &errStrg);
    if( *errStrg != 0 ) printMsg( "could not convert minE to double" );
    MaxEValue = strtod( argv[5], &errStrg );
    if( *errStrg != 0 ) printMsg( "could not convert maxE to double" );
    lMax = strtol( argv[6], &errStrg, 10 );
    if( *errStrg != 0 ) printMsg( "could not convert lMax to int" ); 
    lMaxPlus1 = lMax + 1;
    
    // read group boundaries from an input file
    g = fopen( argv[1], "r" );
    if( g == NULL ) {
        fprintf( stderr, "Could not open input file '%s'", argv[1] );
        exit( EXIT_SUCCESS );
    }
    fscanf(g, "%i", &ngb);
    gb = (double *) myMalloc( ngb * sizeof( double ), "gb" );                                            // AF1, Freed at end of main. 
    float fltP;
    for( i = 0; i < ngb; i++ ){ fscanf(g, "%e", &fltP); gb[i] = fltP;  }
    fclose(g);
    ng = ngb - 1;
    
    calcKernel( massRatio, T_MeV, MinEValue, MaxEValue, gb, ng, &gMax );
   
    fp = fopen ("AveEnergy.out", "w+");
    currentEEMP = 0;
    while( ( l = getNextEEMPline( ) ) == 0 ) {
        if( ( EPrior >= 0. ) && ( EPrior != E ) ) {
            addToMuList( &MuL, &EpL, muPrior );
            convertEMuEp_To_EEpl( &XEEplL, &MuL, lMaxPlus1, EPrior, ng, gb, fp, gMax );
            muPrior = Mu;
            EpL.n = 0;
            addToEpList( &EpL, E, Mu, Ep, prob );
            EPrior = E; }
        else if( Mu != muPrior ) {
            addToMuList( &MuL, &EpL, muPrior );
            EpL.n = 0;
            muPrior = Mu; }
        else {
            addToEpList( &EpL, E, Mu, Ep, prob );
        }
        EPrior = E;
    }
    convertEMuEp_To_EEpl( &XEEplL, &MuL, lMaxPlus1, EPrior, ng, gb, fp, gMax  );
    fclose(fp);


    // group the legendre matrix
    XEEpl = (double *) myMalloc( ( ( lMaxPlus1 + 1 ) * ng * ng ) * sizeof( double ), "XEEpl" );               // AF9, Freed at end of main.
    for( i = 0; i < ( lMaxPlus1 + 1 ) * ng * ng; i++ ) XEEpl[i] = 0.;                                         // XEEpl is Transfer matrix stored as (l,g,h).
    tmpXEplL1.E = XEEplL.XEplL[0].E;
    nTmp = ( lMaxPlus1 + 1 ) * ng;
    tmpXEplL1.XEpl = myMalloc( nTmp * sizeof( double ), "tmpXEplL1.XEpl" );                                 // AF10, Freed at end of main.
    tmpXEplL2.XEpl = myMalloc( nTmp * sizeof( double ), "tmpXEplL2.XEpl" );                                 // AF11, Freed at end of main.
    for( i = 0; i < ( ( lMaxPlus1 + 1 ) * ng ); i++ ) tmpXEplL1.XEpl[i] = 0.;
    priorXEplL = &tmpXEplL1;
    for( j = 0; j < XEEplL.n; j++ ) {
        if( XEEplL.XEplL[j].E > gb[0] ) break;
        priorXEplL = &(XEEplL.XEplL[j]);
    }
    for( i = 0; i < ng; i++ ) if( gb[i + 1] < XEEplL.XEplL[XEEplL.n - 1].E ) iMax = i + 1;
    for( i = 0; i < ng; i++ ) if( priorXEplL->E < gb[i + 1] ) break;
    if( priorXEplL != &tmpXEplL1 ) interpolateXE( ng, lMaxPlus1, gb[i], priorXEplL, &(XEEplL.XEplL[j]), &tmpXEplL1 );
    for(      ; i < iMax; i++ ) {
        while( ( priorXEplL->E < gb[i + 1] ) && ( j < XEEplL.n ) ) {
            if( XEEplL.XEplL[j].E <= gb[i + 1] ) {
                nextXEplL =  &(XEEplL.XEplL[j]); }
            else {
                interpolateXE( ng, lMaxPlus1, gb[i + 1], priorXEplL, &(XEEplL.XEplL[j]), &tmpXEplL1 );
                nextXEplL = &tmpXEplL1;
            }
            for( l = 0; l <= lMaxPlus1; l++ ) {
                for( k = 0; k < ng; k++ ) {
                    XEEpl[( l * ng + i ) * ng + k] += 0.5 * ( nextXEplL->XEpl[k + ng * l] + priorXEplL->XEpl[k + ng * l] ) * ( nextXEplL->E - priorXEplL->E );
                }
            }
            if( nextXEplL == &tmpXEplL1 ) {
                tmpXEplL2.E = tmpXEplL1.E;
                for( k = 0; k < nTmp; k++ ) tmpXEplL2.XEpl[k] = tmpXEplL1.XEpl[k];
                priorXEplL = &tmpXEplL2; }
            else {
                priorXEplL = nextXEplL;
                j++;
            }
        }
    }
    printTransferMatrixAndData( ng, gb, lMaxPlus1, XEEpl, gMax);

    myFree( XEEpl );                                                                                     // AF9, Allocated in main. 
    myFree( tmpXEplL1.XEpl );                                                                            // AF10, Allocated in main.
    myFree( tmpXEplL2.XEpl );                                                                            // AF11, Allocated in main.
    myFree( gb );                                                                                        // AF1, Allocated in main.  
    for( i = 0; i < XEEplL.n; i++ ) myFree( XEEplL.XEplL[i].XEpl );                                      // AF2, Allocated in convertEMuEp_To_EEpl. 
    myFree( XEEplL.XEplL );                                                                              // AF6, Allocated in convertEMuEp_To_EEpl.
    myFree( MuL.EpL );                                                                                   // AF7, Allocated in addToMuList. 
    myFree( EpL.EpX );                                                                                   // AF8, Allocated in addToEpList. 

    exit( EXIT_SUCCESS );
}
/*
============================================================
*/
void printTransferMatrixAndData( int ng, double *gb, int lMax, double *XEEpl, int gMax ) {

    int i, j, k, g, h, l, inc = 6;
    double *p = XEEpl, *d, s, *fluxAve;

    fluxAve = myMalloc( ng * sizeof( double ), "printTransferMatrixAndData:fluxAve" );                    // AF13, free at end of routine. 
    for( i = 0; i < ng; i++ ) fluxAve[i] = gb[i + 1] - gb[i];

    FILE *fp;
    fp = fopen ("upscatterLegendre.out", "w+");
   
    fprintf( fp, "upscatter: version 1 \n");
    fprintf( fp, "maxIncidentEnergyGroup: %d \n",gMax);
    fprintf( fp, "outputLegendreOrder: %d \n",lMax-1);
    fprintf( fp, "numEinBins: %d  \n",ng);
    fprintf( fp, "numEoutBins: %d  \n",ng);
    fprintf( fp, "Integrals, weight = 1: numEinBins = %d  \n",ng);
    
    for( g = ng-1; g >= 0; g-- ) {
        fprintf( fp, " EinBin = %d : numEoutBins = %d \n", ng-g-1, ng  );
        for( h = ng-1; h >= 0; h-- ) {
            for( l = 0; l < lMax; l++ ) {
                p = &(XEEpl[( l + 1 ) * ng * ng - g * ng  - h - 1 ]);
                fprintf( fp, "%e ", ( *p/fluxAve[ng - g - 1] )); 
            }
            fprintf( fp, "\n" );
        }
    }

    fclose(fp);
    
    myFree( fluxAve );                                                                                    // AF13, Allocated at beginning of routine. 
}
/*
============================================================
*/
void interpolateXE( int n, int lMax, double E, XEplList *priorXEplL, XEplList *nextXEplL, XEplList *XEplL ) {

    int i;
    double E1 = priorXEplL->E, E2 = nextXEplL->E, dEinv = ( E2 - E1 );

    if( dEinv != 0. ) dEinv = 1. / dEinv;
    XEplL->E = E;
    for( i = 0; i < n * ( lMax + 1 ); i++ ) XEplL->XEpl[i] = dEinv * ( priorXEplL->XEpl[i] * ( E2 - E ) + nextXEplL->XEpl[i] * ( E - E1 ) );
}
/*
============================================================
*/
void convertEMuEp_To_EEpl( XEEplList *XEEplL, MuList *MuL, int lMax, double E, int ng, double *gb, FILE *fp, int gMax ) {

    int i, j, l, n = 0, npEpX;
    double EpAve, *XMuEp, *pEpX, *p, *XEpl;
    DummyList *d = (DummyList *) XEEplL;
    EpList *EpL;

    EpAve = CalculateEpAve( MuL, E );
    fprintf( fp, "E = %12.5e %12.5e\n", E, EpAve );
    
    checkListSizeToAdd( d, sizeof( XEplList ), 1, 10000 );                                                // AF6, Freed at end of main. 
    XEEplL->XEplL[XEEplL->n].E = E;
    XEEplL->XEplL[XEEplL->n].XEpl = myMalloc( ng * ( lMax + 1 ) * sizeof( double ), "XEEplL->XEplL[XEEplL->n].XEpl" ); // AF2, Freed at end of main. 
    XMuEp = (double *) myMalloc( ng * MuL->n * sizeof( double ), "XMuEp" );                                // AF3, Freed at end of routine.
    for( i = 0; i < MuL->n; i++ ) if( n < MuL->EpL[i].n ) n = MuL->EpL[i].n;
    n += ng + 1;
    if( n < MuL->n ) n = MuL->n;
    n += 2;
    pEpX = (double *) myMalloc( 2 * n * sizeof( double ), "tmp" );                                        // AF4, Freed at end of routine.
    for( i = 0; i < MuL->n; i++ ) {                                                                        // Group the data for E'. 
        EpL = &(MuL->EpL[i]);
        npEpX = mergeXXYLists( ng + 1, gb,  EpL->n, (double *) EpL->EpX, n, pEpX );
        p = pEpX;
        for( j = 0; j < ng; j++ ) {
            XMuEp[i * ng + j] = 0.;
            while( ( *p < gb[j + 1] ) || ( ( npEpX > 1 ) && ( j == ( ng - 1 ) ) ) ) {
                XMuEp[i * ng + j] += 0.5 * ( p[1] + p[3] ) * ( p[2] - p[0] );
                p += 2;
                npEpX--;
            }
        }
        
    }

    for( i = 0; i < MuL->n; i++ ) pEpX[i] = MuL->EpL[i].Mu;
    XEpl = XEEplL->XEplL[XEEplL->n].XEpl;
    for( l = 0; l <= lMax; l++ ) {                                                                        // Calculate Legendre coefficients. 
        for( j = 0; j < ng; j++ ) {
            XEpl[l * ng + j] = 0.;
            for( i = 1; i < MuL->n; i++ ) XEpl[l * ng + j] += muPl( l, pEpX[i - 1], XMuEp[j + ng * ( i -  1) ], pEpX[i], XMuEp[j + ng * i] );
        }
    }
    XEEplL->n++;
    myFree( XMuEp );                                                                                    // AF3, Allocated at beginning of routine. 
    myFree( pEpX );                                                                                        // AF4, Allocated at beginning of routine. 

    for( i = 0; i < MuL->n; i++ ) myFree( MuL->EpL[i].EpX );                                            // AF5, Allocated in addToEpList. 
    MuL->n = 0;
}
/*
============================================================
*/
static void P_ls( int l, double mu, double Pls[5] ) {
//  Returns P_l(mu) for l = (l-2, l-1, l, l+1 and l+2) where P_l is the Legendre polynomial of degree l.

    double Pl = 0.,  Plp1 = 1.,  l_, twoLp1 = 1, Plm1 = 0., Plm2 = 0.;

    for( l_ = 0; l_ < l + 1; l_++, twoLp1 += 2 ) {
        Plm2 = Plm1;
        Plm1 = Pl;
        Pl = Plp1;
        Plp1 = ( twoLp1 * mu * Pl - l_ * Plm1 ) / ( l_ + 1 );
    }
    Pls[0] = Plm2;
    Pls[1] = Plm1;
    Pls[2] = Pl;
    Pls[3] = Plp1;
    Pls[4] = ( twoLp1 * mu * Plp1 - l_ * Pl ) / ( l_ + 1 );
}
/*
============================================================
*/
static double muPl( int l, double mu1, double X1, double mu2, double X2 ) {
/*
   This function integrates_mu1^mu2 dMu P_l(mu) * L(mu), where L(mu) is the linear function ( X2 - X1 ) ( mu - mu1 ) / ( mu2 - mu1 ) + X1.
   Or, integrates_mu1^mu2 dMu P_l(mu) * ( mu * b + a ), where a = ( X2 - X1 ) / ( mu2 - mu1 ) and b = X1 - mu1 * ( X2 - X1 ) / ( mu2 - mu1 ).


   This new muPl is good for all orders in l.
   The integral of P_l(mu) is ( P_{l+1}(mu) - P_{l-1}(mu) ) / ( 2 * l + 1 ). The integral of mu P_l(mu) can be done with "integration by parts" as
   integral( v * u' ) = v * u - integral( v' * u ) with u' = P_l(mu) and v = mu. Integration by parts yields
       integral( mu * P_l ) = ( x * ( P_{l+1} - P_{l-1} ) + 2 * ( 2 * l + 1 ) P_l  / ( 2 * l + 3 ) / ( 2 * l - 1 ) 
                               - P_{l+1} / ( 2 * l + 3 ) - P_{l+1} / ( 2 * l - 1 ) ) / ( 2 * l + 1 ).
*/
    long tlm1 = 2 * l - 1, tlp1 = tlm1 + 2, tlp3 = tlp1 + 2;
    double s, a = ( X2 - X1 ) / ( mu2 - mu1 ), b = X1 - mu1 * ( X2 - X1 ) / ( mu2 - mu1 ), Pls[5];

    P_ls( l, mu2, Pls );
    s  = ( a * mu2 + b ) * ( Pls[3] - Pls[1] ) + a * ( 2. * Pls[2] * tlp1 / tlm1 / tlp3 - Pls[4] / tlp3 - Pls[0] / tlm1 );
    P_ls( l, mu1, Pls );
    s -= ( a * mu1 + b ) * ( Pls[3] - Pls[1] ) + a * ( 2. * Pls[2] * tlp1 / tlm1 / tlp3 - Pls[4] / tlp3 - Pls[0] / tlm1 );
    s /= tlp1;
    return( s );
}
/*
============================================================
*/
double CalculateEpAve( MuList *MuL, double E ) {

    int i, j, notFirstTime = 0;
    EpList *EpL;
    double mu1, mu2, Sum = 0., EpSum = 0., dEp, Sum1, Sum2, EpSum1, EpSum2;

    for( i = 0; i < MuL->n; i++ ) {
        EpL = &(MuL->EpL[i]);
        mu2 = EpL->Mu;
        Sum2 = 0.;
        EpSum2 = 0.;
        for( j = 1; j < EpL->n; j++ ) {
            dEp = ( EpL->EpX[j].Ep - EpL->EpX[j-1].Ep );
            Sum2 += ( EpL->EpX[j].X + EpL->EpX[j-1].X ) * dEp;                                              // / 2. 
            EpSum2 += ( ( EpL->EpX[j].X + EpL->EpX[j-1].X ) * ( EpL->EpX[j].Ep + EpL->EpX[j-1].Ep ) + 
                          EpL->EpX[j].X * EpL->EpX[j].Ep + EpL->EpX[j-1].X * EpL->EpX[j-1].Ep ) * dEp;      // / 6.  
        }
        if( notFirstTime ) {
            Sum += ( Sum2 + Sum1 ) * ( mu2 - mu1 );                                                         // / 2. 
            EpSum += ( EpSum2 + EpSum1 ) * ( mu2 - mu1 );                                                   // / 2. 
        }
        notFirstTime = 1;
        mu1 = mu2;
        Sum1 = Sum2;
        EpSum1 = EpSum2;
    }
    if( Sum <= 0. ) {
        fprintf( stderr, "legendreExpand:CalculateEpAve: Sum = %e <= 0. for E = %e", Sum, E );
        exit( EXIT_FAILURE );
    }
    return( EpSum / ( 3. * Sum ) );
}
/*
============================================================
*/
void addToMuList( MuList *l, EpList *EpL, double Mu ) {

    int i;
    DummyList *d = (DummyList *) l;
    EpList *EpLNew;

    checkListSizeToAdd( d, sizeof( EpList ), 1, 200 );                                                      // AF7, Freed end of main. 
    EpLNew = &(l->EpL[l->n]);
    EpLNew->n = 0;
    EpLNew->nAllocated = 0;
    EpLNew->Mu = Mu;
    EpLNew->EpX = NULL;
    checkListSizeToAdd( (DummyList *) EpLNew, sizeof( EpX ), EpL->n, 1 );                                   // AF5, Freed in convertEMuEp_To_EEpl. 
    for( i = 0; i < EpL->n; i++ ) {
        EpLNew->EpX[i].Ep = EpL->EpX[i].Ep;
        EpLNew->EpX[i].X  = EpL->EpX[i].X;
    }
    EpLNew->n = EpL->n;
    l->n++;
}
/*
============================================================
*/
void addToEpList( EpList *l, double d1, double d2, double d3, double d4 ) {

    DummyList *d = (DummyList *) l;

    checkListSizeToAdd( d, sizeof( EpX ), 1, 1000 );                                                        // AF8, Freed in . 
    l->EpX[l->n].Ep = d3;
    l->EpX[l->n].X  = d4;
    l->n++;
}
/*
============================================================
*/
int mergeXXYLists( int nx, double *x,  int nxy, double *xy, int nmerge, double *mergexxy ) {

    int ix, ixy, n = 0;
    double *priorxy = NULL, *firstx = x;

    while( ( nx > 0 ) || ( nxy > 0 ) ) {
        if( ( nxy > 0 ) && ( ( nx == 0 ) || ( *xy <= *x ) ) ) {
            if( ( priorxy == NULL ) && ( firstx != x ) ) {
                *(mergexxy++) = *xy;
                *(mergexxy++) = 0.;
                n++;
            }
            if( ( nx > 0 ) && ( *xy == *x ) ) {
                x++;
                nx--;
            }
            priorxy = xy;
            *(mergexxy++) = (*xy++);
            *(mergexxy++) = (*xy++);
            nxy--;
            n++; }
        else {
            if( nxy == 0 ) {
                *(mergexxy++) = priorxy[0];
                *(mergexxy++) = priorxy[1];
                nxy--;
            }
            *(mergexxy++) = *x;
            if( ( priorxy == NULL ) || ( nxy < 0 ) ) {
                *(mergexxy++) = 0.; }
            else {
                *(mergexxy++) = ( priorxy[1] * ( xy[0] - *x ) + xy[1] * ( *x - priorxy[0] ) ) / ( xy[0] - priorxy[0] );
            }
            x++;
            nx--;
            n++;
        }
    }
    return( n );
}
/*
============================================================
*/
void checkListSizeToAdd( DummyList *l, size_t s, int nToAdd, int increment ) {

    if( ( l->n + nToAdd ) > l->nAllocated ) {
        if( nToAdd > increment ) {
            l->nAllocated += nToAdd; }
        else {
            l->nAllocated += increment;
        }
        l->v = myRealloc( l->v, l->nAllocated * s, "l->v" );
        if( l->v == NULL ) {
            fprintf( stderr, "legendreExpand:checkListSizeToAdd : realloc failed for size %lu", l->nAllocated * s );
            exit( EXIT_FAILURE );
        }
    }
}
/*
============================================================
*/
void *myMalloc( size_t s, char *m ) {

    void *p = malloc( s );

    if( p == NULL ) {
        fprintf( stderr, "legendreExpand:myMalloc %s", m );
        exit( EXIT_FAILURE );
    }
    myPrintf( "myMalloc:    %p\n", p );
    return( p );
}
/*
============================================================
*/
void *myRealloc( void *p, size_t s, char *m ) {

    if( p != NULL ) myPrintf( "myRealloc f: %p\n", p );
    p = realloc( p, s );
    if( p == NULL ) {
        fprintf( stderr, "legendreExpand:myRealloc %s", m );
        exit( EXIT_FAILURE );
    }
    myPrintf( "myRealloc:   %p\n", p );
    return( p );
}
/*
============================================================
*/
void myFree( void *p ) {

    myPrintf( "freeing    : %p\n", p );
    free( p );
}
/*
============================================================
*/
int getLine( FILE *f, char *fileName ) {

    int i;
    char *c,  buffer[bufferSize];
    
    fileLine++;
    c = fgets( buffer, bufferSize, f );
    if( c == NULL ) {
        if( feof( f ) ) {
            fileLine--;
            return( 1 );
        }
    } else {
        i = sscanf( buffer, "%le %le %le %le", &E, &Mu, &Ep, &prob );
        if( i != 4 ) {
            fprintf( stderr, "legendreExpand:getLine : sscanf failed at line %d : i = %d",fileLine,i );
            exit( EXIT_FAILURE );
        }
    }
    return( 0 );
}
/*
============================================================
*/
void printMsg( char *s ) {

    fprintf( stderr, "\nError from calcUpscatterKernelBoth: %s\n\n", s );
    exit( EXIT_FAILURE );
}
/*
============================================================
*/
void myPrintfNothing( char *s, ... ) {

    return;
}
/*
============================================================
*/
void calcKernel( double A, double T_MeV, double MinEValue,  double MaxEValue , double * gb, int ng, int * gMax) {

    int i, iMu, nMu = 40;
    char *e;
    double E, EMax, fE, x, y, Mu, P;
    double xMin, xMax, xxMax, MuMin = -1., MuMax = 1., dMu;
    double xsecRelMin = 1e-20;
    
    FILE *fp ;
    fp = fopen("upscatterEMuEp.out","w+");
    
    fE = pow( 10., 1. / 100. );
    dMu = ( MuMax - MuMin ) / nMu;

    EMax = 1000. * T_MeV / A;
    if( EMax < EmaxScalingFactor * T_MeV ) EMax = EmaxScalingFactor * T_MeV;    // Large massive target, still go to 3 time neutron v_thermal. 
    if( EMax > .1 ) EMax = .1;                                                  // No more than 0.1 MeV. 
    
    for( i = 1; i < ng; i++ ) {  if( (EMax > gb[i-1] ) && (EMax <= gb[i]) ) *gMax = i-2;  }   // set the Emax to a nearby group boundary so that the substitution will match up 
    if( *gMax < 0 ) *gMax = 0;
    EMax = fE * fE * gb[*gMax+1];       // Go a little above group so while loop below includes gb[*gMax+1].

    E = MinEValue;
    while( E <= EMax ) {
        y = E / T_MeV;
        x = 0.;
        if( A > 1 ) x = ( A - 1. ) / ( A + 1 );     // Obtained by setting alpha + beta = 0 at Mu = 0. 
        if( x < MinEValue / E ) x = MinEValue / E;
        xsecMin = xsecRelMin * xsec( y, x, 0., A );
        for( iMu = 0, Mu = MuMin; iMu <= nMu; iMu++, Mu += dMu ) {
            currentxPData = 0;
            if( fabs( Mu ) < 1e-14 ) Mu = 0.;
            if( Mu > 0.9999 ) Mu = 0.9999;
            xMin = 0.;
            if( ( Mu * Mu + A * A - 1. ) >= 0 ) xMin = ( Mu + sqrt( Mu * Mu + A * A - 1. ) ) / ( A + 1 );    // Obtained by setting alpha + beta = 0. 
            xMin *= xMin;                                           // The above solution is for sqrt( E' / E ), need to square to get x.
            if( xMin < MinEValue / E ) xMin = MinEValue / E;
            xMax = MaxEValue / E;
            x = 0.5 * ( xMin + xMax );
            while( 1 ) {
                P = xsec( y, x, Mu, A );
                if( P <= 0. ) {
                    xMax = x; }
                else if( P > 1.1 * xsecMin ) {
                    xMin = x; }
                else {
                    break;
                }
                x = 0.5 * ( xMin + xMax );
                if( x == xMin ) break;
            }
            xMax = x;
            xMin = MinEValue / E;
            if( ( P = xsec( y, xMin, Mu, A ) ) == 0 ) {
                xxMax = xMax;
                while( 1 ) {
                    x = 0.5 * ( xMin + xxMax );
                    P = xsec( y, x, Mu, A );
                    if( P <= 0 ) {
                        xMin = x; }
                    else if( P > 1.1 * xsecMin ) {
                        xxMax = x; }
                    else {
                        break;
                    }
                }
                xMin = x;
            }
            insertPoint( -1, xMin, xsec( y, xMin, Mu, A ) );
            insertPoint( 0, xMax, xsec( y, xMax, Mu, A ) );
            checkPoint( 0, 1, y, Mu, A );
            if( Mu > 0.9995 ) Mu = 1.;
            for( i = 0; i != -1; i = xPDatas[i].Next ) fprintf( fp, "%14.6e %11.3e %14.6e %14.6e\n", E, Mu, E *  xPDatas[i].x,  xPDatas[i].P );
            for( i = 0; i != -1; i = xPDatas[i].Next ) { insertEEMPpoint( E, E*xPDatas[i].x, Mu,  xPDatas[i].P ); }    
            
            clearPoints( );
        }
        E *= fE;
    }
    EEMPdata = realloc( EEMPdata, (currentEEMP+1) * sizeof( EEpMuP ) );

    fclose(fp);
    
}
/*
============================================================
*/
int checkPoint( int indexMin, int indexMax, double y, double Mu, double A ) {

    int i;
    double x = 0.5 * ( xPDatas[indexMin].x + xPDatas[indexMax].x ), P = xsec( y, x, Mu, A ), PMid, xx, xMin, xMax;

    if( P == 0 ) {
        xMin = xPDatas[indexMin].x;
        xMax = x;
        while( 1 ) {
            xx = 0.5 * ( xMin + xMax );
            P = xsec( y, xx, Mu, A );
            if( P == 0 ) {
                xMax = xx; }
            else if( P > 1.1 * xsecMin ) {
                xMin = xx; }
            else {
                break;
            }
        }
        i = insertPoint( indexMin, xx, P );
        checkPoint( indexMin, i, y, Mu, A );

        xMin = x;
        xMax = xPDatas[indexMax].x;
        while( 1 ) {
            xx = 0.5 * ( xMin + xMax );
            P = xsec( y, xx, Mu, A );
            if( P == 0 ) {
                xMin = xx; }
            else if( P > 1.1 * xsecMin ) {
                xMax = xx; }
            else {
                break;
            }
        }
        i = insertPoint( i, xx, P );
        checkPoint( i, indexMax, y, Mu, A );
        return( 0 );
    }
    if( P * fraction >= fabs( P - 0.5 * ( xPDatas[indexMin].P + xPDatas[indexMax].P ) ) ) {
        PMid = xsec( y, 0.5 * ( xPDatas[indexMin].x + x ), Mu, A );
        if(  PMid * fraction >= fabs( PMid - 0.5 * ( xPDatas[indexMin].P + P ) ) ) {
            PMid = xsec( y, 0.5 * ( x + xPDatas[indexMax].x ), Mu, A );
            if(  PMid * fraction >= fabs( PMid - 0.5 * ( P + xPDatas[indexMax].P ) ) ) return( 1 );
        }
    }
    i = insertPoint( indexMin, x, P );
    checkPoint( indexMin, i, y, Mu, A );
    checkPoint( i, indexMax, y, Mu, A );
    return( 1 );
}
/*
============================================================
*/
int insertPoint( int indexPrior, double x, double P ) {

    if( nxPData <= currentxPData ) {
        nxPData += 1000;
        xPDatas = realloc( xPDatas, nxPData * sizeof( xPData ) );
        if( xPDatas == NULL ) printMsg( "realloc memory for xPDatas" );
    }
    xPDatas[currentxPData].Index = currentxPData;
    xPDatas[currentxPData].x = x;
    xPDatas[currentxPData].P = P;
    if( currentxPData == 0 ) {
        xPDatas[currentxPData].Prior = -1;
        xPDatas[currentxPData].Next  = -1; }
    else {
        xPDatas[currentxPData].Next = xPDatas[indexPrior].Next;
        if( xPDatas[currentxPData].Next != -1 ) xPDatas[xPDatas[currentxPData].Next].Prior = currentxPData;
        xPDatas[currentxPData].Prior = indexPrior;
        xPDatas[indexPrior].Next = currentxPData;
    }
    currentxPData++;
    return( currentxPData - 1 );
}
/*
============================================================
*/
void insertEEMPpoint(  double E, double Ep, double Mu, double P ) {


    if( nEEMP <= currentEEMP ) {
        nEEMP += 1000;
        EEMPdata = realloc( EEMPdata, nEEMP * sizeof( EEpMuP ) );
        if( EEMPdata == NULL ) printMsg( "realloc memory for EEMPdata" );
    }
    EEMPdata[currentEEMP].Index = currentEEMP;
    EEMPdata[currentEEMP].E = E;
    EEMPdata[currentEEMP].Ep = Ep;
    EEMPdata[currentEEMP].Mu = Mu;
    EEMPdata[currentEEMP].P = P;
    //if (currentEEMP%100000 == 0 ) printf( "adding to EEMP : index=%d, current=%d, E=%14.6e Mu=%11.3e Ep=%14.6e Prob=%14.6e\n", EEMPdata[currentEEMP].Index, currentEEMP, E, Mu, Ep, P  ) ;  
    
    currentEEMP++;
    
}
/*
============================================================
*/
int getNextEEMPline(  ) {

    
    if( nEEMP < currentEEMP ) { return( 1 ); }
    if(EEMPdata[currentEEMP].Index != currentEEMP) {  
        // printf( "bad  EEMP access request: NumEEMPpts=%d, index=%d, current=%d, E=%14.6e Mu=%11.3e Ep=%14.6e Prob=%14.6e\n", nEEMP , EEMPdata[currentEEMP].Index, currentEEMP, E, Mu, Ep, prob  ) ; 
        return (1); 
    }
    
    E = EEMPdata[currentEEMP].E;
    Ep = EEMPdata[currentEEMP].Ep;
    Mu = EEMPdata[currentEEMP].Mu;
    prob = EEMPdata[currentEEMP].P;
    //if (currentEEMP%100000 == 0 ) printf( "getting from EEMP : index=%d, current=%d, E=%14.6e Mu=%11.3e Ep=%14.6e Prob=%14.6e\n", EEMPdata[currentEEMP].Index, currentEEMP, E, Mu, Ep, prob  ) ;  
    
    currentEEMP++; 
    
    return(0);
   
}
/*
============================================================
*/
void clearPoints( void ) {

    int i = 0;

    for( i = 0; i < currentxPData; i++ ) {
        xPDatas[i].Prior = 0;
        xPDatas[i].Next = 0;
        xPDatas[i].Index = 0;
        xPDatas[i].x = 0.;
        xPDatas[i].P = 0.;
    }
}
/*
============================================================
*/
double xsec( double y, double x, double Mu, double A ) {

    double sqrtx = sqrt( x ), xs, a = y * ( 1. + x - 2. * Mu * sqrtx ) / A, b = y * ( x - 1. );

    b = a + b;
    b *= b;
    xs = exp( -b / ( 4. * a ) );
    xs *= sqrtx / ( 2. * sqrt( 4. * Pi * a ) ) / T_MeV;
    if( xs < xsecMin ) xs = 0.;
    return( xs );
}
