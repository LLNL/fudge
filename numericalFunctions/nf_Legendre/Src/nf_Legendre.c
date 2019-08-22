/*
# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>
*/

#include "nf_Legendre.h"

struct nf_Legendre_from_ptwXY_callback_s {
    int l;
    double mu1, mu2, f1, f2;
};

static nfu_status nf_Legendre_to_ptwXY2( double mu, double *P, void *argList );
static nfu_status nf_Legendre_from_ptwXY_callback( double mu, double *f, void *argList );
/*
************************************************************
*/
nf_Legendre *nf_Legendre_new( int initialSize, int maxOrder, double *Cls, nfu_status *status ) {

    int l;
    nf_Legendre *Legendre = (nf_Legendre *) nfu_malloc( sizeof( nf_Legendre ) );

    *status = nfu_mallocError;
    if( Legendre == NULL ) return( NULL );
    if( ( *status = nf_Legendre_setup( Legendre, initialSize, maxOrder ) ) != nfu_Okay ) {
        nfu_free( Legendre );
        return( NULL );
    }
    for( l = 0; l <= Legendre->maxOrder; l++ ) Legendre->Cls[l] = Cls[l];
    return( Legendre );
}
/*
************************************************************
*/
nfu_status nf_Legendre_setup( nf_Legendre *Legendre, int initialSize, int maxOrder ) {

    memset( Legendre, 0, sizeof( nf_Legendre ) );
    if( maxOrder < 0 ) maxOrder = -1;
    if( maxOrder > nf_Legendre_maxMaxOrder ) maxOrder = nf_Legendre_maxMaxOrder;
    Legendre->maxOrder = maxOrder;
    if( initialSize < ( maxOrder + 1 ) ) initialSize = maxOrder + 1;
    return( nf_Legendre_reallocateCls( Legendre, initialSize, 0 ) );
}
/*
************************************************************
*/
nfu_status nf_Legendre_release( nf_Legendre *Legendre ) {

    if( Legendre->allocated > 0 ) nfu_free( Legendre->Cls );
    memset( Legendre, 0, sizeof( nf_Legendre ) );
    return( nfu_Okay );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_free( nf_Legendre *Legendre ) {

    nf_Legendre_release( Legendre );
    nfu_free( Legendre );
    return( NULL );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_clone( nf_Legendre *nfL, nfu_status *status ) {

    return( nf_Legendre_new( 0, nfL->maxOrder, nfL->Cls, status ) );
}
/*
************************************************************
*/
nfu_status nf_Legendre_reallocateCls( nf_Legendre *Legendre, int size, int forceSmallerResize ) {

    nfu_status status = nfu_Okay;

    if( size < nf_Legendre_minMaxOrder ) size = nf_Legendre_minMaxOrder;
    if( size > ( nf_Legendre_maxMaxOrder + 1 ) ) size = nf_Legendre_maxMaxOrder + 1;
    if( size != Legendre->allocated ) {
        if( size > Legendre->allocated ) {
            Legendre->Cls = (double *) nfu_realloc( size * sizeof( double ), Legendre->Cls ); }
        else {
            if( size < ( Legendre->maxOrder + 1 ) ) size = Legendre->maxOrder + 1;
            if( ( Legendre->allocated > 2 * size ) || forceSmallerResize ) {
                    Legendre->Cls = (double *) nfu_realloc( size * sizeof( double ), Legendre->Cls ); } 
            else {
                size = Legendre->allocated;
            }
        }
        if( Legendre->Cls == NULL ) {
            size = 0;
            status = nfu_mallocError;
        }
        Legendre->allocated = size;
    }
    return( status );
}
/*
************************************************************
*/
int nf_Legendre_maxOrder( nf_Legendre *Legendre ) {

    return( Legendre->maxOrder );
}
/*
************************************************************
*/
int nf_Legendre_allocated( nf_Legendre *Legendre ) {

    return( Legendre->allocated );
}
/*
************************************************************
*/
double nf_Legendre_getCl( nf_Legendre *Legendre, int l, nfu_status *status ) {

    *status = nfu_Okay;
    if( ( l < 0 ) || ( l > Legendre->maxOrder ) ) {
        *status = nfu_badIndex;
        return( 0. );
    }
    return( Legendre->Cls[l] );
}
/*
************************************************************
*/
nfu_status nf_Legendre_setCl( nf_Legendre *Legendre, int l, double Cl ) {

    nfu_status status;

    if( ( l < 0 ) || ( l > ( Legendre->maxOrder + 1 ) ) ) return( nfu_badIndex );
    if( Legendre->allocated <= l ) {
        if( ( status = nf_Legendre_reallocateCls( Legendre, l + nf_Legendre_sizeIncrement, 0 ) ) != nfu_Okay ) return( status );
    }
    if( l > Legendre->maxOrder ) Legendre->maxOrder = l;
    Legendre->Cls[l] = Cl;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status nf_Legendre_normalize( nf_Legendre *Legendre ) {

    int l;
    double norm;

    if( Legendre->maxOrder >= 0 ) {
        if( ( norm = Legendre->Cls[0] ) == 0 ) return( nfu_divByZero );
        for( l = 0; l <= Legendre->maxOrder; l++ ) Legendre->Cls[l] /= norm;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
double nf_Legendre_evauluateAtMu( nf_Legendre *Legendre, double mu, nfu_status *status ) {

    int l;
    double P = 0.;

    *status = nfu_XOutsideDomain;
    if( ( mu >= -1. ) && ( mu <= 1. ) ) {
        *status = nfu_Okay;
        for( l = 0; l <= Legendre->maxOrder; l++ ) P += ( l + 0.5 ) * Legendre->Cls[l] * nf_Legendre_PofL_atMu( l, mu );
    }
    return( P );
}
/*
************************************************************
*/
double nf_Legendre_PofL_atMu( int l, double mu ) {

    int l_, twoL_plus1;
    double Pl_minus1, Pl, Pl_plus1;

    if( l == 0 ) {
        return( 1. ); }
    else if( l == 1 ) {
        return( mu ); }
/*
    else if( l <= 9 ) {
        double mu2 = mu * mu;
        if ( l == 2 ) {
            return(                1.5 * mu2 - 0.5 ); }
        else if( l == 3 ) {
            return(                2.5 * mu2 -          1.5 ) * mu; }
        else if( l == 4 ) {
            return(              4.375 * mu2 -         3.75 ) * mu2 +       0.375; }
        else if( l == 5 ) {
            return( (            7.875 * mu2 -         8.75 ) * mu2 +       1.875 ) * mu; }
        else if( l == 6 ) {
            return( (          14.4375 * mu2 -      19.6875 ) * mu2 +      6.5625 ) * mu2 -      0.3125; }
        else if( l == 7 ) {
            return( ( (        26.8125 * mu2 -      43.3125 ) * mu2 +     19.6875 ) * mu2 -      2.1875 ) * mu; }
        else if( l == 8 ) {
            return( ( (     50.2734375 * mu2 -     93.84375 ) * mu2 +   54.140625 ) * mu2 -     9.84375 ) * mu2 + 0.2734375; }
        else {
            return( ( ( (   94.9609375 * mu2 -    201.09375 ) * mu2 +  140.765625 ) * mu2 -    36.09375 ) * mu2 + 2.4609375 ) * mu;
        }
    }
*/

    Pl = 0.;
    Pl_plus1 = 1.;
    for( l_ = 0, twoL_plus1 = 1; l_ < l; l_++, twoL_plus1 += 2 ) {
        Pl_minus1 = Pl;
        Pl = Pl_plus1;
        Pl_plus1 = ( twoL_plus1 * mu * Pl - l_ * Pl_minus1 ) / ( l_ + 1 );
    }
    return( Pl_plus1 );
}
/*
************************************************************
*/
ptwXYPoints *nf_Legendre_to_ptwXY( nf_Legendre *Legendre, double accuracy, int biSectionMax, int checkForRoots, nfu_status *status ) {

    int i, n = 1;
    double dx, xs[1000];
    void *argList = (void *) Legendre;

    *status = nfu_Okay;
    xs[0] = -1;
    if( Legendre->maxOrder > 1 ) {
        n = Legendre->maxOrder - 1;
        if( n > 249 ) n = 249;
        n = 4 * n + 1;
        dx = 2. / n;
        for( i = 1; i < n; i++ ) xs[i] = xs[i-1] + dx;
    }
    xs[n] = 1.;
    return( ptwXY_createFromFunction( n + 1, xs, nf_Legendre_to_ptwXY2, (void *) argList, accuracy, checkForRoots, biSectionMax, status ) );
}
/*
************************************************************
*/
static nfu_status nf_Legendre_to_ptwXY2( double mu, double *P, void *argList ) {

    nfu_status status;      /* Set by nf_Legendre_evauluateAtMu. */

    *P = nf_Legendre_evauluateAtMu( (nf_Legendre *) argList, mu, &status );
    return( status );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_from_ptwXY( ptwXYPoints *ptwXY, int maxOrder, nfu_status *status ) {

    int l, i, n = (int) ptwXY_length( ptwXY );
    nf_Legendre *Legendre;
    double mu1, mu2, f1, f2, Cl, Cls[1] = { 0 }, integral;
    struct nf_Legendre_from_ptwXY_callback_s argList;

    if( ( *status = ptwXY_getStatus( ptwXY ) ) != nfu_Okay ) return( NULL );

    ptwXY_getXYPairAtIndex( ptwXY, 0, &mu1, &f1 );
    if( mu1 < -1 ) {
        *status = nfu_XOutsideDomain;
        return( NULL );
    }

    ptwXY_getXYPairAtIndex( ptwXY, 0, &mu2, &f2 );
    if( mu2 > 1 ) {
        *status = nfu_XOutsideDomain;
        return( NULL );
    }

    if( ( Legendre = nf_Legendre_new( maxOrder + 1, -1, Cls, status ) ) == NULL ) return( NULL );

    if( maxOrder > nf_Legendre_maxMaxOrder ) maxOrder = nf_Legendre_maxMaxOrder;
    for( l = 0; l <= maxOrder; l++ ) {
        ptwXY_getXYPairAtIndex( ptwXY, 0, &mu1, &f1 );
        argList.l = l;
        for( i = 1, Cl = 0; i < n; i++ ) {
            ptwXY_getXYPairAtIndex( ptwXY, i, &mu2, &f2 );
            argList.mu1 = mu1;
            argList.f1 = f1;
            argList.mu2 = mu2;
            argList.f2 = f2;
            if( ( *status = nf_Legendre_GaussianQuadrature( l + 1, mu1, mu2, nf_Legendre_from_ptwXY_callback, (void *) &argList, &integral ) ) != nfu_Okay )
                goto err;
            Cl += integral;
            mu1 = mu2;
            f1 = f2;
        }
        if( ( *status = nf_Legendre_setCl( Legendre, l, Cl ) ) != nfu_Okay ) goto err;
    }
    return( Legendre );

err:
    nf_Legendre_free( Legendre );
    return( NULL );    
}
/*
************************************************************
*/
static nfu_status nf_Legendre_from_ptwXY_callback( double mu, double *f, void *argList ) {

    struct nf_Legendre_from_ptwXY_callback_s *args = (struct nf_Legendre_from_ptwXY_callback_s *) argList;

    *f = ( args->f1 * ( args->mu2 - mu ) + args->f2 * ( mu - args->mu1 ) ) / ( args->mu2 - args->mu1 );
    *f *= nf_Legendre_PofL_atMu( args->l, mu );
    return( nfu_Okay );
}
