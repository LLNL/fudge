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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>

#define size 1001

static int verbose = 0;

double getDouble( const char const *s );
ptwXYPoints *randomUV( double accuracy, int biSectionMax );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int doRandom = 0, iarg, echo = 0, biSectionMax = 10;
    int64_t i, j, n;
    ptwXYPoints *u, *y, *e, *p;
    double xyPoints[2*2], x, dx, z, accuracy = 1e-3;
    ptwXYPoint *xy1, *xy2;
    nfu_status status;
    FILE *ff;
    char fmt[] = "%.14e %.14e\n";

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-r", argv[iarg] ) == 0 ) {
            doRandom = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );
    
    nfu_setMemoryDebugMode( 0 );

    xyPoints[0] = 0.0;
    xyPoints[1] = 1.0;
    xyPoints[2] = 1.0;
    xyPoints[3] = -0.2;
    if( argc == 5 ) {
        xyPoints[0] = getDouble( argv[1] );
        xyPoints[1] = getDouble( argv[2] );
        xyPoints[2] = getDouble( argv[3] );
        xyPoints[3] = getDouble( argv[4] );
    }

    ff = fopen( "info.dat", "w" );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# biSectionMax = %d\n", biSectionMax );
    fclose( ff );

    if( doRandom ) {
        u = randomUV( accuracy, biSectionMax ); }
    else {
        if( ( u = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, 2, xyPoints, &status, 0 ) ) == NULL ) 
            printMsg( "u creation: status = %d: %s", status, nfu_statusMessage( status ) );
    }
    if( ( y = ptwXY_clone( u, &status ) ) == NULL ) printMsg("y creation: status = %d: %s", status, nfu_statusMessage( status ));

    ff = fopen( "curve_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) u->length );
    ptwXY_simpleWrite( u, ff, fmt );
    fclose( ff );

    if( ( status = ptwXY_exp( y, 1 ) ) != nfu_Okay ) printMsg( "exp( z ): status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_simpleCoalescePoints( y ) ) != nfu_Okay ) printMsg( "coalescing z: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "exp_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) y->length );
    ptwXY_simpleWrite( y, ff, fmt );
    fclose( ff );

    if( ( e = ptwXY_new( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, &status, 0 ) ) == NULL ) printMsg( "e creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );
    n = ptwXY_length( u );
    xy1 = ptwXY_getPointAtIndex( u, 0 );
    for( i = 1; i < n; i++ ) {
        xy2 = ptwXY_getPointAtIndex( u, i );
        x = xy1->x;
        dx = ( xy2->x - xy1->x ) / 20;
        for( j = 0; j < 20; j++ ) {
            ptwXY_getValueAtX( u, x, &z );
            ptwXY_setValueAtX( e, x, exp( z ) );
            x += dx;
        }
        xy1 = xy2;
    }
    ptwXY_getValueAtX( u, xy2->x, &z );
    ptwXY_setValueAtX( e, x, exp( z ) );
    ff = fopen( "exactExp_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) e->length );
    ptwXY_simpleWrite( e, ff, fmt );
    fclose( ff );

    if( ( status = ptwXY_exp( u, -1 ) ) != nfu_Okay ) printMsg( "exp( -y ): status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( p = ptwXY_mul2_ptwXY( u, y, &status ) ) == NULL ) printMsg( "exp( y ) * exp( -y ): status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "exp_u_times_exp_minusU.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) p->length );
    ptwXY_simpleWrite( p, ff, fmt );
    fclose( ff );

    ptwXY_free( y );
    ptwXY_free( u );
    ptwXY_free( e );
    ptwXY_free( p );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
double getDouble( const char const *s ) {

    double d;
    char *e;

    errno = 0;
    d = strtod( s, &e );
    if( ( *e != 0 ) || ( errno != 0 ) ) printMsg( "could not convert '%s' to double, err = %d, e = %s", s, errno, e );
    return( d );
}
/*
****************************************************************
*/
ptwXYPoints *randomUV( double accuracy, int biSectionMax ) {

    int64_t i;
    double x, y;
    nfu_status status;
    ptwXYPoints *f;

    if( ( f = ptwXY_new( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, &status, 0 ) ) == NULL ) printMsg( "f new: status = %d: %s",
            status, nfu_statusMessage( status ) );
    for( i = 0, x = 0, y = 0; i < size; i++ ) {
        x += drand48( );
        y += drand48( ) - 0.5;
        if( ( status = ptwXY_setValueAtX( f, x, y ) ) != nfu_Okay ) printMsg( "ptwXY_setValueAtX( f, x, y ) failed: status = %d, %s", status,
            nfu_statusMessage( status ) );
    }
    return( f );
}
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
