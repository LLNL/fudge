/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";
static double yMin = -0.9, yMax = 0.4;

static int checkClipping( ptwXYPoints *data );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XY, *expXY, *mulXY;
    double x, expXYs[4];

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( XY = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i - .5, yMax + i - .1 ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkClipping( XY );
    ptwXY_neg( XY );
    errCount += checkClipping( XY );
    ptwXY_neg( XY );

    for( ; i < 2 * nSame + 1; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i - .5, yMin - i - .1 ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkClipping( XY );
    ptwXY_neg( XY );
    errCount += checkClipping( XY );
    ptwXY_neg( XY );

    if( ( status = ptwXY_clear( XY ) ) != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: clear, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    for( i = 0; i < 501; i++ ) {
        x = i * M_PI / 50;
        if( ( status = ptwXY_setValueAtX( XY, x, sin( x ) ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 3, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkClipping( XY );

    expXYs[0] = ptwXY_getXMin( XY );
    expXYs[1] = 0.;
    expXYs[2] = ptwXY_getXMax( XY );
    expXYs[3] = 1.;
    if( ( expXY = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 100, 10, 2, expXYs, &status, 0 ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: expXYs create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( expXY );
    if( ( status = ptwXY_exp( expXY, 1. ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_exp, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( expXY );
    if( ( mulXY = ptwXY_mul_ptwXY( XY, expXY, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY_mul_ptwXY, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkClipping( mulXY );

    ptwXY_free( XY );
    ptwXY_free( expXY );
    ptwXY_free( mulXY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkClipping( ptwXYPoints *data ) {

    int i, errCount = 0;
    ptwXYPoints *clipped, *u;
    ptwXYPoint *point;
    nfu_status status;
    double x1, y1, x2, y2, x, y, s;

    if( verbose ) {
        printf( "# yMin = %.14e\n", yMin );
        printf( "# yMax = %.14e\n", yMax );
    }
    printIfVerbose( data );
    if( ( clipped = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data clone, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    ptwXY_clip( clipped, yMin, yMax );
    printIfVerbose( clipped );

    if( ( u = ptwXY_union( clipped, data, &status, ptwXY_union_fill ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: u, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < u->length; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( u, i );
        if( point->y < yMin ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g < yMin = %g", __FILE__, point->x, point->y, yMin );
            errCount++; }
        else if( point->y > yMax ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g > yMax = %g", __FILE__, point->x, point->y, yMax );
            errCount++;
        }
        if( ( status = ptwXY_getValueAtX( data, point->x, &y ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_getValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
        if( y < yMin ) {
            if( point->y != yMin ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMin = %g", __FILE__, point->x, y, point->y, yMin );
                errCount++;
            }
        }
        else if( y > yMax ) {
            if( point->y != yMax ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMax = %g", __FILE__, point->x, y, point->y, yMax );
                errCount++;
            }
        }
    }

    for( i = 0; i < data->length; i++ ) {       /* test insures that points are added at yMin and yMax crossings. */
        point = ptwXY_getPointAtIndex_Unsafely( data, i );
        x2 = point->x;
        y2 = point->y;
        if( i > 0 ) {
            if( ( ( y1 - yMin ) * ( y2 - yMin ) ) < 0 ) {
                x = ( x1 * ( y2 - yMin ) + x2 * ( yMin - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( clipped, x, &y ) ) != nfu_Okay )
                        nfu_printErrorMsg( "ERROR %s: ptwXY_getValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
                s = 0.5 * ( fabs( y ) + fabs( yMin ) );
                if( fabs( y - yMin ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMin = %g", __FILE__, i, x, y, yMin );
                    errCount++;
                }
            }
            if( ( ( y1 - yMax ) * ( y2 - yMax ) ) < 0 ) {
                x = ( x1 * ( y2 - yMax ) + x2 * ( yMax - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( clipped, x, &y ) ) != nfu_Okay )
                        nfu_printErrorMsg( "ERROR %s: ptwXY_getValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
                s = 0.5 * ( fabs( y ) + fabs( yMax ) );
                if( fabs( y - yMax ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMax = %g", __FILE__, i, x, y, yMax );
                    errCount++;
                }
            }
        }
        x1 = x2;
        y1 = y2;
    }
    ptwXY_free( u );
    ptwXY_free( clipped );

    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
