/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <math.h>
#include <float.h>

#include "ptwXY.h"

typedef nfu_status (*interpolation_func)( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth );

static double ptwXY_flatInterpolationToLinear_eps( double px, double eps );
static nfu_status ptwXY_toOtherInterpolation2( ptwXYPoints *desc, ptwXYPoints *src, interpolation_func func );
static nfu_status ptwXY_LogLogToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth );
static nfu_status ptwXY_LinLogToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth );
static nfu_status ptwXY_LogLinToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth );
/*
************************************************************
*/
nfu_status ptwXY_interpolatePoint( ptwXY_interpolation interpolation, double x, double *y, double x1, double y1, double x2, double y2 ) {

    nfu_status status = nfu_Okay;

    if( interpolation == ptwXY_interpolationOther ) return( nfu_otherInterpolation );
    if( ( x1 > x2 ) || ( x < x1 ) || ( x > x2 ) ) return( nfu_invalidInterpolation );
    if( y1 == y2 ) {
        *y = y1; }
    else if( x1 == x2 ) {
        *y = 0.5 * ( y1  + y2 ); }
    else if( x == x1 ) {
        *y = y1; }
    else if( x == x2 ) {
        *y = y2; }
    else {
        switch( interpolation ) {
        case ptwXY_interpolationLinLin :
            *y = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
            break;
        case ptwXY_interpolationLogLin :
            if( ( x <= 0. ) || ( x1 <= 0. ) || ( x2 <= 0. ) ) return( nfu_invalidInterpolation );
            *y = ( y1 * log( x2 / x ) + y2 * log( x / x1 ) ) / log( x2 / x1 );
            break;
        case ptwXY_interpolationLinLog :
            if( ( y1 <= 0. ) || ( y2 <= 0. ) ) return( nfu_invalidInterpolation );
            *y = exp( ( log( y1 ) * ( x2 - x ) + log( y2 ) * ( x - x1 ) ) / ( x2 - x1 ) );
            break;
        case ptwXY_interpolationLogLog :
            if( ( x <= 0. ) || ( x1 <= 0. ) || ( x2 <= 0. ) ) return( nfu_invalidInterpolation );
            if( ( y1 <= 0. ) || ( y2 <= 0. ) ) return( nfu_invalidInterpolation );
            *y = exp( ( log( y1 ) * log( x2 / x ) + log( y2 ) * log( x / x1 ) ) / log( x2 / x1 ) );
            break;
        case ptwXY_interpolationFlat :
            *y = y1;
            break;
        default :
            status = nfu_invalidInterpolation;
        }
    }
    return( status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_flatInterpolationToLinear( ptwXYPoints *ptwXY, double lowerEps, double upperEps, nfu_status *status ) {

    int64_t i, length;
    double x;
    ptwXYPoints *n;
    ptwXYPoint *p1 = NULL, *p2 = NULL, *p3;

#define minEps 5e-16

    if( ( *status = ptwXY_simpleCoalescePoints( ptwXY ) ) != nfu_Okay ) return( NULL );
    *status = nfu_invalidInterpolation;
    if( ptwXY->interpolation != ptwXY_interpolationFlat ) return( NULL );
    *status = nfu_badInput;
    if( ( lowerEps < 0 ) || ( upperEps < 0 ) || ( ( lowerEps == 0 ) && ( upperEps == 0 ) ) ) return( NULL );
    if( ( lowerEps != 0 ) && ( lowerEps < minEps ) ) lowerEps = minEps;
    if( ( upperEps != 0 ) && ( upperEps < minEps ) ) upperEps = minEps;

    length = ptwXY->length * ( 1 + ( lowerEps == 0 ? 0 : 1 ) + ( lowerEps == 0 ? 0 : 1 ) );
    if( ( n = ptwXY_new( ptwXY_interpolationLinLin, ptwXY->biSectionMax, ptwXY->accuracy, length, ptwXY->overflowLength, status, ptwXY->userFlag ) ) == NULL ) return( NULL );

    p3 = ptwXY->points;
    if( ptwXY->length > 0 ) ptwXY_setValueAtX( n, p3->x, p3->y );
    for( i = 0; i < ptwXY->length; i++, p3++ ) {
        if( i > 1 ) {
            if( lowerEps > 0 ) {
                x = ptwXY_flatInterpolationToLinear_eps( p2->x, -lowerEps );
                if( x <= p1->x ) goto xErr;
                if( ( *status = ptwXY_setValueAtX( n, x, p1->y ) ) != nfu_Okay ) goto Err;
            }
            if( lowerEps == 0 ) if( ( *status = ptwXY_setValueAtX( n, p2->x, p1->y ) ) != nfu_Okay ) goto Err;
            if( upperEps == 0 ) if( ( *status = ptwXY_setValueAtX( n, p2->x, p2->y ) ) != nfu_Okay ) goto Err;
            if( upperEps > 0 ) {
                x = ptwXY_flatInterpolationToLinear_eps( p2->x, upperEps );
                if( x >= p3->x ) goto xErr;
                if( ( *status = ptwXY_setValueAtX( n, x, p2->y ) ) != nfu_Okay ) goto Err;
            }
        }
        p1 = p2;
        p2 = p3;
    }
    if( ptwXY->length > 1 ) {
        if( ( lowerEps != 0 ) && ( p1->y != p2->y ) ) {
            x = ptwXY_flatInterpolationToLinear_eps( p2->x, -lowerEps );
            if( x <= p1->x ) goto xErr;
            if( ( *status = ptwXY_setValueAtX( n, x, p1->y ) ) != nfu_Okay ) goto Err;
        }
        if( ( *status = ptwXY_setValueAtX( n, p2->x, p2->y ) ) != nfu_Okay ) goto Err;
    }

    return( n );

xErr:
    *status = nfu_XNotAscending;
Err:
    ptwXY_free( n );
    return( NULL );

#undef minEps
}
/*
************************************************************
*/
static double ptwXY_flatInterpolationToLinear_eps( double px, double eps ) {

    double x;

    if( px < 0 ) {
        x = ( 1 - eps ) * px; }
    else if( px > 0 ) {
        x = ( 1 + eps ) * px; }
    else {
        x = eps;
    }
    return( x );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_toOtherInterpolation( ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, double accuracy, nfu_status *status ) {

    ptwXYPoints *n;
    interpolation_func func = NULL;

    if( ( *status = ptwXY->status ) != nfu_Okay ) return( NULL );
/* ?????? need to test interpolation also */
    if( ptwXY->interpolation == interpolation ) {
        *status = nfu_Okay;
        return( ptwXY_clone( ptwXY, status ) ); }
    else {
        switch( ptwXY->interpolation ) {
        case ptwXY_interpolationLogLog :
            func = ptwXY_LogLogToLinLin; break;
        case ptwXY_interpolationLinLog :
            func = ptwXY_LinLogToLinLin; break;
        case ptwXY_interpolationLogLin :
            func = ptwXY_LogLinToLinLin; break;
        case ptwXY_interpolationLinLin :        /* Stops compilers from complaining. */
        case ptwXY_interpolationFlat :
        case ptwXY_interpolationOther :
            break;
        }
    }
    *status = nfu_unsupportedInterpolationConversion;
    if( func == NULL ) return( NULL );

    *status = nfu_Okay;
    if( ( n = ptwXY_clone( ptwXY, status ) ) == NULL ) return( NULL );
    n->interpolation = interpolation;
    if( accuracy < ptwXY->accuracy ) accuracy = ptwXY->accuracy;
    n->accuracy = accuracy;

    *status = ptwXY_toOtherInterpolation2( n, ptwXY, func );
    if( *status != nfu_Okay ) n = ptwXY_free( n );
    return( n );
}
/*
************************************************************
*/
static nfu_status ptwXY_toOtherInterpolation2( ptwXYPoints *desc, ptwXYPoints *src, interpolation_func func ) {

    nfu_status status;
    int64_t i;
    double x1, y1, x2, y2;

    if( ( status = ptwXY_simpleCoalescePoints( src ) ) != nfu_Okay ) return( status );

    x1 = src->points[0].x;
    y1 = src->points[0].y;
    for( i = 1; i < src->length; i++ ) {
        x2 = src->points[i].x;
        y2 = src->points[i].y;
        if( ( x1 != x2 ) && ( y1 != y2 ) ) {
            if( ( status = func( desc, x1, y1, x2, y2, 0 ) ) != nfu_Okay ) break;
        }
        x1 = x2;
        y1 = y2;
    }
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_LogLogToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x, y, u, u2 = x2 / x1, v2 = y2 / y1, logYXs, logXs = log( u2 ), logYs = log( v2 ), vLin, vLog, w;

    logYXs = logYs / logXs;

    if( depth > 16 ) return( nfu_Okay );
    if( fabs( logYXs  - 1 ) < 1e-5 ) {
        u = 0.5 * ( 1 + u2 );
        w = ( logYXs  - 1 ) * logXs;
        vLog = u * ( 1. + w * ( 1 + 0.5 * w ) ); }
    else {
        u = logYXs * ( u2 - v2 ) / ( ( 1 - logYXs ) * ( v2 - 1 ) );
        vLog = pow( u, logYXs );
    }
    vLin = ( u2 - u + v2 * ( u - 1 ) ) / ( u2 - 1 );
    if( fabs( vLog - vLin ) <= ( vLog * desc->accuracy ) ) return( status );
    x = x1 * u;
    y = y1 * vLog;
    if( ( status = ptwXY_setValueAtX( desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LogLogToLinLin( desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LogLogToLinLin( desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_LinLogToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x, y, logYs = log( y2 / y1 ), yLinLin;

    if( depth > 16 ) return( nfu_Okay );
    x = ( x2 - x1 ) / ( y2 - y1 ) * ( ( y2 - y1 ) / logYs - y1 ) + x1;
    y = y1 * exp( logYs / ( x2 - x1 ) * ( x - x1 ) );
    yLinLin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
    if( fabs( y - yLinLin ) <= ( y * desc->accuracy ) ) return( status );
    if( ( status = ptwXY_setValueAtX( desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LinLogToLinLin( desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LinLogToLinLin( desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
static nfu_status ptwXY_LogLinToLinLin( ptwXYPoints *desc, double x1, double y1, double x2, double y2, int depth ) {

    nfu_status status = nfu_Okay;
    double x = sqrt( x2 * x1 ), y, logXs = log( x2 / x1 ), yLinLin;

    if( depth > 16 ) return( nfu_Okay );
#if 0 /* The next line is very unstable at determineing x. Initial x must be chosen better. */
    x = ( y1 * x2 - y2 * x1 ) / ( y1 * logXs + ( y2 - y1 ) * ( log( x / x1 ) - 1 ) );
#endif
    y = ( y2 - y1 ) * log( x / x1 ) / logXs + y1;
    yLinLin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 );
    if( fabs( y - yLinLin ) <= ( y * desc->accuracy ) ) return( status );
    if( ( status = ptwXY_setValueAtX( desc, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_LogLinToLinLin( desc, x1, y1, x, y, depth + 1 ) ) != nfu_Okay ) return( status );
    return( ptwXY_LogLinToLinLin( desc, x, y, x2, y2, depth + 1 ) );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_toUnitbase( ptwXYPoints *ptwXY, nfu_status *status ) {

    int64_t i;
    ptwXYPoints *n;
    ptwXYPoint *p;
    double xMin, xMax, dx, inverseDx;

    *status = nfu_tooFewPoints;
    if( ptwXY->length < 2 ) return( NULL );
    if( ( n = ptwXY_clone( ptwXY, status ) ) == NULL ) return( NULL );

    xMin = n->points[0].x;
    xMax = n->points[n->length-1].x;
    dx = xMax - xMin;
    inverseDx = 1. / dx;
    for( i = 0, p = n->points; i < n->length; i++, p++ ) {
        p->x = ( p->x - xMin ) * inverseDx;
        p->y = p->y * dx;
    }
    n->points[n->length-1].x = 1.;                          /* Make sure last point is realy 1. */
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_fromUnitbase( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status ) {

    int64_t i, length;
    ptwXYPoints *n;
    ptwXYPoint *p, *p2;
    double dx, inverseDx, xLast = 0.;

    *status = nfu_tooFewPoints;
    if( ptwXY->length < 2 ) return( NULL );
    if( ( n = ptwXY_clone( ptwXY, status ) ) == NULL ) return( NULL );

    dx = xMax - xMin;
    inverseDx = 1. / dx;
    length = n->length;
    for( i = 0, p2 = p = n->points; i < length; ++i, ++p ) {
        p2->x = p->x * dx + xMin;
        if( i > 0 ) {
            if( fabs( p2->x - xLast ) <= 10. * DBL_EPSILON * ( fabs( p2->x ) + fabs( xLast ) ) ) {
                --(n->length);
                continue;
            }
        }
        p2->y = p->y * inverseDx;
        xLast = p2->x;
        ++p2;
    }
    n->points[n->length-1].x = xMax;                          /* Make sure last point is realy xMax. */
    return( n );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_unitbaseInterpolate( double w, double w1, ptwXYPoints *ptwXY1, double w2, ptwXYPoints *ptwXY2, nfu_status *status ) {
/* 
*   Should we not be checking the interpolation members???????
*/
    int64_t i;
    ptwXYPoints *n1 = NULL, *n2 = NULL, *a = NULL, *r = NULL;
    ptwXYPoint *p;
    double f, g, xMin, xMax;

    *status = nfu_XOutsideDomain;
    if( w <= w1 ) {
        if( w < w1 ) return( NULL );
        return( ptwXY_clone( ptwXY1, status ) );
    }
    if( w >= w2 ) {
        if( w > w2 ) return( NULL );
        return( ptwXY_clone( ptwXY2, status ) );
    }
    if( ( n1 = ptwXY_toUnitbase( ptwXY1, status ) ) == NULL ) return( NULL );
    if( ( n2 = ptwXY_toUnitbase( ptwXY2, status ) ) == NULL ) goto Err;
    f = ( w - w1 ) / ( w2 - w1 );
    g = 1. - f;
    for( i = 0, p = n1->points; i < n1->length; i++, p++ ) p->y *= g;
    for( i = 0, p = n2->points; i < n2->length; i++, p++ ) p->y *= f;
    if( ( a = ptwXY_add_ptwXY( n1, n2, status ) ) == NULL ) goto Err;

    xMin = g * ptwXY1->points[0].x + f * ptwXY2->points[0].x;
    xMax = g * ptwXY1->points[ptwXY1->length-1].x + f * ptwXY2->points[ptwXY2->length-1].x;
    if( ( r = ptwXY_fromUnitbase( a, xMin, xMax, status ) ) == NULL ) goto Err;
    ptwXY_free( n1 );
    ptwXY_free( n2 );
    ptwXY_free( a );
    return( r );

Err:
    if( n1 != NULL ) ptwXY_free( n1 );
    if( n2 != NULL ) ptwXY_free( n2 );
    if( a  != NULL ) ptwXY_free( a );
    return( NULL );
}
