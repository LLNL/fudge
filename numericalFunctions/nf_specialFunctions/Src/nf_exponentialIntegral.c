/*********************************************************************
   Returns the exponential integral function

   E_n(x) = int_1^infinity e^( -x * t ) / t^n dt,     for x > 0.

   C.A. Bertulani        May/15/2000
*********************************************************************/

#include "nf_specialFunctions.h"

#define EULER 0.57721566490153286   /* Euler's constant gamma */
#define MAXIT 100                   /* Maximum allowed number of iterations. */
#define FPMIN 1.0e-300              /* close to the smallest representable floting-point number. */
#define EPS 1.0e-15                 /* Desired relative error, not smaller than the machine precision. */

/*
************************************************************
*/
nfu_status nf_exponentialIntegral( statusMessageReporting *smr, int n, double x, double *value ) {

    int i, ii, nm1;
    double a, b, c, d, del, fact, h, psi;
    double ans = 0.0;
    nfu_status status = nfu_Okay;

    *value = 0.;

    if( !isfinite( x ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid input (x = %e).", x );
        return( nfu_badInput );
    }

    nm1 = n - 1;
    if( ( n < 0 ) || ( x < 0.0 ) || ( ( x == 0.0 ) && ( ( n == 0 ) || ( n == 1 ) ) ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid input (n = %d, x = %e).", n, x );
        return( nfu_badInput ); }
    else {
        if( n == 0 ) {
            ans = exp( -x ) / x; }                  /* Special case */
        else if( x == 0.0 ) {
            ans = 1.0 / nm1; }                      /* Another special case */
        else if( x > 1.0 ) {                        /* Lentz's algorithm */
            b = x + n;
            c = 1.0 / FPMIN;
            d = 1.0 / b;
            h = d;
            for( i = 1; i <= MAXIT; i++ ) {
                a = -i * ( nm1 + i );
                b += 2.0;
                d = 1.0 / ( a * d + b );            /* Denominators cannot be zero */
                c = b + a / c;
                del = c * d;
                h *= del;
                if( fabs( del - 1.0 ) < EPS ) {
                    *value = h * exp( -x );
                    return( nfu_Okay );
                }
            }
            status = nfu_failedToConverge; }
        else {
            ans = ( nm1 != 0 ) ? 1.0 / nm1 : -log(x) - EULER;   /* Set first term */
            fact = 1.0;
            for( i = 1; i <= MAXIT; i++ ) {
                fact *= -x / i;
                if( i != nm1 ) {
                    del = -fact / ( i - nm1 ); }
                else {
                    psi = -EULER;                   /* Compute psi(n) */
                    for( ii = 1; ii <= nm1; ii++ ) psi += 1.0 / ii;
                    del = fact * ( -log( x ) + psi );
                }
                ans += del;
                if( fabs( del ) < fabs( ans ) * EPS ) {
                    *value = ans;
                    return( nfu_Okay );
                }
            }
            status = nfu_failedToConverge;
        }
    }
    if( status != nfu_Okay ) smr_setReportError2p( smr, nfu_SMR_libraryID, status, "Failed to converge." );

    *value = ans;
    return( status );
}
