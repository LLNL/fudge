/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef specialFunctions_h_included
#define specialFunctions_h_included

#include <math.h>
#include <float.h>
#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

double nf_polevl( double x, double coef[], int N );
double nf_p1evl( double x, double coef[], int N );
double nf_exponentialIntegral( int n, double x, nfu_status *status );
double nf_gammaFunction( double x, nfu_status *status );
double nf_logGammaFunction( double x, nfu_status *status );
double nf_incompleteGammaFunction( double a, double x, nfu_status *status );
double nf_incompleteGammaFunctionComplementary( double a, double x, nfu_status *status );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwXY_h_included. */
