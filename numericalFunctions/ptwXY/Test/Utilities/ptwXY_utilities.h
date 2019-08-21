/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef nf_utility_h_included
#define nf_utility_h_included

#include <nf_utilities.h>
#include <ptwXY.h>

#if defined __cplusplus
    extern "C" {
#endif

int nfu_ptwXY_cmp( ptwXYPoints *p1, ptwXYPoints *p2, int verbose, double frac );

#if defined __cplusplus
    }
#endif

#endif              /* End of nf_utility_h_included. */
