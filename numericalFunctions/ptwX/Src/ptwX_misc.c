/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include "ptwX.h"

/*
************************************************************
*/
nfu_status ptwX_simpleWrite( statusMessageReporting *smr, ptwXPoints const *ptwX, FILE *f, char const *format ) {

    int64_t i1;
    double *p1 = ptwX->points;

    for( i1 = 0; i1 < ptwX->length; ++i1, ++p1 ) fprintf( f, format, *p1 );
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_simplePrint( statusMessageReporting *smr, ptwXPoints const *ptwX, char const *format ) {

    return( ptwX_simpleWrite( smr, ptwX, stdout, format ) );
}
