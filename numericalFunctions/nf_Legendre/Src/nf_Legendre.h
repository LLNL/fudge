/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
*/

#ifndef nf_Legendre_h_included
#define nf_Legendre_h_included

#include <nf_utilities.h>
#include <ptwXY.h>

#if defined __cplusplus
    extern "C" {
#endif

#define nf_Legendre_minMaxOrder 4
#define nf_Legendre_maxMaxOrder 128
#define nf_Legendre_sizeIncrement 8

typedef struct nf_Legendre_s nf_Legendre;

struct nf_Legendre_s {
    nfu_status status;
    int maxOrder;
    int allocated;          /* Will never be less than nf_Legendre_minMaxOrder. */
    double *Cls;
};

typedef nfu_status (*nf_Legendre_GaussianQuadrature_callback)( double x, double *y, void *argList );

/*
* Methods in nf_Legendre.c
*/
nf_Legendre *nf_Legendre_new( statusMessageReporting *smr, int initialSize, int maxOrder, double *Cls );
nfu_status nf_Legendre_initialize( statusMessageReporting *smr, nf_Legendre *nfL, int initialSize, int maxOrder );
nfu_status nf_Legendre_release( statusMessageReporting *smr, nf_Legendre *nfL );
nf_Legendre *nf_Legendre_free( nf_Legendre *nfL );
nf_Legendre *nf_Legendre_clone( statusMessageReporting *smr, nf_Legendre *nfL );
nfu_status nf_Legendre_reallocateCls( statusMessageReporting *smr, nf_Legendre *Legendre, int size, int forceSmallerResize );
nfu_status nf_Legendre_maxOrder( statusMessageReporting *smr, nf_Legendre *Legendre, int *maxOrder );
nfu_status nf_Legendre_allocated( statusMessageReporting *smr, nf_Legendre *Legendre, int *allocated );
nfu_status nf_Legendre_getCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double *Cl );
nfu_status nf_Legendre_setCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double Cl );
nfu_status nf_Legendre_normalize( statusMessageReporting *smr, nf_Legendre *Legendre );
nfu_status nf_Legendre_evauluateAtMu( statusMessageReporting *smr, nf_Legendre *nfL, double mu, double *P );
double nf_Legendre_PofL_atMu( int l, double mu );
ptwXYPoints *nf_Legendre_to_ptwXY( statusMessageReporting *smr, nf_Legendre *nfL, double accuracy, int biSectionMax, 
        int checkForRoots );
nf_Legendre *nf_Legendre_from_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY, int maxOrder );

/*
* Methods in nf_Legendre_GaussianQuadrature.c
*/
nfu_status nf_Legendre_GaussianQuadrature( int degree, double x1, double x2, nf_Legendre_GaussianQuadrature_callback func, void *argList, double *integral );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_Legendre_h_included. */
