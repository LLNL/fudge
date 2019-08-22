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

#include <Python.h>
#include <pyport.h>
#include "structmember.h"
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>

#include <nf_specialFunctions.h>

#if( PY_VERSION_HEX < 0x02050000 )
#define Py_ssize_t int64_t 
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define ssizeobjargproc intobjargproc
#define ssizessizeobjargproc intintobjargproc
#endif

static PyObject *nf_exponentialIntegral_C( PyObject *self, PyObject *args );
static PyObject *nf_gamma_C( PyObject *self, PyObject *args );
static PyObject *nf_incompleteGamma_C( PyObject *self, PyObject *args, PyObject *keywords );
static PyObject *nf_amc_wigner_3j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_wigner_6j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_wigner_9j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_racah_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_clebsh_gordan_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_z_coefficient_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_zbar_coefficient_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_reduced_matrix_element_C( PyObject *self, PyObject *args );
static PyObject *nf_erf_C( PyObject *self, PyObject *args, PyObject *keywords );
static PyObject *nf_specialFunctions_C_SetPyErrorExceptionReturnNull( const char *s, ... );

DL_EXPORT( void ) initnf_specialFunctions_C( void );
/*
************************************************************
*/
static PyObject *nf_exponentialIntegral_C( PyObject *self, PyObject *args ) {

    int n;
    double En, x;
    nfu_status status_nf;

    if( !PyArg_ParseTuple( args, "id", &n, &x ) ) return( NULL );

    if( n < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid n: n = %d < 0", n ) );
    if( x < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e < 0", x ) );
    if( ( n < 2 ) && ( x == 0. ) ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid n for x = 0.: n = %d < 2", n ) );

    En = nf_exponentialIntegral( n, x, &status_nf );
    if( status_nf != nfu_Okay ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Error from nf_exponentialIntegral: %s", 
        nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", En ) );
}
/*
************************************************************
*/
static PyObject *nf_gamma_C( PyObject *self, PyObject *args ) {

    double gamma, x;
    nfu_status status_nf;

    if( !PyArg_ParseTuple( args, "d", &x ) ) return( NULL );

    if( x <= 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e <= 0", x ) );

    gamma = nf_gammaFunction( x, &status_nf );
    if( status_nf != nfu_Okay ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Error from nf_gammaFunction: %s", 
        nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", gamma ) );
}
/*
************************************************************
*/
static PyObject *nf_incompleteGamma_C( PyObject *self, PyObject *args, PyObject *keywords ) {

    int doComplementary = 0;
    double gamma, x, s;
    nfu_status status_nf;
    static char *kwlist[] = { "s", "x", "complementary", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "dd|i", kwlist, &s, &x, &doComplementary ) ) return( NULL );

    if( x < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e < 0", x ) );

    if( doComplementary ) {
        gamma = nf_incompleteGammaFunctionComplementary( s, x, &status_nf ); }
    else {
        gamma = nf_incompleteGammaFunction( s, x, &status_nf );
    }
    if( status_nf != nfu_Okay ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Error from nf_incompleteGammaFunction*: %s", 
        nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", gamma ) );
}



/*
************************************************************
*/
static PyObject *nf_amc_wigner_3j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_3j( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_wigner_6j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_6j( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_wigner_9j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6, j7, j8, j9;

    if( !PyArg_ParseTuple( args, "iiiiiiiii", &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_9j( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_racah_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_racah( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_clebsh_gordan_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, m1, m2;

    if( !PyArg_ParseTuple( args, "iiiii", &j1, &j2, &m1, &m2, &j3 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_clebsh_gordan( j1, j2, m1, m2, j3 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_zbar_coefficient_C( PyObject *self, PyObject *args ) {

    int j1, j2, l1, l2, ll, s;

    if( !PyArg_ParseTuple( args, "iiiiii", &l1, &j1, &l2, &j2, &s, &ll ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_zbar_coefficient( l1, j1, l2, j2, s, ll ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_z_coefficient_C( PyObject *self, PyObject *args ) {

    int j1, j2, l1, l2, ll, s;

    if( !PyArg_ParseTuple( args, "iiiiii", &l1, &j1, &l2, &j2, &s, &ll ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_z_coefficient( l1, j1, l2, j2, s, ll ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_reduced_matrix_element_C( PyObject *self, PyObject *args ) {

    int lt, st, jt, l0, j0, l1, j1;

    if( !PyArg_ParseTuple( args, "iiiiiii", &lt, &st, &jt, &l0, &j0, &l1, &j1 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_reduced_matrix_element( lt, st, jt, l0, j0, l1, j1 ) ) );
}




/*
************************************************************
*/
static PyObject *nf_erf_C( PyObject *self, PyObject *args, PyObject *keywords ) {

    int doComplementary = 0;
    double erf_, x;
    static char *kwlist[] = { "x", "complementary", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "d|i", kwlist, &x, &doComplementary ) ) return( NULL );

    if( doComplementary ) {
        erf_ = erfc( x ); }
    else {
        erf_ = erf( x );
    }
    return( Py_BuildValue( "d", erf_ ) );
}
/*
************************************************************
*/
static PyObject *nf_specialFunctions_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

    va_list args;
    char Str[1024];

    va_start( args, s );
    vsnprintf( Str, sizeof( Str ), s, args );
    Str[sizeof( Str ) - 1] = 0;
    PyErr_SetString( PyExc_Exception, Str );
    va_end( args );
    return( NULL );
}
/*
************************************************************
*/
static PyMethodDef nf_specialFunctions_C_MiscPyMethods[] = {

    { "exponentialIntegral", (PyCFunction) nf_exponentialIntegral_C, METH_VARARGS, 
        "exponentialIntegral( n, x )\n\n" \
        "Returns the value of the exponential integral defined as\n\n" \
        "     - infinity    \n" \
        "    |              \n" \
        "    |     -x t     \n" \
        "    |    e         \n" \
        "    |  ------ dt   \n" \
        "    |     n        \n" \
        "    |    t         \n" \
        "    |              \n" \
        "   _  1            \n" \
        "\n" \
        "\nArguments are:\n" \
        "    n  an integer that must be greater than or equal to 0\n" \
        "    x  a float\n" },
    { "gamma", (PyCFunction) nf_gamma_C, METH_VARARGS,
        "gamma( s )\n\n" \
        "Returns the value of the gamma function defined as\n\n" \
        "      -| infinity   \n" \
        "     |              \n" \
        "     |   s-1  -t    \n" \
        "     |  t    e   dt \n" \
        "     |              \n" \
        "     |              \n" \
        "   |_  0            \n" \
        "\n" \
        "\nArguments are:\n" \
        "    s  a float that must be greater than 0\n" },
    { "incompleteGamma", (PyCFunction) nf_incompleteGamma_C, METH_VARARGS | METH_KEYWORDS,
        "incompleteGamma( s, x, complementary = False )\n\n" \
        "Returns the value of the incomplete gamma function defined as\n\n" \
        "      -| x          \n" \
        "     |              \n" \
        "     |   s-1  -t    \n" \
        "     |  t    e   dt \n" \
        "     |              \n" \
        "     |              \n" \
        "   |_  0            \n" \
        "\n" \
        "If the third arugment (i.e., complementary) is False (default), the integral is as shown from 0\n" \
        "to x, otherwise, the integral is from x to infinite.\n" \
        "\nArguments are:  ([o] implies optional argument)\n" \
        "    s                  a float that must be greater than 0,\n" \
        "    x                  a float that must be greater than 0,\n" \
        "    complementary  [o] an integer [default False].\n" },
    { "erf", (PyCFunction) nf_erf_C, METH_VARARGS | METH_KEYWORDS,
        "erf( x, complementary = False )\n\n" \
        "Returns the value of the error function defined as\n\n" \
        "          -| x       \n" \
        "         |           \n" \
        "    2    |     2     \n" \
        "  -----  |   -t      \n" \
        "    __   |  e   dt   \n" \
        "  \\/pi   |           \n" \
        "         |           \n" \
        "       |_  0         \n" \
        "\n" \
        "If the second arugment (i.e., complementary) is False (default), the integral is as shown from 0\n" \
        "to x, otherwise, the integral is from x to infinite.\n" \
        "\nArguments are:  ([o] implies optional argument)\n" \
        "    x                  a float,\n" \
        "    complementary  [o] an integer [default False].\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
static PyMethodDef nf_specialFunctions_C_AngularMomentumCouplingMethods[] = {
/*
Testing issues:
    wigner_3j( 6, 6, 4, 0, 0, 0 ) / wigner_3j( 4, 6, 6, 0, 0, 0 ) - 1 = -4.392597396929432e-12
*/
    { "wigner_3j", (PyCFunction) nf_amc_wigner_3j_C, METH_VARARGS,
        "wigner_3j( j1, j2, j3, j4, j5, j6 )\n\n" \
        "Wigner's 3J symbol (similar to Clebsh-Gordan)\n\n" \
        "    = / j1 j2 j3 \\\n" \
        "      \\ j4 j5 j6 /\n\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "wigner_6j", (PyCFunction) nf_amc_wigner_6j_C, METH_VARARGS,
        "wigner_6j( j1, j2, j3, j4, j5, j6 )\n\n" \
        "Wigner's 6J symbol (similar to Racah)\n\n" \
        "        = { j1 j2 j3 }\n" \
        "          { j4 j5 j6 }\n\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "wigner_9j", (PyCFunction) nf_amc_wigner_9j_C, METH_VARARGS,
        "wigner_9j( j1, j2, j3, j4, j5, j6, j7, j8, j9 )\n\n" \
        "Wigner's 9J symbol\n\n" \
        "      / j1 j2 j3 \\\n" \
        "    = | j4 j5 j6 |\n" \
        "      \\ j7 j8 j9 /\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "racah", (PyCFunction) nf_amc_racah_C, METH_VARARGS,
        "racah( j1, j2, l2, l1, j3, l3 )\n\n" \
        "Racah coefficient\n\n" \
        "    = W(j1, j2, l2, l1 ; j3, l3)\n\n" \
        "    = (-1)^(j1+j2+l1+l2) * { j1 j2 j3 }\n" \
        "                           { l1 l2 l3 }\n\n" \
        "\nArguments are: \n" \
        "    j1, j2, l2, l1, j3, l3     all integers, 2x angular momenta\n\n"\
        "We use the Racah coefficient definition in Edmonds (A.R. Edmonds, 'Angular Momentum in Quantum Mechanics', Princeton (1980).\n"\
        "Note that the call signature of W(...) appears jumbled, but hey, that's the convention.\n"\
        "This convention is exactly that used by Blatt-Biedenharn (Rev. Mod. Phys. 24, 258 (1952)) too."
    },
    { "clebsh_gordan", (PyCFunction) nf_amc_clebsh_gordan_C, METH_VARARGS,
        "clebsh_gordan( j1, j2, j3, m1, m2)\n\n" \
        "Clebsh-Gordan coefficient\n\n" \
        "    = <j1,j2,m1,m2|j3,m1+m2>\n\n" \
        "    = (-)^(j1-j2+m1+m2) * sqrt(2*j3+1) * / j1 j2   j3   \\\n" \
        "                                         \\ m1 m2 -m1-m2 /\n\n" \
        "Note: Last value m3 is preset to m1+m2.  Any other value will evaluate to 0.0.\n\n"\
        "\nArguments are: \n" \
        "    j1, j2, j3     integers, 2x the angular momentum (so for j=1/2 use j=1)\n" \
        "    m1, m2, m3     integers, 2x the projection of the angular momentum onto the z axis\n" },
    { "z_coefficient", (PyCFunction) nf_amc_z_coefficient_C, METH_VARARGS,
        "z_coefficient( j1, j2, l1, l2, ll, s )\n\n" \
        "Biedenharn's Z-coefficient coefficient\n\n" \
        "    =  Z(l1  j1  l2  j2 | S L )\n\n" \
        "\nArguments are: \n" \
        "    j1, j2     integers, 2x the total angular momentum (so for j=1/2 use j=1)\n" \
        "    l1, l2     integers, 2x the orbital angular momentum\n" \
        "    ll         integer, 2x the orbital angular momentum that l1 and l2 couple up to,\n" \
        "    s          integer, 2x the projection of the ll onto the z axis\n" },
    { "zbar_coefficient", (PyCFunction) nf_amc_zbar_coefficient_C, METH_VARARGS,
        "zbar_coefficient( j1, j2, l1, l2, ll, s )\n\n" \
        "Lane & Thomas's Zbar-coefficient coefficient\n\n" \
        "    =  ZBar(l1  j1  l2  j2 | S L )\n\n" \
        "    = (-i)^( -l1 + l2 + ll ) * Z(l1  j1  l2  j2 | S L )\n\n"
        "From Lane & Thomas Rev. Mod. Phys. 30, 257-353 (1958).\n"
        "Note, Lane & Thomas define this because they did not like the different phase convention in Blatt & Biedenharn's Z coefficient.  They changed it to get better time-reversal behavior.\n"
        "Froehner uses Lane & Thomas convention as does T. Kawano.\n\n"
        "\nArguments are: \n" \
        "    j1, j2     integers, 2x the total angular momentum (so for j=1/2 use j=1)\n" \
        "    l1, l2     integers, 2x the orbital angular momentum\n" \
        "    ll         integer, 2x the orbital angular momentum that l1 and l2 couple up to,\n" \
        "    s          integer, 2x the projection of the ll onto the z axis\n" },
    { "reduced_matrix_element", (PyCFunction) nf_amc_reduced_matrix_element_C, METH_VARARGS,
        "reduced_matrix_element( lt, st, jt, l0, j0, l1, j1 )\n\n" \
        "Reduced Matrix Element for Tensor Operator\n\n" \
        "    = < l1j1 || T(YL,sigma_S)J || l0j0 >\n\n" \
        "From M.B.Johnson, L.W.Owen, G.R.Satchler, Phys. Rev. 142, 748 (1966)\n\n" \
        "Note: definition differs from JOS by the factor sqrt(2j1+1)\n\n" \
        "\nArguments are: \n" \
        "    lt     integer, 2x the orbital angular momentum of the tensor operator,\n" \
        "    st     integer, 2x the projection of lt onto the z axis,\n" \
        "    jt     integer, 2x the total angular momentum of the tensor operator,\n" \
        "    l0     integer, 2x the orbital angular momentum of the ket,\n" \
        "    j0     integer, 2x the total angular momentum of the ket,\n" \
        "    l1     integer, 2x the orbital angular momentum of the bra,\n" \
        "    j1     integer, 2x the total angular momentum of the bra.\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
DL_EXPORT( void ) initspecialFunctions( void ) {

    PyObject *m;

    if( ( m = Py_InitModule3( "specialFunctions", nf_specialFunctions_C_MiscPyMethods, 
        "A module that contains special math functions not in the python math module." ) ) == NULL ) return;
    if( ( m = Py_InitModule3( "angularMomentumCoupling", nf_specialFunctions_C_AngularMomentumCouplingMethods, 
        "A module that contains special math functions for angular momentum coupling in physics." ) ) == NULL ) return;
}
