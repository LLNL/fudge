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

#include <nf_integration.h>
#include <nf_Legendre.h>

#if( PY_VERSION_HEX < 0x02050000 )
#define Py_ssize_t int64_t 
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define ssizeobjargproc intobjargproc
#define ssizessizeobjargproc intintobjargproc
#endif

typedef struct GnG_parameters_s {
    int degree;
    PyObject *func;
    PyObject *argList;
} GnG_parameters;

static PyObject *nf_GnG_adaptiveQuadrature_C( PyObject *self, PyObject *args );
static nfu_status nf_GnG_adaptiveQuadrature_C_callback( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, double x1,
        double x2, double *integral );
static nfu_status nf_GnG_adaptiveQuadrature_C_integrandFunction( double x, double *y, void *argList );
static PyObject *nf_integration_C_SetPyErrorExceptionReturnNull( const char *s, ... );

DL_EXPORT( void ) initintegration( void );
/*
************************************************************
*/
static PyObject *nf_GnG_adaptiveQuadrature_C( PyObject *self, PyObject *args ) {

    nfu_status status_nf;
    int maxDepth = nf_GnG_adaptiveQuadrature_MaxMaxDepth, degree;
    long int evaluations = 0;
    double x1, x2, tolerance, integral;
    GnG_parameters argList;
    PyObject *integrandFunction_Py, *argList_Py;

    if( !PyArg_ParseTuple( args, "iOOddd|i", &degree, &integrandFunction_Py, &argList_Py, &x1, &x2, &tolerance, &maxDepth ) ) return( NULL );
    argList.degree = degree;
    argList.func = integrandFunction_Py;
    argList.argList  = argList_Py;

    status_nf = nf_GnG_adaptiveQuadrature( nf_GnG_adaptiveQuadrature_C_callback, nf_GnG_adaptiveQuadrature_C_integrandFunction, 
        (void *) &argList, x1, x2, maxDepth, tolerance, &integral, &evaluations );
    if( status_nf != nfu_Okay ) return( nf_integration_C_SetPyErrorExceptionReturnNull( "Error from nf_GnG_adaptiveQuadrature: %s", nfu_statusMessage( status_nf ) ) );

    return( Py_BuildValue( "(d,i)", integral, (int) evaluations ) );
}
/*
************************************************************
*/
static nfu_status nf_GnG_adaptiveQuadrature_C_callback( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, 
        double x1, double x2, double *integral ) {

    GnG_parameters *parameters = argList;

    return( nf_Legendre_GaussianQuadrature( parameters->degree, x1, x2, integrandFunction, argList, integral ) );
}
/*
************************************************************
*/
static nfu_status nf_GnG_adaptiveQuadrature_C_integrandFunction( double x, double *y, void *argList ) {

    nfu_status status_nf = nfu_Okay;
    GnG_parameters *parameters = (GnG_parameters *) argList;
    PyObject *result;

    result = PyEval_CallFunction( parameters->func, "(d,O)", x, parameters->argList );
    if( result == NULL ) return( nfu_badInput );

    if( PyFloat_Check( result ) ) {
        *y = PyFloat_AsDouble( result ); }
    else {
        status_nf = nfu_badInput;
    }

    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *nf_integration_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

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
static PyMethodDef nf_integration_CMiscPyMethods[] = {

    { "adaptiveQuadrature_GnG", (PyCFunction) nf_GnG_adaptiveQuadrature_C, METH_VARARGS, 
        "This function uses an adaptive quadrature function to integrate an integrand from x1 to x2 to tolerance (i.e., integral_x1^x2 dx y(x))." \
        " The adaptive quadrature continuously bisects the domain until each region is small enough so the either" \
        " tolerance is met or maxDepth bisects have occured. The user must supply an integrand function that takes" \
        " x and argList as its arguments and returns the value of the integrand at xi (i.e., y(x)).\n" \
        "\nArguments are: [o] are optional arguments,\n" \
        "   degree                  the degree of the Legendre quadrature function to use (typically in the range 2 to 6).\n" \
        "   integrand               a python function that calculates the value of the integrand give x and argList.\n" \
        "   argList                 a python object passed to the integrand function.\n" \
        "   x1                      the lower limit of integration.\n" \
        "   x2                      the upper limit of integration.\n" \
        "   tolerance               the tolerance for the integration.\n" \
        "   maxDepth            [o] maximum recursive depth.\n" \
        "   .\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
DL_EXPORT( void ) initintegration( void ) {

    PyObject *m;

    if( ( m = Py_InitModule3( "integration", nf_integration_CMiscPyMethods, "A module that contains the quadrature integrators." ) ) == NULL ) return;
}
