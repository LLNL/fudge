/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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

static PyObject *nf_exponentialIntegral_C( PyObject *self, PyObject *args );
static PyObject *nf_gamma_C( PyObject *self, PyObject *args );
static PyObject *nf_incompleteGamma_C( PyObject *self, PyObject *args, PyObject *keywords );
static PyObject *nf_erf_C( PyObject *self, PyObject *args, PyObject *keywords );
static PyObject *nf_specialFunctions_C_SetPyErrorExceptionReturnNull( const char *s, ... );
static void nf_specialFunctions_C_SetPyErrorExceptionFromSMR( PyObject *type, statusMessageReporting *smr );

/*
************************************************************
*/
static PyObject *nf_exponentialIntegral_C( PyObject *self, PyObject *args ) {

    int n;
    double En, x;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTuple( args, "id", &n, &x ) ) return( NULL );

    if( n < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid n: n = %d < 0", n ) );
    if( x < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e < 0", x ) );
    if( ( n < 2 ) && ( x == 0. ) ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid n for x = 0.: n = %d < 2", n ) );

    if( nf_exponentialIntegral( &smr, n, x, &En ) != nfu_Okay ) {
        nf_specialFunctions_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", En ) );
}
/*
************************************************************
*/
static PyObject *nf_gamma_C( PyObject *self, PyObject *args ) {

    double gamma, x;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTuple( args, "d", &x ) ) return( NULL );

    if( x <= 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e <= 0", x ) );

    if( nf_gammaFunction( &smr, x, &gamma ) != nfu_Okay ) {
        nf_specialFunctions_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", gamma ) );
}
/*
************************************************************
*/
static PyObject *nf_incompleteGamma_C( PyObject *self, PyObject *args, PyObject *keywords ) {

    int doComplementary = 0;
    double gamma, x, s;
    static char *kwlist[] = { "s", "x", "complementary", NULL };
    statusMessageReporting smr;
    nfu_status status_nf;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "dd|i", kwlist, &s, &x, &doComplementary ) ) return( NULL );

    if( x < 0 ) return( nf_specialFunctions_C_SetPyErrorExceptionReturnNull( "Invalid x: x = %e < 0", x ) );

    if( doComplementary ) {
        status_nf = nf_incompleteGammaFunctionComplementary( &smr, s, x, &gamma ); }
    else {
        status_nf = nf_incompleteGammaFunction( &smr, s, x, &gamma );
    }
    if( status_nf != nfu_Okay ) {
        nf_specialFunctions_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", gamma ) );
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
static void nf_specialFunctions_C_SetPyErrorExceptionFromSMR( PyObject *type, statusMessageReporting *smr ) {

    PyErr_SetString( type, smr_getMessage( smr_firstReport( smr ) ) );
    smr_release( smr );
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

static char const doc[] = "A module that contains special math functions not in the python math module.";

#if PY_MAJOR_VERSION < 3

    PyMODINIT_FUNC initspecialFunctions( void ) {

        Py_InitModule3( "specialFunctions", nf_specialFunctions_C_MiscPyMethods, doc );
    }

#else

    static struct PyModuleDef specialFunctions_CModule = {
        PyModuleDef_HEAD_INIT,
        "specialFunctions",         /* name of module */
        doc,                        /* module documentation, may be NULL */
        -1,                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        nf_specialFunctions_C_MiscPyMethods
    };

    PyMODINIT_FUNC PyInit_specialFunctions( void ) {

        return( PyModule_Create( &specialFunctions_CModule ) );
    }

#endif
