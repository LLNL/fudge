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

#include <nf_integration.h>
#include <nf_Legendre.h>

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

    return( Py_BuildValue( "di", integral, (int) evaluations ) );
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

    result = PyObject_CallFunction( parameters->func, "(d,O)", x, parameters->argList );
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

static char const doc[] = "A module that contains the quadrature integrators.";

#if PY_MAJOR_VERSION < 3

    PyMODINIT_FUNC initintegration( void ) {

        Py_InitModule3( "integration", nf_integration_CMiscPyMethods, doc );
    }

#else

    static struct PyModuleDef integration_CModule = {
        PyModuleDef_HEAD_INIT,
        "integration",         /* name of module */
        doc,                        /* module documentation, may be NULL */
        -1,                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        nf_integration_CMiscPyMethods
    };  
        
    PyMODINIT_FUNC PyInit_integration( void ) {
    
        return( PyModule_Create( &integration_CModule ) );
    }
        
#endif
