/*
# <<BEGIN-copyright>>
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
static PyMethodDef nf_specialFunctions_CMiscPyMethods[] = {

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
DL_EXPORT( void ) initnf_specialFunctions_C( void ) {

    PyObject *m;

    if( ( m = Py_InitModule3( "nf_specialFunctions_C", nf_specialFunctions_CMiscPyMethods, 
        "A module that contains special math functions not in the python math module." ) ) == NULL ) return;
}
