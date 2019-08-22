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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>
*/

#include <Python.h>
#include <pyport.h>
#include "structmember.h"
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include <pointwiseXY_C.h>
#include <nf_Legendre.h>

#if( PY_VERSION_HEX < 0x02050000 )
#define Py_ssize_t int64_t 
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define ssizeobjargproc intobjargproc
#define ssizessizeobjargproc intintobjargproc
#endif

staticforward PyTypeObject nf_Legendre_CPyType;

typedef struct nf_Legendre_CPy_s {
    PyObject_HEAD
    nf_Legendre *nfL;
} nf_Legendre_CPy;

#define is_nf_Legendre_CPyObject( v ) ((v)->ob_type == &nf_Legendre_CPyType)

static char nf_Legendre_C__doc__[] = 
    "The Legendre class stores the coefficients for a Legendre series.\n" \
    "\n" \
    "Constructor arguments are:\n" \
    "   Cls             the list of Legendre coefficients,\n" \
    "   initialSize     the initial size of the memory allocated for storing the Legendre coefficients (default = 0),\n";

static nf_Legendre_CPy *nf_Legendre_CNewInitialize( void );
static int nf_Legendre_C__init__( nf_Legendre_CPy *self, PyObject *args, PyObject *keywords );
static void nf_Legendre_C_dealloc( PyObject *self );
static PyObject *nf_Legendre_C__repr__( nf_Legendre_CPy *self );
static Py_ssize_t nf_Legendre_C__len__( nf_Legendre_CPy *self );
static PyObject *nf_Legendre_C__getitem__( nf_Legendre_CPy *self, Py_ssize_t index );
static int nf_Legendre_C__setitem__( nf_Legendre_CPy *self, Py_ssize_t index, PyObject *value );
static PyObject *nf_Legendre_C_getMaxOrder( nf_Legendre_CPy *self );
static PyObject *nf_Legendre_C_normalize( nf_Legendre_CPy *self );
static PyObject *nf_Legendre_C_toString( nf_Legendre_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *nf_Legendre_C_toString2( nf_Legendre_CPy *self, char *format, char *sep );
static char *nf_Legendre_C_toString_isFormatForDouble( char *format, int returnFormat );
static int nf_Legendre_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds );
static int nf_Legendre_C_PyNumberToFloat( PyObject *n, double *d );
static PyObject *nf_Legendre_C_SetPyErrorExceptionReturnNull( const char *s, ... );
static int nf_Legendre_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... );

static PyObject *nf_Legendre_C_from_pointwiseXY_C( PyObject *self, PyObject *args );
static PyObject *nf_Legendre_C_getMaxMaxOrder( PyObject *self );

DL_EXPORT( void ) initLegendre( void );
/*
******************** nf_Legendre_CNewInitialize ************************
*/
static nf_Legendre_CPy *nf_Legendre_CNewInitialize( void ) {

    nf_Legendre_CPy *self = (nf_Legendre_CPy *) PyObject_New( nf_Legendre_CPy, &nf_Legendre_CPyType );

    if( self ) {
        self->nfL = NULL;
    }
    return( self );
}
/*
************************ __init__ ************************************
*/
static int nf_Legendre_C__init__( nf_Legendre_CPy *self, PyObject *args, PyObject *keywords ) {

    int initialSize = 0, maxOrder = -1;
    double *Cls = NULL;
    nfu_status status_nf;
    static char *kwlist[] = { "Cls", "initialSize", NULL };
    PyObject *ClsPy = NULL;

    self->nfL = NULL;
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|Oi", kwlist, &ClsPy, &initialSize ) ) return( -1 );

    if( ClsPy != NULL ) {   /* Check that ClsPy is list of doubles and put them into Cls. ?????????? */
        if( ( maxOrder = nf_Legendre_C_pythonDoubleListToCList( ClsPy, &Cls ) ) < 0 ) return( -1 );
        maxOrder--;
    }
    if( ( self->nfL = nf_Legendre_new( 0, maxOrder, Cls, &status_nf ) ) == NULL ) {
        if( Cls != NULL ) free( Cls );
        PyErr_NoMemory( );
        return( -1 );
    }
    if( Cls != NULL ) free( Cls );
    return( 0 );
}
/*
************************************************************
*/
static void nf_Legendre_C_dealloc( PyObject *self ) {

    nf_Legendre_CPy *nfL = (nf_Legendre_CPy *) self;

    if( nfL->nfL != NULL ) {
        nf_Legendre_free( nfL->nfL );
    }
    self->ob_type->tp_free( self );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C__repr__( nf_Legendre_CPy *self ) {

    return( nf_Legendre_C_toString2( self, "%16.8e", " " ) );
}
/*
************************************************************
*/
static Py_ssize_t nf_Legendre_C__len__( nf_Legendre_CPy *self ) {

    return( (Py_ssize_t) nf_Legendre_maxOrder( self->nfL ) + 1 );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C__getitem__( nf_Legendre_CPy *self, Py_ssize_t l ) {

    nfu_status status_nf;
    double d;

    if( l < 0 ) l += nf_Legendre_maxOrder( self->nfL ) + 1;
    d = nf_Legendre_getCl( self->nfL, (int) l, &status_nf );

    if( status_nf != nfu_Okay ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( NULL );
    }

    return( (PyObject *) Py_BuildValue( "d", d ) );
}
/*
************************************************************
*/
static int nf_Legendre_C__setitem__( nf_Legendre_CPy *self, Py_ssize_t l_, PyObject *value ) {
/*
*   This routine can be called either when "self[l_] = number" or "del self[l_]".
*/
    int l = (int) l_;
    double Cl;

    if( l < 0 ) l += nf_Legendre_maxOrder( self->nfL ) + 1;

    if( value == NULL ) {                   /* Currently, deleting a coefficient is not supported. */
        return( 0 ); }
    else {
        if( nf_Legendre_C_PyNumberToFloat( value, &Cl ) < 0 ) return( -1 );
        if( nf_Legendre_setCl( self->nfL, l, Cl ) != nfu_Okay )
            return( nf_Legendre_C_SetPyErrorExceptionReturnMinusOne( "l = %d out of range", l_ ) );
    }
    return( 0 );
}
/*
************************************************************
*/
static PySequenceMethods nf_Legendre_CPy_sequence = {
    (lenfunc) nf_Legendre_C__len__,                 /* sq_length */
    0,                                              /* sq_concat */
    0,                                              /* sq_repeat */
    (ssizeargfunc) nf_Legendre_C__getitem__,        /* sq_item */
    (ssizessizeargfunc) 0,                          /* sq_slice */
    (ssizeobjargproc) nf_Legendre_C__setitem__,     /* sq_ass_item */
    (ssizessizeobjargproc) 0,                       /* sq_ass_slice */
    0,                                              /* sq_contains */
    0,                                              /* sq_inplace_concat */
    0                                               /* sq_inplace_repeat */
};
/*
************************************************************
*/
static PyObject *nf_Legendre_C_getMaxOrder( nf_Legendre_CPy *self ) {

    return( Py_BuildValue( "i", nf_Legendre_maxOrder( self->nfL ) ) );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_normalize( nf_Legendre_CPy *self ) {

    nf_Legendre *nfL;
    nfu_status status_nf;
    nf_Legendre_CPy *nfL_Py;

    if( ( nfL = nf_Legendre_clone( (nf_Legendre *) self->nfL, &status_nf ) ) == NULL )
        return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "Error from nf_Legendre_clone: %s", nfu_statusMessage( status_nf ) ) );
    if( ( status_nf = nf_Legendre_normalize( nfL ) ) != nfu_Okay ) {
        nf_Legendre_free( nfL );
        return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "Error from nf_Legendre_clone: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nfL_Py = nf_Legendre_CNewInitialize( ) ) == NULL ) return( NULL );
    nfL_Py->nfL = nfL;
    return( (PyObject *) nfL_Py );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_toPointwiseLinear( nf_Legendre_CPy *self, PyObject *args, PyObject *keywords ) {

    int biSectionMax = 16, checkForRoots = 0, infill = 1, safeDivide = 1;
    double accuracy;
    static char *kwlist[] = { "accuracy", "biSectionMax", "checkForRoots", "infill", "safeDivide", NULL };
    nfu_status status_nf;
    ptwXYPoints *ptwXY;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "d|iiii", kwlist, &accuracy, &biSectionMax, &checkForRoots, &infill, &safeDivide ) ) return( NULL );

    if( ( ptwXY = nf_Legendre_to_ptwXY( self->nfL, accuracy, biSectionMax, checkForRoots, &status_nf ) ) == NULL ) 
        return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "Error from nf_Legendre_to_ptwXY: %s", nfu_statusMessage( status_nf ) ) );

    return( pointwiseXY_C_factory_create( ptwXY, infill, safeDivide ) );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_toString( nf_Legendre_CPy *self, PyObject *args, PyObject *keywords ) {

    char *format = "%16.8e", *sep = " ";
    static char *kwlist[] = { "format", "sep", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ss", kwlist, &format, &sep ) ) return( NULL );
    return( nf_Legendre_C_toString2( self, format, sep ) );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_toString2( nf_Legendre_CPy *self, char *format, char *sep ) {

    int l, extraSpace = 3, dummyLen;     /* three more for reserve. */
    int length, strLength, sepLength = (int) strlen( sep );
    double d;
    nf_Legendre *nfL = self->nfL;
    char *s, *e, *p, dummy[1024];
    PyObject *str;
    nfu_status status_nf;

    if( ( length = nf_Legendre_maxOrder( nfL ) ) < 0 ) return( PyString_FromString( "" ) );
    length++;

    s = nf_Legendre_C_toString_isFormatForDouble( format, 0 );
    if( ( s = nf_Legendre_C_toString_isFormatForDouble( s, 1 ) ) == NULL ) 
        return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "invalid format = '%s' for converting double to a string", format ) );

    dummyLen = snprintf( dummy, 0, format, -3.1238972718972345e-21 );
    if( dummyLen < 0 ) return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "Invalid format string = '%s'", format ) );
    strLength = length * ( dummyLen + sepLength + extraSpace );

    if( ( s = (char *) malloc( strLength * sizeof( char ) ) ) == NULL ) {
        PyErr_NoMemory( );
        return( NULL );
    }
    s[0] = 0;
    for( l = 0, e = s; l < length; l++ ) {
        d = nf_Legendre_getCl( nfL, l, &status_nf );
        dummyLen = sprintf( e, format, d );
        e += dummyLen;
        if( l < length - 1 ) {
            for( p = sep; *p; p++, e++ ) *e = *p;
            *e = 0;
        }
    }
    str = PyString_FromString( s );
    free( s );
    return( str );
}
/*
************************************************************
*/
static char *nf_Legendre_C_toString_isFormatForDouble( char *format, int returnFormat ) {

    size_t i;
    char *s;

    if( format == NULL ) return( format );
    for( s = format; *s; s++ ) {
        if( *s == '%' ) {
            if( s[1] != '%' ) break;
            s++;
        }
    }
    if( returnFormat ) return( ( *s == 0 ) ? format : NULL );
    if( *s == 0 ) return( NULL );
    s++;
    i = strspn( s, "+-0123456789." );
    s = &(s[i]);
    if( ( *s != 'e' ) && ( *s != 'E' ) && ( *s != 'f' ) && ( *s != 'F' ) && ( *s != 'g' ) && ( *s != 'G' ) ) return( NULL );
    s++;
    return( s );
}
/*
************************************************************
*/
static int nf_Legendre_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds ) {

    int status = 0, length, i;
    double *d;
    PyObject *item, *iterator;

    *ds = NULL;
    if( ( iterator = PyObject_GetIter( PyDoubleList ) ) == NULL ) return( -1 );
    if( ( length = (int) PySequence_Size( PyDoubleList ) ) != (int) 0 ) {
        if( ( *ds = (double *) malloc( length * sizeof( double ) ) ) == NULL ) {
            PyErr_NoMemory( );
            Py_DECREF( iterator );
            return( -1 );
        }
        for( i = 0, d = *ds, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i++, d++ ) {
            if( ( status = nf_Legendre_C_PyNumberToFloat( item, d ) ) != 0 )
                nf_Legendre_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i );
            Py_DECREF( item );
        }
    }
    Py_DECREF( iterator );
    if( status != 0 ) {
        length = -1;
        free( *ds );
        *ds = NULL;
    }
    return( length );
}
/*
************************************************************
*/
static int nf_Legendre_C_PyNumberToFloat( PyObject *n, double *d ) {

    *d = PyFloat_AsDouble( n );
    if( PyErr_Occurred( ) != NULL ) return( -1 );
    return( 0 );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

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
static int nf_Legendre_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... ) {

    va_list args;
    char Str[1024];

    va_start( args, s );
    vsnprintf( Str, sizeof( Str ), s, args );
    Str[sizeof( Str ) - 1] = 0;
    PyErr_SetString( PyExc_Exception, Str );
    va_end( args );
    return( -1 );
}
/*
************************************************************
*/
static PyMethodDef nf_Legendre_CPyMethods[] = {

    { "maxOrder", (PyCFunction) nf_Legendre_C_getMaxOrder, METH_NOARGS, 
        "Returns the maxOrder of self or -1 if no coefficients are defined.\n" \
        "\nArguments are: (this method does not take any arguments).\n" },
    { "normalize", (PyCFunction) nf_Legendre_C_normalize, METH_NOARGS, 
        "Returns a new Legendre instance that is the a clone of self, except that it is normalized to 1.\n" \
        "\nArguments are: (this method does not take any arguments).\n" },
    { "toPointwiseLinear", (PyCFunction) nf_Legendre_C_toPointwiseLinear, METH_VARARGS | METH_KEYWORDS, 
        "Returns a pointwiseXY_C instance of self by evaluating self at enough mu values to represent the Legendre series\n" \
        "as the function f(mu) for -1 <= mu <= 1 to accuracy. Bisection is used to fill in the domain to the desired accuracy.\n" \
        "The number of bisections is controlled by the biSectionMax argument. That is, the bisecting is stopped when either\n" \
        "accuracy or biSectionMax is reached. However, if checkForRoots is True and f(mu) is determined to have a zero between\n" \
        "two mu-values, then an additional point is added at the zero of f(mu).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   accuracy        the desired accuracy of the pointwise representation of the Legendre series,\n" \
        "   biSectionMax    [o] the maximum number of bisections (default = 16),\n" \
        "   checkForRoots   [o] if True, check for roots (default = False),\n" \
        "   infill          [o] see module pointwiseXY_C for meaning (default = True),\n" \
        "   safeDivide      [o] see module pointwiseXY_C for meaning (default = True).\n" },
    { "toString", (PyCFunction) nf_Legendre_C_toString, METH_VARARGS | METH_KEYWORDS, 
        "Returns a string representation of the Legendre Series (i.e., coefficients).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   format      [o] the format specifier for a coefficients (default = '%16.8e'),\n" \
        "   sep         [o] a string to use as the separator between coefficients (default = ' ').\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
static PyTypeObject nf_Legendre_CPyType = {
    PyObject_HEAD_INIT( NULL )
    0,                                          /* ob_size        */
    "Legendre.Series",                          /* tp_name        */
    sizeof( nf_Legendre_CPy ),                  /* tp_basicsize   */
    0,                                          /* tp_itemsize    */
    /* methods */ 
    (destructor) nf_Legendre_C_dealloc,         /* tp_dealloc     */
    0,                                          /* tp_print       */
    0,                                          /* tp_getattr     */
    0,                                          /* tp_setattr     */
    0,                                          /* tp_compare     */
    (reprfunc) nf_Legendre_C__repr__,           /* tp_repr        */
    0,                                          /* tp_as_number   */
    & nf_Legendre_CPy_sequence,                 /* tp_as_sequence */
    0,                                          /* tp_as_mapping  */
    0,                                          /* tp_hash        */
    0,                                          /* tp_call        */
    0,                                          /* tp_str         */
    0,                                          /* tp_getattro    */
    0,                                          /* tp_setattro    */
    0,                                          /* tp_as_buffer   */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, /* tp_flags     */
    nf_Legendre_C__doc__,                       /* tp_doc         */
    0,                                          /* tp_traverse    */
    0,                                          /* tp_clear       */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter        */
    0,                                          /* tp_iternext    */
    nf_Legendre_CPyMethods,                     /* tp_methods     */
    0,                                          /* tp_members     */
    0,                                          /* tp_getset      */
    0,                                          /* tp_base        */
    0,                                          /* tp_dict        */
    0,                                          /* tp_descr_get   */
    0,                                          /* tp_descr_set   */
    0,                                          /* tp_dictoffset  */
    (initproc) nf_Legendre_C__init__,           /* tp_init        */
    0,                                          /* tp_alloc       */
    0                                           /* tp_new         */
};
/*
************************************************************
*/
static PyObject *nf_Legendre_C_from_pointwiseXY_C( PyObject *self, PyObject *args ) {

    int maxOrder;
    nfu_status status_nf;
    ptwXYPoints *ptwXY;
    PyObject *ptwXY_Py;
    nf_Legendre *nfL;
    nf_Legendre_CPy *nfL_Py;

    if( !PyArg_ParseTuple( args, "Oi", &ptwXY_Py, &maxOrder ) ) return( NULL );
    if( ( ptwXY = pointwiseXY_C_factory_get_ptwXYPoints( ptwXY_Py ) ) == NULL ) return( NULL );

    if( ( nfL = nf_Legendre_from_ptwXY( ptwXY, maxOrder, &status_nf ) ) == NULL )
        return( nf_Legendre_C_SetPyErrorExceptionReturnNull( "Error from nf_Legendre_from_ptwXY: %s", nfu_statusMessage( status_nf ) ) );

    if( ( nfL_Py = nf_Legendre_CNewInitialize( ) ) == NULL ) {
        nf_Legendre_free( nfL ); }
    else {
        nfL_Py->nfL = nfL;
    }
    return( (PyObject *) nfL_Py );
}
/*
************************************************************
*/
static PyObject *nf_Legendre_C_getMaxMaxOrder( PyObject *self ) {

    return( Py_BuildValue( "i", nf_Legendre_maxMaxOrder ) );
}
/*
************************************************************
*/
static PyMethodDef nf_Legendre_CMiscPyMethods[] = {

    { "from_pointwiseXY_C", (PyCFunction) nf_Legendre_C_from_pointwiseXY_C, METH_VARARGS, 
        "This function returns a Legendre class representation of a pointwiseXY_C instance.\n" \
        "\nArguments are:\n" \
        "   ptwXYs  A pointwiseXY_C instance.\n" \
        "   order   The maximum order for the Legendre series representation of ptwXYs.\n" },
    { "maxMaxOrder", (PyCFunction) nf_Legendre_C_getMaxMaxOrder, METH_NOARGS, 
        "The Legendre class limits the maximum Legendre order that an instance can have.\n" \
        "This function returns the largest allowed Legendre order.\n" \
        "\nArguments are: (this method does not take any arguments).\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
DL_EXPORT( void ) initLegendre( void ) {

    PyObject *m;

    nf_Legendre_CPyType.tp_new = PyType_GenericNew;
    if( PyType_Ready( &nf_Legendre_CPyType ) < 0 ) return;

    if( ( m = Py_InitModule3( "Legendre", nf_Legendre_CMiscPyMethods, "A module that contains the Legendre class." ) ) == NULL ) return;

    if( import_pointwiseXY_C( ) < 0 ) return;

    Py_INCREF( &nf_Legendre_CPyType );
    PyModule_AddObject( m, "Series", (PyObject *) &nf_Legendre_CPyType );
}
