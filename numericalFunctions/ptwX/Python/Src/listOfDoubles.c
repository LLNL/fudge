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

#include <ptwX.h>

#if( PY_VERSION_HEX < 0x02050000 )
#define Py_ssize_t int64_t 
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define ssizeobjargproc intobjargproc
#define ssizessizeobjargproc intintobjargproc
#endif

staticforward PyTypeObject listOfDoubles_CPyType;

typedef struct listOfDoubles_CPy_s {
    PyObject_HEAD
    ptwXPoints *ptwX;
} listOfDoubles_CPy;

#define is_listOfDoubles_CPyObject( v ) ((v)->ob_type == &listOfDoubles_CPyType)

static char listOfDoubles_C__doc__[] = 
    "The listOfDoubles_C class stores and manipulates a list of C double values.\n" \
    "\n" \
    "Constructor:\n" \
    "listOfDoubles_C( data = [], initialSize = 100 )\n\n" \
    "Constructor arguments are ([o] implies optional argument):\n" \
    "   data            [o] an iterable containing a list of values convertible to a C double.\n" \
    "   initialSize     [o] the initial size of the primary data cache (default = 100).\n";

static listOfDoubles_CPy *listOfDoubles_CNewInitialize( void );
static int listOfDoubles_C__init__( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static void listOfDoubles_C_dealloc( PyObject *self );
static PyObject *listOfDoubles_C__repr__( listOfDoubles_CPy *self );
static Py_ssize_t listOfDoubles_C__len__( listOfDoubles_CPy *self );
static PyObject *listOfDoubles_C__getitem__( listOfDoubles_CPy *self, Py_ssize_t index );
static void listOfDoubles_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 );
static PyObject *listOfDoubles_C__getslice__( listOfDoubles_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_ );
static int listOfDoubles_C__setitem__( listOfDoubles_CPy *self, Py_ssize_t index, PyObject *value );
static int listOfDoubles_C__setslice__( listOfDoubles_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *value );
static int listOfDoubles_C_contains( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C__add__( PyObject *self, PyObject *other );
static PyObject *listOfDoubles_C__iadd__( PyObject *self, PyObject *other );
static int listOfDoubles_C_addListNTimes( ptwXPoints *ptwX, int numberOfTimes, int64_t length, double *xs );
static PyObject *listOfDoubles_C__mul__( PyObject *self, PyObject *other );
static PyObject *listOfDoubles_C__imul__( PyObject *self, PyObject *other );
static PyObject *listOfDoubles_C__abs__( PyObject *self );
static PyObject *listOfDoubles_C__neg__( PyObject *self );
static PyObject *listOfDoubles_C_allocatedSize( listOfDoubles_CPy *self );
static PyObject *listOfDoubles_C_copy( listOfDoubles_CPy *self );
static PyObject *listOfDoubles_C_reallocatePoints( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_append( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C_count( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C_extend( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C_index( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_insert( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C_pop( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_remove( listOfDoubles_CPy *self, PyObject *args );
static PyObject *listOfDoubles_C_reverse( listOfDoubles_CPy *self );
static PyObject *listOfDoubles_C_setData( listOfDoubles_CPy *self, PyObject *args );
static int listOfDoubles_C_setData2( listOfDoubles_CPy *self, PyObject *PyXList );
static PyObject *listOfDoubles_C_sort( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_toString( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_toString2( listOfDoubles_CPy *self, int valuesPerLine, char *format, char *valueSeparator );
static char *listOfDoubles_C_toString_isFormatForDouble( char *format, int returnFormat );
static PyObject *listOfDoubles_C_unique( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *listOfDoubles_C_richCompare( PyObject *self, PyObject *other, int op );

static PyObject *listOfDoubles_C_createFromString( listOfDoubles_CPy *self, PyObject *args );

static int listOfDoubles_C_PyNumberToFloat( PyObject *n, double *d );
static int64_t listOfDoubles_C_doubleListToCList( PyObject *PyDoubleList, double **ds );
static listOfDoubles_CPy *listOfDoubles_C_copyPy( listOfDoubles_CPy *self, int asNew );
static void listOfDoubles_C_Get_listOfDoubles_CAsSelf( PyObject *self, PyObject *other, PyObject **s, PyObject **o );
static PyObject *listOfDoubles_C_GetNone( void );
static PyObject *listOfDoubles_C_SetPyErrorExceptionReturnNull( const char *s, ... );
static int listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... );
static int listOfDoubles_C_checkStatus( listOfDoubles_CPy *self, char const *func );

DL_EXPORT( void ) initlistOfDoubles_C( void );
/*
******************** listOfDoubles_CNewInitialize ************************
*/
static listOfDoubles_CPy *listOfDoubles_CNewInitialize( void ) {

    listOfDoubles_CPy *self = (listOfDoubles_CPy *) PyObject_New( listOfDoubles_CPy, &listOfDoubles_CPyType );

    if( self ) self->ptwX = NULL;
    return( self );
}
/*
************************ __init__ ************************************
*/
static int listOfDoubles_C__init__( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    int initialSize = 100;
    nfu_status status_nf;
    static char *kwlist[] = { "data", "initialSize" };
    ptwXPoints *ptwX;
    PyObject *dataPy = NULL;

    self->ptwX = NULL;
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|Oi", kwlist, &dataPy, &initialSize ) ) return( -1 );

    if( ( ptwX = ptwX_new( initialSize, &status_nf ) ) == NULL ) {
        PyErr_NoMemory( );
        return( -1 );
    }
    self->ptwX = ptwX;

    if( dataPy != NULL ) {
        if( listOfDoubles_C_setData2( self, dataPy ) != 0 ) {
            ptwX_free( ptwX );
            self->ptwX = NULL;
            return( -1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static void listOfDoubles_C_dealloc( PyObject *self ) {

    listOfDoubles_CPy *ptwX = (listOfDoubles_CPy *) self;

    if( ptwX->ptwX != NULL ) {
        ptwX_free( ptwX->ptwX );
    }
    self->ob_type->tp_free( self );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__repr__( listOfDoubles_CPy *self ) {

    if( listOfDoubles_C_checkStatus( self, "__repr__" ) != 0 ) return( NULL );
    return( listOfDoubles_C_toString2( self, 1, "%16.8e", " " ) );
}
/*
************************************************************
*/
static Py_ssize_t listOfDoubles_C__len__( listOfDoubles_CPy *self ) {

    if( listOfDoubles_C_checkStatus( self, "__len__" ) != 0 ) return( -1 );
    return( (Py_ssize_t) ptwX_length( self->ptwX ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__getitem__( listOfDoubles_CPy *self, Py_ssize_t index_ ) {

    int64_t index = (int64_t) index_;
    double *point;

    if( listOfDoubles_C_checkStatus( self, "__getitem__" ) != 0 ) return( NULL );
    if( index < 0 ) index += self->ptwX->length;
    if( ( point = ptwX_getPointAtIndex( self->ptwX, index ) ) == NULL ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( NULL );
    }
    return( (PyObject *) Py_BuildValue( "d", *point ) );
}
/*
************************************************************
*/
static void listOfDoubles_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 ) {

    if( length == 0 ) {
        *index1 = 0;
        *index2 = 0; }
    else {
        if( *index1 < 0 ) {
            *index1 += length;
            if( *index1 < 0 ) *index1 = 0;
        }
        if( *index1 > length ) *index1 = length;
        if( *index2 < 0 ) {
            *index2 = *index2 + length;
            if( *index2 < 0 ) *index2 = 0;
        }
        if( *index2 > length ) *index2 = length;
        if( *index1 > *index2 ) *index2 = *index1;
    }
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__getslice__( listOfDoubles_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_ ) {

    int64_t index1 = (int64_t) index1_, index2 = (int64_t) index2_;
    listOfDoubles_CPy *newPy;
    ptwXPoints *n;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( self, "__getslice__" ) != 0 ) return( NULL );

    listOfDoubles_C_getSliceIndices( self->ptwX->length, &index1, &index2 );
    if( ( newPy = listOfDoubles_CNewInitialize( ) ) == NULL ) return( NULL );
    if( ( n = ptwX_slice( self->ptwX, index1, index2, &status_nf ) ) == NULL ) {
        Py_DECREF( newPy );
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from __getslice__: %s", nfu_statusMessage( status_nf ) ) );
    }
    newPy->ptwX = n;
    return( (PyObject *) newPy );
}
/*
************************************************************
*/
static int listOfDoubles_C__setitem__( listOfDoubles_CPy *self, Py_ssize_t index_, PyObject *valuePy ) {
/*
*   This routine can be called either when "self[index_] = number" or "del self[index_]".
*/
    int status = 0;
    int64_t index = (int64_t) index_;

    if( listOfDoubles_C_checkStatus( self, "__setitem__" ) != 0 ) return( -1 );
    if( ( index < 0 ) || ( index >= self->ptwX->length ) ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( -1 );
    }
    if( valuePy == NULL ) {
        ptwX_deletePoints( self->ptwX, index, index + 1 ); }    /* Have already check for possible errors above, no need to here. */
    else {
        double value = 0;

        if( ( status = listOfDoubles_C_PyNumberToFloat( valuePy, &value ) ) == 0 ) {
            ptwX_setPointAtIndex( self->ptwX, index, value );       /* index verified above. */
        }
    }
    return( status );
}
/*
************************************************************
*/
static int listOfDoubles_C__setslice__( listOfDoubles_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *PyDoubleList ) {
/*
*   This routine can be called either when "self[index1_:index2_] = list of numbers" or "del self[index1_:index2_]".
*/
    int64_t index1 = (int64_t) index1_, index2 = (int64_t) index2_, length;
    double *xs;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( self, "__setslice_" ) != 0 ) return( -1 );

    listOfDoubles_C_getSliceIndices( self->ptwX->length, &index1, &index2 );

    ptwX_deletePoints( self->ptwX, index1, index2 );        /* Have already check for possible errors above, no need to here. */
    if( PyDoubleList != NULL ) {
        if( ( length = listOfDoubles_C_doubleListToCList( PyDoubleList, &xs ) ) < 0 ) return( -1 );
        if( length > 0 ) {
            status_nf = ptwX_insertPointsAtIndex( self->ptwX, index1, length, xs );
            free( xs );
            if( status_nf != nfu_Okay )
                return( listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( "Error from ptwX_insertPointsAtIndex: %s", nfu_statusMessage( status_nf ) ) );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static int listOfDoubles_C_contains( listOfDoubles_CPy *self, PyObject *args ) {

    int status = 0;
    int64_t index;
    double value = 0, difference;
    nfu_status status_nf;

    if( listOfDoubles_C_PyNumberToFloat( args, &value ) ) return( -1 );

    if( ( status_nf = ptwX_closesDifference( self->ptwX, value, &index, &difference ) ) != nfu_Okay ) {
        return( listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( "Error from ptwX_closesDifference: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( difference == 0. ) status = 1;
    return( status );
}
/*
************************************************************
*/
static PySequenceMethods listOfDoubles_CPy_sequence = {
    (lenfunc) listOfDoubles_C__len__,                       /* sq_length */
    0,                                                      /* sq_concat */
    0,                                                      /* sq_repeat */
    (ssizeargfunc) listOfDoubles_C__getitem__,              /* sq_item */
    (ssizessizeargfunc) listOfDoubles_C__getslice__,        /* sq_slice */
    (ssizeobjargproc) listOfDoubles_C__setitem__,           /* sq_ass_item */
    (ssizessizeobjargproc) listOfDoubles_C__setslice__,     /* sq_ass_slice */
    (objobjproc) listOfDoubles_C_contains,                  /* sq_contains */
    0,                                                      /* sq_inplace_concat */
    0                                                       /* sq_inplace_repeat */
};
/*
************************************************************
*/
static PyObject *listOfDoubles_C__add__( PyObject *self, PyObject *other ) {

    int status;
    int64_t length;
    double *xs = NULL;
    listOfDoubles_CPy *nPy = NULL;

    if( listOfDoubles_C_checkStatus( (listOfDoubles_CPy *) self, "__add__" ) != 0 ) return( NULL );
    if( ( length = listOfDoubles_C_doubleListToCList( other, &xs ) ) < 0 ) return( NULL );
    if( ( nPy = listOfDoubles_C_copyPy( (listOfDoubles_CPy *) self, 0 ) ) == NULL ) {
        free( xs );
        return( NULL );
    }

    if( length > 0 ) {
        status = listOfDoubles_C_addListNTimes( nPy->ptwX, 1, length, xs );
        free( xs );
        if( status != 0 ) {
            Py_DECREF( nPy );
            return( NULL );
        }
    }

    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__iadd__( PyObject *self, PyObject *other ) {

    int status;
    int64_t length;
    double *xs = NULL;
    ptwXPoints *ptwX = ((listOfDoubles_CPy *) self)->ptwX;

    if( listOfDoubles_C_checkStatus( (listOfDoubles_CPy *) self, "__iadd__" ) != 0 ) return( NULL );
    if( ( length = listOfDoubles_C_doubleListToCList( other, &xs ) ) < 0 ) return( NULL );

    if( length > 0 ) {
        status = listOfDoubles_C_addListNTimes( ptwX, 1, length, xs );
        free( xs );
        if( status != 0 ) return( NULL );
    }

    Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static int listOfDoubles_C_addListNTimes( ptwXPoints *ptwX, int numberOfTimes, int64_t length, double *xs ) {

    int64_t i1, size = numberOfTimes * length;
    nfu_status status_nf;

    if( size < 1 ) return( 0 );
    for( i1 = 0; i1 < numberOfTimes; i1++ ) {
        status_nf = ptwX_insertPointsAtIndex( ptwX, ptwX->length, length, xs );
        if( status_nf != nfu_Okay ) {
            return( listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( "Error from ptwX_insertPointsAtIndex: %s", nfu_statusMessage( status_nf ) ) );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__mul__( PyObject *self, PyObject *other ) {

    int status;
    listOfDoubles_CPy *nPy = NULL, *oPy;
    long multiplier;
    PyObject *self2, *other2;

    listOfDoubles_C_Get_listOfDoubles_CAsSelf( self, other, &self2, &other2 );
    if( listOfDoubles_C_checkStatus( (listOfDoubles_CPy *) self2, "__mul__" ) != 0 ) return( NULL );
    oPy = (listOfDoubles_CPy *) self2;
    multiplier = PyInt_AsLong( other2 );
    if( multiplier == -1 ) {
        if( PyErr_Occurred( ) != NULL ) return( NULL );
    }
    if( multiplier < 0 ) multiplier = 0;

    if( ( nPy = listOfDoubles_C_copyPy( (listOfDoubles_CPy *) self2, 1 ) ) == NULL ) return( NULL );

    if( oPy->ptwX->length > 0 ) {       /* BRB; the next line needs to be fixed so that multiplier is checked before casting or something???????? */
        status = listOfDoubles_C_addListNTimes( nPy->ptwX, (int) multiplier, oPy->ptwX->length, oPy->ptwX->points );
        if( status != 0 ) {
            Py_DECREF( nPy );
            return( NULL );
        }
    }

    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__imul__( PyObject *self, PyObject *other ) {

    listOfDoubles_CPy *oPy = (listOfDoubles_CPy *) self;
    long multiplier;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( (listOfDoubles_CPy *) self, "__imul__" ) != 0 ) return( NULL );
    multiplier = PyInt_AsLong( other );
    if( multiplier == -1 ) {
        if( PyErr_Occurred( ) != NULL ) return( NULL );
    }

    if( oPy->ptwX->length > 0 ) {
        if( multiplier < 1 ) {
            ptwX_clear( oPy->ptwX ); }
        else if( multiplier > 1 ) {
            if( ( status_nf = ptwX_reallocatePoints( oPy->ptwX, multiplier * oPy->ptwX->length, 0 ) ) != nfu_Okay )
                return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_reallocatePoints: %s", nfu_statusMessage( status_nf ) ) );
                    /* BRB; the next line needs to be fixed so that multiplier is checked before casting or something???????? */
            if( listOfDoubles_C_addListNTimes( oPy->ptwX, (int) multiplier - 1, oPy->ptwX->length, oPy->ptwX->points ) != 0 ) return( NULL );
        }
    }
    Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__abs__( PyObject *self ) {

    nfu_status status_nf;
    listOfDoubles_CPy *nPy;

    if( ( nPy = listOfDoubles_C_copyPy( (listOfDoubles_CPy *) self, 0 ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwX_abs( ((listOfDoubles_CPy *) nPy)->ptwX ) ) != nfu_Okay ) {
        Py_DECREF( nPy );
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_abs: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C__neg__( PyObject *self ) {

    nfu_status status_nf;
    listOfDoubles_CPy *nPy;

    if( ( nPy = listOfDoubles_C_copyPy( (listOfDoubles_CPy *) self, 0 ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwX_neg( ((listOfDoubles_CPy *) nPy)->ptwX ) ) != nfu_Okay ) {
        Py_DECREF( nPy );
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_neg: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyNumberMethods listOfDouble_CPy_number = {
    (binaryfunc) listOfDoubles_C__add__,    /* nb_add */
    0,                                      /* nb_subtract */
    (binaryfunc) listOfDoubles_C__mul__,    /* nb_multiply */
    0,                                      /* nb_divide */
    0,                                      /* nb_remainder */
    0,                                      /* binaryfunc nb_divmod */
    0,                                      /* nb_power */
    (unaryfunc) listOfDoubles_C__neg__,     /* nb_negative */
    0,                                      /* unaryfunc nb_positive */
    (unaryfunc) listOfDoubles_C__abs__,     /* nb_absolute */
    0,                                      /* inquiry nb_nonzero */
    0,                                      /* unaryfunc nb_invert */
    0,                                      /* binaryfunc nb_lshift */
    0,                                      /* binaryfunc nb_rshift */
    0,                                      /* binaryfunc nb_and */
    0,                                      /* binaryfunc nb_xor */
    0,                                      /* binaryfunc nb_or */
    0,                                      /* nb_coerce */
    0,                                      /* nb_int */
    0,                                      /* nb_long */
    0,                                      /* nb_float */
    0,                                      /* unaryfunc nb_oct */
    0,                                      /* unaryfunc nb_hex */
    (binaryfunc) listOfDoubles_C__iadd__,   /* nb_inplace_add */            /* Added in release 2.0 */
    0,                                      /* nb_inplace_subtract */
    (binaryfunc) listOfDoubles_C__imul__,   /* nb_inplace_multiply */
    0,                                      /* nb_inplace_divide */
    (binaryfunc) 0,                         /* nb_inplace_remainder */
    (ternaryfunc) 0,                        /* nb_inplace_power */
#if 0
    binaryfunc nb_inplace_lshift
    binaryfunc nb_inplace_rshift
    binaryfunc nb_inplace_and
    binaryfunc nb_inplace_xor
    binaryfunc nb_inplace_or
    binaryfunc nb_floor_divide         /* The following added in release 2.2, and require the Py_TPFLAGS_HAVE_CLASS flag */
    binaryfunc nb_true_divide
    binaryfunc nb_inplace_floor_divide
    binaryfunc nb_inplace_true_divide
    unaryfunc nb_index                 /* The following added in release 2.5 */
#endif
};

/*
************************************************************
*/
static PyObject *listOfDoubles_C_allocatedSize( listOfDoubles_CPy *self ) {

    return( (PyObject *) Py_BuildValue( "l", self->ptwX->allocatedSize ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_copy( listOfDoubles_CPy *self ) {

    nfu_status status_nf;
    ptwXPoints *n;
    listOfDoubles_CPy *nPy;

    n = ptwX_clone( self->ptwX, &status_nf );
    if( n == NULL ) return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from copy: %s", nfu_statusMessage( status_nf ) ) );
    if( ( nPy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( n );
        return( NULL );
    }
    nPy->ptwX = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_reallocatePoints( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    int size;
    int forceSmallerResize = 1;
    static char *kwlist[] = { "size", "forceSmaller", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "i|i", kwlist, &size, &forceSmallerResize ) ) return( NULL );

    if( ( status_nf = ptwX_reallocatePoints( self->ptwX, size, forceSmallerResize ) ) != nfu_Okay )
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from reallocatePoints: %s", nfu_statusMessage( status_nf ) ) );

    return( listOfDoubles_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_append( listOfDoubles_CPy *self, PyObject *args ) {

    double value;

    if( listOfDoubles_C_checkStatus( self, "append" ) != 0 ) return( NULL );
    if( !PyArg_ParseTuple( args, "d", &value ) ) return( NULL );

    ptwX_setPointAtIndex( self->ptwX, self->ptwX->length, value );
    return( listOfDoubles_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_count( listOfDoubles_CPy *self, PyObject *args ) {

    int count;
    double value;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( self, "count" ) != 0 ) return( NULL );
    if( !PyArg_ParseTuple( args, "d", &value ) ) return( NULL );

    status_nf = ptwX_countOccurrences( self->ptwX, value, &count );
    return( Py_BuildValue( "i", (int) count ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_extend( listOfDoubles_CPy *self, PyObject *args ) {

    int64_t length;
    PyObject *PyDoubleList;
    double *xs;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( self, "extend" ) != 0 ) return( NULL );
    if( !PyArg_ParseTuple( args, "O", &PyDoubleList ) ) return( NULL );

    if( ( length = listOfDoubles_C_doubleListToCList( PyDoubleList, &xs ) ) < 0 ) return( NULL );

    if( length  > 0 ) {
        status_nf = ptwX_insertPointsAtIndex( self->ptwX, self->ptwX->length, length, xs );
        free( xs );
        if( status_nf != nfu_Okay ) 
            return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_insertPointsAtIndex: %s", nfu_statusMessage( status_nf ) ) );
    }

    return( listOfDoubles_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_index( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    int iStart = 0, iStop = (int) self->ptwX->length;       /* BRB; this needs to be fixed so that  length is checked before casting???????? */
    int64_t index;
    double value, difference;
    static char *kwlist[] = { "start", "stop", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "d|ii", kwlist, &value, &iStart, &iStop ) ) return( NULL );
    ptwX_closesDifferenceInRange( self->ptwX, iStart, iStop, value, &index, &difference );
    if( difference != 0 ) {
        PyErr_SetString( PyExc_ValueError, "value not in self" );
        return( NULL );
    }
    return( Py_BuildValue( "i", (int) index ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_insert( listOfDoubles_CPy *self, PyObject *args ) {

    int index;
    double value;
    nfu_status status_nf;

    if( listOfDoubles_C_checkStatus( self, "insert" ) != 0 ) return( NULL );
    if( !PyArg_ParseTuple( args, "id", &index, &value ) ) return( NULL );

    if( index < 0 ) index += self->ptwX->length;    /* Next few lines attempt to replicate weird logic of python's [] class. */
    if( index < 0 ) index = 0;
    if( index > self->ptwX->length ) index = (int) self->ptwX->length; /* BRB; this needs to be fixed: see above??????? */

    status_nf = ptwX_insertPointsAtIndex( self->ptwX, index, 1, &value );
    if( status_nf != nfu_Okay ) 
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_insertPointsAtIndex: %s", nfu_statusMessage( status_nf ) ) );
    return( listOfDoubles_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_pop( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    int index = -1;
    static char *kwlist[] = { "index", NULL };
    double value;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &index ) ) return( NULL );
    if( index < 0 ) index += self->ptwX->length;
    if( ( index < 0 ) || ( index >= self->ptwX->length ) ) {
        PyErr_SetString( PyExc_IndexError, "pop index out of range" );
        return( NULL );
    }

    value = *ptwX_getPointAtIndex( self->ptwX, index );
    ptwX_deletePoints( self->ptwX, index, index + 1 );

    return( Py_BuildValue( "d", value ) );
}   
/*
************************************************************
*/
static PyObject *listOfDoubles_C_remove( listOfDoubles_CPy *self, PyObject *args ) {

    int64_t index;
    double value, difference;

    if( !PyArg_ParseTuple( args, "d", &value ) ) return( NULL );
    ptwX_closesDifference( self->ptwX, value, &index, &difference );

    if( difference != 0 ) {
        PyErr_SetString( PyExc_ValueError, "value not in self" );
        return( NULL );
    }

    ptwX_deletePoints( self->ptwX, index, index + 1 );
    return( listOfDoubles_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_reverse( listOfDoubles_CPy *self ) {

    nfu_status status_nf;
    listOfDoubles_CPy *reversedPy;
    ptwXPoints *reversed;

    if( ( reversed = ptwX_clone( self->ptwX, &status_nf ) ) == NULL ) {
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_clone: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( ( status_nf = ptwX_reverse( reversed ) ) != nfu_Okay ) {
        ptwX_free( reversed );
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_reverse: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( ( reversedPy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( reversed );
        return( NULL );
    }
    reversedPy->ptwX = reversed;
    return( (PyObject *) reversedPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_setData( listOfDoubles_CPy *self, PyObject *args ) {

    PyObject *status = NULL, *PyXList;

    if( !PyArg_ParseTuple( args, "O", &PyXList ) ) return( NULL );
    if( listOfDoubles_C_setData2( self, PyXList ) == 0 ) status = listOfDoubles_C_GetNone( );
    return( status );

}
/*
************************************************************
*/
static int listOfDoubles_C_setData2( listOfDoubles_CPy *self, PyObject *PyXList ) {

    nfu_status status_nf;
    int64_t length;
    double *xs;

    ptwX_setData( self->ptwX, 0, NULL ); 
    if( ( length = listOfDoubles_C_doubleListToCList( PyXList, &xs ) ) < 0 ) return( -1 );
    if( length > 0 ) {
        status_nf = ptwX_setData( self->ptwX, length, xs );
        free( xs );
        if( status_nf != nfu_Okay )
            return( listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( "Error from ptwX_setData: %s", nfu_statusMessage( status_nf ) ) );
    }

    return( 0 );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_sort( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    int descendingOrder = 1;
    static char *kwlist[] = { "reverse", NULL };
    ptwXPoints *sorted;
    listOfDoubles_CPy *sortedPy;
    enum ptwX_sort_order order = ptwX_sort_order_ascending;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &descendingOrder ) ) return( NULL );

    if( ( sorted = ptwX_clone( self->ptwX, &status_nf ) ) == NULL ) {
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_clone: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( descendingOrder ) order = ptwX_sort_order_descending;
    if( ( status_nf = ptwX_sort( sorted, order ) ) != nfu_Okay ) {
        ptwX_free( sorted );
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_sort: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( sortedPy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( sorted );
        return( NULL );
    }
    sortedPy->ptwX = sorted;
    return( (PyObject *) sortedPy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_toString( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    int valuesPerLine = 1;
    char *format = "%16.8e", *valueSeparator = " ";
    static char *kwlist[] = { "valuesPerLine", "format", "valueSeparator", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|iss", kwlist, &valuesPerLine, &format, &valueSeparator ) ) return( NULL );
    return( listOfDoubles_C_toString2( self, valuesPerLine, format, valueSeparator ) );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_toString2( listOfDoubles_CPy *self, int valuesPerLine, char *format, char *valueSeparator ) {

    int lineFeed, extraSpace = 4, dummyLen;     /* '\n' plus three more for reserve. */
    ptwXPoints *ptwX = self->ptwX;
    int64_t index, strLength, valueSeparatorLength = strlen( valueSeparator );
    char *s, *e, *p, dummy[1024];
    PyObject *str;

    if( listOfDoubles_C_checkStatus( self, "toString2" ) != 0 ) return( NULL );

    s = listOfDoubles_C_toString_isFormatForDouble( format, 0 );
    if( s == NULL ) return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "invalid format = '%s' for converting a double to a string", format ) );
    s = listOfDoubles_C_toString_isFormatForDouble( s, 0 );
    if( s != NULL ) return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "invalid format = '%s' for converting a double to a string", format ) );

    dummyLen = snprintf( dummy, 0, format, -3.1238972718972345e-21, -6.1432987987397648974e-33 );
    if( dummyLen < 0 ) return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Invalid format string = '%s'", format ) );
    strLength = ptwX->length * ( dummyLen + valueSeparatorLength + extraSpace );

    if( ( s = (char *) malloc( (size_t) strLength * sizeof( char ) ) ) == NULL ) {
        PyErr_NoMemory( );
        return( NULL );
    }
    s[0] = 0;
    for( index = 0, e = s, lineFeed = valuesPerLine; index < ptwX->length; index++ ) {
        sprintf( e, format, ptwX->points[index] );
        while( *e ) e++;
        if( index < ptwX->length - 1 ) {
            for( p = valueSeparator; *p; p++, e++ ) *e = *p;
            *e = 0;
        }
        lineFeed--;
        if( lineFeed <= 0 ) {
            *e = '\n';
            e++;
            *e = 0;
            lineFeed = valuesPerLine;
        }
    }
    if( lineFeed != valuesPerLine ) {
        *e = '\n';
        e++;
        *e = 0;
    }
    str = PyString_FromString( s );
    free( s );
    return( str );
}
/*
************************************************************
*/
static char *listOfDoubles_C_toString_isFormatForDouble( char *format, int returnFormat ) {

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
static PyObject *listOfDoubles_C_unique( listOfDoubles_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    int order = 0;
    static char *kwlist[] = { "order", NULL };
    ptwXPoints *unique;
    listOfDoubles_CPy *uniquePy;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &order ) ) return( NULL );

    if( ( unique = ptwX_unique( self->ptwX, order, &status_nf ) ) == NULL ) {
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_unique: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( uniquePy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( unique );
        return( NULL );
    }
    uniquePy->ptwX = unique;
    return( (PyObject *) uniquePy );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_richCompare( PyObject *self, PyObject *other, int op ) {

    int comparison;
    listOfDoubles_CPy *self2 = (listOfDoubles_CPy *) self, *other2 = (listOfDoubles_CPy *) other;
    nfu_status status_nf;
    PyObject *result = Py_False;

    if( PyObject_TypeCheck( other, &listOfDoubles_CPyType ) ) {

        if( ( status_nf = ptwX_compare( self2->ptwX, other2->ptwX, &comparison ) ) != nfu_Okay )
            return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_compare: %s", nfu_statusMessage( status_nf ) ) );

        switch( op ) {
        case Py_LT :
            if( comparison < 0 ) result = Py_True;
            break;
        case Py_LE :
            if( comparison <= 0 ) result = Py_True;
            break;
        case Py_EQ :
            if( comparison == 0 ) result = Py_True;
            break;
        case Py_NE :
            if( comparison != 0 ) result = Py_True;
            break;
        case Py_GE :
            if( comparison >= 0 ) result = Py_True;
            break;
        case Py_GT :
            if( comparison > 0 ) result = Py_True;
            break;
        }
    }

    Py_XINCREF(result);
    return( result );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_createFromString( listOfDoubles_CPy *self, PyObject *args ) {

    nfu_status status_nf;
    ptwXPoints *ptwX;
    listOfDoubles_CPy *ptwXPy;
    char *str, *endCharacter;

    if( !PyArg_ParseTuple( args, "s", &str ) ) return( NULL );

    if( ( ptwX = ptwX_fromString( str, &endCharacter, &status_nf ) ) == NULL )
        return( listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from ptwX_fromString: %s", nfu_statusMessage( status_nf ) ) );

    if( ( ptwXPy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( ptwX );
        return( NULL );
    }
    ptwXPy->ptwX = ptwX;
    return( Py_BuildValue( "(O,s)", ptwXPy, endCharacter ) );
}
/*
************************************************************
*/
static int listOfDoubles_C_PyNumberToFloat( PyObject *valuePy, double *value ) {

    PyObject *floatPy = PyNumber_Float( valuePy );

    if( floatPy == NULL ) {
        return( listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( "could not convert object to float" ) );
    }
    *value = PyFloat_AsDouble( floatPy );
    return( 0 );
}
/*
************************************************************
*/
static int64_t listOfDoubles_C_doubleListToCList( PyObject *PyDoubleList, double **xs ) {

    int status = 0;
    int64_t length, i1;
    double *d1;
    PyObject *item, *iterator;

    *xs = NULL;
    if( ( iterator = PyObject_GetIter( PyDoubleList ) ) == NULL ) return( -1 );
    if( ( length = (int64_t) PySequence_Size( PyDoubleList ) ) != (int64_t) 0 ) {
        if( ( *xs = (double *) malloc( (size_t) length * sizeof( double ) ) ) == NULL ) {
            PyErr_NoMemory( );
            Py_DECREF( iterator );
            return( -1 );
        }
        for( i1 = 0, d1 = *xs, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i1++, d1++ ) {
            if( ( status = listOfDoubles_C_PyNumberToFloat( item, d1 ) ) != 0 )
                listOfDoubles_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i1 );
            Py_DECREF( item );
        }
    }
    Py_DECREF( iterator );
    if( status != 0 ) {
        length = -1;
        free( *xs );
        *xs = NULL;
    }

    return( length );
}
/*
************************************************************
*/
static listOfDoubles_CPy *listOfDoubles_C_copyPy( listOfDoubles_CPy *self, int asNew ) {

    nfu_status status_nf;
    ptwXPoints *n;
    listOfDoubles_CPy *nPy;

    if( asNew ) {
        n = ptwX_new( self->ptwX->length, &status_nf ); }
    else {
        n = ptwX_clone( self->ptwX, &status_nf );
    }
    if( n == NULL ) return( (listOfDoubles_CPy *) listOfDoubles_C_SetPyErrorExceptionReturnNull( "Error from copy: %s", nfu_statusMessage( status_nf ) ) );
    if( ( nPy = listOfDoubles_CNewInitialize( ) ) == NULL ) {
        ptwX_free( n );
        return( NULL );
    }
    nPy->ptwX = n;
    return( nPy );
}
/*
************************************************************
*/
static void listOfDoubles_C_Get_listOfDoubles_CAsSelf( PyObject *self, PyObject *other, PyObject **s, PyObject **o ) {

    *s = self;
    *o = other;
    if( !PyObject_TypeCheck( self, &listOfDoubles_CPyType ) ) {
        *s = other;
        *o = self;
    }
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_GetNone( void ) {

    Py_INCREF( Py_None );
    return( Py_None );
}
/*
************************************************************
*/
static PyObject *listOfDoubles_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

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
static int listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... ) {

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
static int listOfDoubles_C_checkStatus( listOfDoubles_CPy *self, char const *func ) {

    if( self->ptwX->status == nfu_Okay ) return( 0 );
    listOfDoubles_C_SetPyErrorExceptionReturnMinusOne( 
        "listOfDoubles_C object had a prior error and is not usable in %s: status = %d", func, self->ptwX->status );
    return( -1 );
}
/*
************************************************************
*/
static PyMethodDef listOfDoubles_CPyMethods[] = {

    { "allocatedSize", (PyCFunction) listOfDoubles_C_allocatedSize, METH_NOARGS, "Returns the size of memory allocated in the points region." },
    { "append", (PyCFunction) listOfDoubles_C_append, METH_VARARGS, "appends the argument to self.\n" \
        "\nArguments are:\n" \
        "   value   float value to append." },
    { "copy", (PyCFunction) listOfDoubles_C_copy, METH_NOARGS, "Returns a copy of self." },
    { "count", (PyCFunction) listOfDoubles_C_count, METH_VARARGS, "Returns the number of times value occurs in self.\n" \
        "\nArguments are:\n" \
        "   value   float value to count." },
    { "extend", (PyCFunction) listOfDoubles_C_extend, METH_VARARGS, "Adds the list of float to the end of self.\n" \
        "\nArguments are:\n" \
        "   list    list of floats." },
    { "index", (PyCFunction) listOfDoubles_C_index, METH_VARARGS | METH_KEYWORDS, "Returns the index of the first occurs of value.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   value       the float value to search for.\n" \
        "   start [o]   the starting index for the search [default is 0].\n" \
        "   stop [o]    the last index for the search [default is -1]." },
    { "insert ", (PyCFunction) listOfDoubles_C_insert, METH_VARARGS, "Insert the float before position index.\n" \
        "\nArguments are:\n" \
        "   index       value is insert before this index.\n" \
        "   value       float value to insert." },
    { "pop", (PyCFunction) listOfDoubles_C_pop, METH_NOARGS, "Removes the float in index and returns it.\n" \
        "   index       index of value to pop (default is -1).\n" },
    { "reallocatePoints", (PyCFunction) listOfDoubles_C_reallocatePoints, METH_VARARGS | METH_KEYWORDS, \
        "self.reallocatePoints( size, forceSmaller = True )\n\n" \
        "Adjusts the memory allocated for primary points to the maximum of size and the current length of self.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   size              the desired allocated size of self (actual size will be larger if length is greater than size),\n" \
        "   forceSmaller  [o] if False action is only taken if allocated is significantly greater than size (default is True).\n" },
    { "remove", (PyCFunction) listOfDoubles_C_remove, METH_NOARGS, "Removes the first occurs of value in self.\n" \
        "\nArguments are:\n" \
        "   value       float value to remove." },
    { "reverse", (PyCFunction) listOfDoubles_C_reverse, METH_NOARGS, "Retuns a listOfDoubles_C instance whose list is the reversed of self." },
    { "setData", (PyCFunction) listOfDoubles_C_setData, METH_VARARGS, "Replaces the data in self with the frist argument. This argument must be a list. " \
        "Each item of the list must contain two floats, or objects that can be convert to float (e.g., [ [ 1, 2 ], [ 2, 4 ], [ 4, 0.5 ] ]" },
    { "sort", (PyCFunction) listOfDoubles_C_sort, METH_VARARGS | METH_KEYWORDS, "Retuns a listOfDoubles_C instance whose list is self's list sorted. " \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   reverse [o]    if True, list is in descending order, otherwise it is in ascending order (default is False)" },
    { "toString", (PyCFunction) listOfDoubles_C_toString, METH_VARARGS | METH_KEYWORDS, "Returns a string representation of self." \
        " This method has three keyword parameters:\nvaluePerLine, format and valueSeparator which are defined as,\n" \
        "    valuesPerLine      the number of values to put on each line\n" \
        "    format             a valid format to convert a value (i.e., a floats) into a string (e.g. format = ' %.3f')\n" \
        "    valueSeparator     a string to put between every value (e.g, to put a comma to separate values use valueSeparator = ',')" },
    { "unique", (PyCFunction) listOfDoubles_C_unique, METH_VARARGS | METH_KEYWORDS, "Retuns a listOfDoubles_C instance whose list contains self's list sorted. " \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   order [o]      if < 0 order is descending, if > 0 order is ascending, otherwise, order is the same as in self (default is 0)" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
static PyTypeObject listOfDoubles_CPyType = {
    PyObject_HEAD_INIT( NULL )
    0,                                          /* ob_size        */
    "listOfDoubles_C.listOfDoubles_C",          /* tp_name        */
    sizeof( listOfDoubles_CPy ),                /* tp_basicsize   */
    0,                                          /* tp_itemsize    */
    /* methods */ 
    (destructor) listOfDoubles_C_dealloc,       /* tp_dealloc     */
    0,                                          /* tp_print       */
    0,                                          /* tp_getattr     */
    0,                                          /* tp_setattr     */
    0,                                          /* tp_compare     */
    (reprfunc) listOfDoubles_C__repr__,         /* tp_repr        */
    &listOfDouble_CPy_number,                   /* tp_as_number   */
    &listOfDoubles_CPy_sequence,                /* tp_as_sequence */
    0,                                          /* tp_as_mapping  */
    0,                                          /* tp_hash        */
    0,                                          /* tp_call        */
    0,                                          /* tp_str         */
    0,                                          /* tp_getattro    */
    0,                                          /* tp_setattro    */
    0,                                          /* tp_as_buffer   */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, /* tp_flags     */
    listOfDoubles_C__doc__,                     /* tp_doc         */
    0,                                          /* tp_traverse    */
    0,                                          /* tp_clear       */
    listOfDoubles_C_richCompare,                /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter        */
    0,                                          /* tp_iternext    */
    listOfDoubles_CPyMethods,                   /* tp_methods     */
    0,                                          /* tp_members     */
    0,                                          /* tp_getset      */
    0,                                          /* tp_base        */
    0,                                          /* tp_dict        */
    0,                                          /* tp_descr_get   */
    0,                                          /* tp_descr_set   */
    0,                                          /* tp_dictoffset  */
    (initproc) listOfDoubles_C__init__,         /* tp_init        */
    0,                                          /* tp_alloc       */
    0                                           /* tp_new         */
};
/*
************************************************************
*/

static PyMethodDef listOfDoubles_CMiscPyMethods[] = {

    { "createFromString", (PyCFunction) listOfDoubles_C_createFromString, METH_VARARGS,
        "createFromString( str )\n\n" \
        "Returns a tuple of two elements. The first element is a listOfDoubles_C instance representing the float values\n" \
        "translated from 'str'. The second element is the portion of 'str' not translated\n" \
        "\nArguments are:\n" \
        "   str           The string containing a list of floats to be converted.\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
DL_EXPORT( void ) initlistOfDoubles_C( void ) {

    PyObject *m;

    listOfDoubles_CPyType.tp_new = PyType_GenericNew;
    if( PyType_Ready( &listOfDoubles_CPyType ) < 0 ) return;

    if( ( m = Py_InitModule3( "listOfDoubles_C", listOfDoubles_CMiscPyMethods, "A module that contains the class listOfDoubles_C." ) ) == NULL ) return;

    Py_INCREF( &listOfDoubles_CPyType );
    PyModule_AddObject( m, "listOfDoubles_C", (PyObject *) &listOfDoubles_CPyType );

}
