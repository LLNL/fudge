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

#define POINTWISEXY_C_MODULE 1
#include <ptwXY.h>
#include <pointwiseXY_C.h>

#if( PY_VERSION_HEX < 0x02050000 )
#define Py_ssize_t int64_t 
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#define ssizeobjargproc intobjargproc
#define ssizessizeobjargproc intintobjargproc
#endif

#define pyObject_Unsupported 0
#define pyObject_Number 1
#define pyObject_List 2
#define pyObject_Tuple 3
#define pyObject_ptwXY 4

enum e_interpolationType { e_interpolationTypeInvalid = -1, e_interpolationTypeLinear, e_interpolationTypeLog, e_interpolationTypeFlat, e_interpolationTypeOther };

typedef nfu_status (*ptwXY_ptwXY_d_func)( ptwXYPoints *, double );
typedef nfu_status (*ptwXY_ptwXY_d_i_func)( ptwXYPoints *, double, int );
typedef ptwXYPoints *(*ptwXY_ptwXY_ptwXY_func)( ptwXYPoints *, ptwXYPoints *, nfu_status * );

static double _defaultAccuracy = 1e-3;

staticforward PyTypeObject pointwiseXY_CPyType;

typedef struct pointwiseXY_CPy_s {
    PyObject_HEAD
    ptwXYPoints *ptwXY;
    int infill;
    int safeDivide;
} pointwiseXY_CPy;

#define is_pointwiseXY_CPyObject( v ) ((v)->ob_type == &pointwiseXY_CPyType)

static char pointwiseXY_C__doc__[] = 
    "The pointwiseXY_C class stores and manipulates a list of XY points (i.e., [x, y] pairs).\n" \
    "Methods to add, substract, multiply and divide a pointwiseXY_C object with a scaler\n" \
    "(i.e., a number) or other pointwiseXY_C object are provided.\n" \
    "\n" \
    "Constructor:\n" \
    "pointwiseXY_C( data = [], dataForm = 'xys', initialSize = 100, overflowSize = 10, accuracy = defaultAccuracy( ), biSectionMax = 3." \
    ", interpolation = 'linear,linear', infill = True, safeDivide = False, userFlag = 0 )\n\n" \
    "Constructor arguments are ([o] implies optional argument):\n" \
    "   data            [o] the [x_i, y_i] pairs given as described by the dataForm argument,\n" \
    "   dataForm        [o] can be one of three strings (case is ignored) that describes the form of data:\n" \
    "                          'XYs'      data are a sequence of [x_i, y_i] pairs with x_i ascending (i.e., x_i < x_{i+1})\n" \
    "                                          (e.g., data = [ [ 1, 1 ], [ 2.3, 2 ], [ 3.4, 6 ], [ 5.1, 4.3 ] ])\n" \
    "                          'XsAndYs'  data are given as [xs, ys] where xs the list of x-values and ys the matching list of y-values,\n" \
    "                                          (e.g., data = [ [ 1, 2.3, 3.4, 5.1 ], [ 1, 2, 6, 4.3 ] ]),\n" \
    "                          'list'     data are given as [x_0, y_0, x_1, y_1, ..., x_n, y_n],\n" \
    "                                          (e.g., data = [ 1, 1, 2.3, 2, 3.4, 6, 5.1, 4.3 ]),\n" \
    "   initialSize     [o] the initial size of the primary data cache (default = 100),\n" \
    "   overflowSize    [o] the initial size of the secondary (overflow) data cache (default = 10),\n" \
    "   accuracy        [o] the accuracy of the data for the given interpolation (default is value returned by function defaultAccuracy)\n" \
    "   biSectionMax    [o] at times (e.g., multiplication), points may need to be added to maintain the give accuracy. In this case, a region\n" \
    "                       is continuously divided into two until the accuracy is met or biSectionMax divisions have occurred (default = 3),\n" \
    "   interpolation   [o] can be one of the following strings:\n" \
    "                          'flat'              for the domain [x_i,x_{i+1}), the y-value is y_i\n" \
    "                          'linear,flat'       same as 'flat'\n" \
    "                          'log,flat'          same as 'flat'\n" \
    "                          'linear,linear'\n" \
    "                          'linear,log'\n" \
    "                          'log,linear'\n" \
    "                          'log,log'\n" \
    "                          'other'             The rest are equivalent to this sting.\n" \
    "                          'other,*'           the y interpolation can be any string.\n" \
    "                          'linear,other'\n" \
    "                          'log,other'\n" \
    "   infill          [o] if True, multiplication will continuously divide a region until accuracy or biSectionMax is met (default = True),\n" \
    "   safeDivide      [o] if True, safe division is used (default = False)\n" \
    "   userFlag        [o] an integer (of type C int) used to store a user defined flag. (see getUserFlag and setUserFlag).\n";

static pointwiseXY_CPy *pointwiseXY_CNewInitialize( int infill, int safeDivide );
static int pointwiseXY_C__init__( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static void pointwiseXY_C_dealloc( PyObject *self );
static PyObject *pointwiseXY_C__repr__( pointwiseXY_CPy *self );
static Py_ssize_t pointwiseXY_C__len__( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C__getitem__( pointwiseXY_CPy *self, Py_ssize_t index );
static PyObject *pointwiseXY_C__getslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_ );
static int pointwiseXY_C__setitem__( pointwiseXY_CPy *self, Py_ssize_t index, PyObject *value );
static int pointwiseXY_C__setslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *value );
static void pointwiseXY_C_Get_pointwiseXY_CAsSelf( PyObject *self, PyObject *other, PyObject **s, PyObject **o );
static PyObject *pointwiseXY_C__add__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__iadd__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__sub__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__isub__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__mul__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__imul__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__div__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__idiv__( PyObject *self, PyObject *other );
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status_nf );
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status_nf );
static PyObject *pointwiseXY_C__mod__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C_add_sub_mul_div_number( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func, const char *oper, nfu_status *status_nf );
static PyObject *pointwiseXY_C_add_sub_mul_div_number_insitu( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func, char const *oper );
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( pointwiseXY_CPy *self, pointwiseXY_CPy *other, ptwXY_ptwXY_ptwXY_func func, const char *oper,
    nfu_status *status_nf );
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( pointwiseXY_CPy *self, pointwiseXY_CPy *other, ptwXY_ptwXY_ptwXY_func func, char const *oper );
static PyObject *pointwiseXY_C__pow__( PyObject *self, PyObject *other, PyObject *dummy );
static PyObject *pointwiseXY_C__neg__( PyObject *self );
static PyObject *pointwiseXY_C__abs__( PyObject *self );
static PyObject *pointwiseXY_C_allocatedSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_applyFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static nfu_status pointwiseXY_C_applyFunction2( ptwXYPoint *ptwXY, void *argList );
static PyObject *pointwiseXY_C_changeInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_changeInterpolationIfNeeded( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_changeInterpolation2( pointwiseXY_CPy *self, ptwXY_interpolation interpolation, double accuracy, double lowerEps, double upperEps );
static PyObject *pointwiseXY_C_clip( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_coalescePoints( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_convolute( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_copy( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_copyDataToXYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_copyDataToXsAndYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_dullEdges( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_exp( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_getAccuracy( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getBiSectionMax( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getInfill( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getInterpolation( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getSecondaryCacheSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getSafeDivide( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getValue( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_getUserFlag( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_integrate( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_integrateWithWeight_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_integrateWithWeight_sqrt_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_normalize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_groupOneFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupTwoFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupThreeFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupFunctionsCommon( pointwiseXY_CPy *f1, PyObject *f2, PyObject *f3, PyObject *groupBoundariesPy, PyObject *normPy );
static PyObject *pointwiseXY_C_areDomainsMutual( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_mutualify( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_overflowAllocatedSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_overflowLength( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_plot( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_reallocatePoints( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_reallocateOverflowPoints( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setData( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setData2( pointwiseXY_CPy *self, PyObject *PyXYList );
static int pointwiseXY_C_setDataFromPtwXY( pointwiseXY_CPy *self, pointwiseXY_CPy *otherPY );
static PyObject *pointwiseXY_C_setAccuracy( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setBiSectionMax( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setSecondaryCacheSize( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setDataFromList( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setDataFromList2( pointwiseXY_CPy *self, PyObject *list );
static PyObject *pointwiseXY_C_setDataFromXsAndYs( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setDataFromXsAndYs2( pointwiseXY_CPy *self, PyObject *PyXs, PyObject *PyYs );
static PyObject *pointwiseXY_C_setInfill( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setSafeDivide( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setValue( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setUserFlag( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_showInteralStructure( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_thicken( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_thin( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_trim( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_toString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_toString2( pointwiseXY_CPy *self, int pairsPerLine, char *format, char *pairSeparator );
static char *pointwiseXY_C_toString_isFormatForDouble( char *format, int returnFormat );
static PyObject *pointwiseXY_C_union(  pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_mergeClosePoints(  pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_getDomainGrid( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_xMin( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_xMax( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_xSlice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_yMin( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_yMax( pointwiseXY_CPy *self );

static PyObject *pointwiseXY_C_createFromFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static nfu_status pointwiseXY_C_createFromFunction2( double x, double *y, void *argList );
static PyObject *pointwiseXY_C_createFromString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );

static PyObject *pointwiseXY_C_ptwXY_interpolatePoint( PyObject *self, PyObject *args );

static PyObject *pointwiseXY_C_gaussian( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_basicGaussian( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_unitbaseInterpolate( pointwiseXY_CPy *self, PyObject *args );

static void pointwiseXY_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 );
static int pointwiseXY_C_PythonXYPairToCPair( PyObject *XYPairPy, double *x, double *y, int64_t index );
static int64_t pointwiseXY_C_PythonXYListToCList( PyObject *PyXYList, double **xys );
static int64_t pointwiseXY_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds, int ascending );
static ptwXPoints *pointwiseXY_C_PyFloatList_to_ptwXPoints( PyObject *PyFloatList );
static PyObject *pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwXPoints *ptwX );
static int pointwiseXY_C_PyNumberToFloat( PyObject *n, double *d );
static int pyObject_NumberOrPtwXY(  PyObject *other );
static int pointwiseXY_C_addedItemToPythonList( PyObject *list, PyObject *item );
static int isOkayAndHasData( ptwXYPoints *ptwXY );
static int pointwiseXY_C_getInterpolationFromString( char *interpolationStr, ptwXY_interpolation *interpolation, int allowOther );
static int pointwiseXY_C_getInterpolationFromTupleOfTwoStrings( PyObject *twoInterpolations, ptwXY_interpolation *interpolation );
static enum e_interpolationType pointwiseXY_C_getInterpolationTypeOfObject( PyObject *o );
static PyObject *pointwiseXY_C_GetNone( void );
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( const char *s, ... );
static int pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... );

DL_EXPORT( void ) initpointwiseXY_C( void );
/*
******************** pointwiseXY_CNewInitialize ************************
*/
static pointwiseXY_CPy *pointwiseXY_CNewInitialize( int infill, int safeDivide ) {

    pointwiseXY_CPy *self = (pointwiseXY_CPy *) PyObject_New( pointwiseXY_CPy, &pointwiseXY_CPyType );

    if( self ) {
        self->ptwXY = NULL;
        self->infill = infill;
        self->safeDivide = safeDivide;
    }
    return( self );
}
/*
************************ __init__ ************************************
*/
static int pointwiseXY_C__init__( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int infill = 1, safeDivide = 0, status = 0, userFlag = 0;
    int initialSize = 100, overflowSize = 10;
    double accuracy = _defaultAccuracy, biSectionMax = 3.;
    nfu_status status_nf;
    static char *kwlist[] = { "data", "dataForm", "initialSize", "overflowSize", "accuracy", "biSectionMax", "interpolation", "infill", 
        "safeDivide", "userFlag", NULL };
    ptwXYPoints *ptwXY;
    PyObject *dataPy = NULL, *dataFormPy = NULL, *xsPy = NULL, *ysPy = NULL, *theEnd = NULL, *iterator;
    char *interpolationStr = NULL, *dataForm, dataFormXYs[] = "xys", dataFormXsAndYs[] = "xsandys", dataFormList[] = "list", dataFormToLower[12], *c;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;

    self->ptwXY = NULL;
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|OOiiddsiii", kwlist, &dataPy, &dataFormPy, &initialSize, &overflowSize, &accuracy, 
        &biSectionMax, &interpolationStr, &infill, &safeDivide, &userFlag ) ) return( -1 );
    self->infill = infill;
    self->safeDivide = safeDivide;

    if( dataFormPy == NULL ) {
        dataForm = dataFormXYs; }
    else {
        if( !PyString_Check( dataFormPy ) ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "dataForm must be a string" ) );
        dataForm = PyString_AsString( dataFormPy );
        if( strlen( dataForm ) > strlen( dataFormXsAndYs ) ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid dataForm = '%s'", dataForm ) );
        strcpy( dataFormToLower, dataForm );
        for( c = dataFormToLower; *c != 0; c++ ) *c = (unsigned char ) tolower( *c );
        if( strcmp( dataFormToLower, dataFormXYs ) == 0 ) {
            dataForm = dataFormXYs; }
        else if( strcmp( dataFormToLower, dataFormXsAndYs ) == 0 ) {
            dataForm = dataFormXsAndYs ; }
        else if( strcmp( dataFormToLower, dataFormList ) == 0 ) {
            dataForm = dataFormList; }
        else {
            return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid dataForm = '%s'", dataForm ) );
        }
    }

    if( pointwiseXY_C_getInterpolationFromString( interpolationStr, &interpolation, 1 ) != 0 ) return( -1 );

    if( ( ptwXY = ptwXY_new( interpolation, biSectionMax, accuracy, initialSize, overflowSize, &status_nf, userFlag ) ) == NULL ) {
        PyErr_NoMemory( );
        return( -1 );
    }
    if( status_nf != nfu_Okay ) {
        free( ptwXY );
        PyErr_NoMemory( );
        return( -1 );
    }
    self->ptwXY = ptwXY;
    if( dataPy != NULL ) {
        if( dataForm == dataFormXYs ) {
            if( is_pointwiseXY_CPyObject( dataPy ) ) {
                if( pointwiseXY_C_setDataFromPtwXY( self, (pointwiseXY_CPy *) dataPy ) != 0 ) return( -1 ); }
            else {
                if( pointwiseXY_C_setData2( self, dataPy ) != 0 ) return( -1 );
            } }
        else if( dataForm == dataFormXsAndYs ) {
            if( ( iterator = PyObject_GetIter( dataPy ) ) == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs'" \
                ", data must be a list of length 2" ) );
            xsPy = PyIter_Next( iterator );
            if( xsPy != NULL ) ysPy = PyIter_Next( iterator );
            if( ysPy != NULL ) theEnd = PyIter_Next( iterator );
            if( theEnd != NULL ) {
                status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs', data must be a list length 2, it is longer" ); }
            else {
                if( ysPy == NULL ) {
                    status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs', data must be a list length 2, not 1" ); }
                else {
                    status = pointwiseXY_C_setDataFromXsAndYs2( self, xsPy, ysPy );
                }
            }
            Py_DECREF( iterator );
            if( xsPy != NULL ) { Py_DECREF( xsPy ); }
            if( ysPy != NULL ) { Py_DECREF( ysPy ); }
            if( theEnd != NULL ) { Py_DECREF( theEnd ); } }
        else {
            status = pointwiseXY_C_setDataFromList2( self, dataPy );
        }
    }
    return( status );
}
/*
************************************************************
*/
static void pointwiseXY_C_dealloc( PyObject *self ) {

    pointwiseXY_CPy *ptwXY = (pointwiseXY_CPy *) self;

    if( ptwXY->ptwXY != NULL ) {
        ptwXY_free( ptwXY->ptwXY );
    }
    self->ob_type->tp_free( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__repr__( pointwiseXY_CPy *self ) {

    ptwXYPoints *ptwXY = self->ptwXY;

    if( ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull("pointwiseXY_C object had a prior error and is not usable") );
    return( pointwiseXY_C_toString2( self, 1, " %16.8e %16.8e", "" ) );
}
/*
************************************************************
*/
static Py_ssize_t pointwiseXY_C__len__( pointwiseXY_CPy *self ) {

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne("pointwiseXY_C object had a prior error and is not usable") );
    return( (Py_ssize_t) ptwXY_length( self->ptwXY ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__getitem__( pointwiseXY_CPy *self, Py_ssize_t index_ ) {

    int64_t index = (int64_t) index_;
    ptwXYPoint *point;

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pointwiseXY_C object had a prior error and is not usable" ) );
    if( ( point = ptwXY_getPointAtIndex( self->ptwXY, index ) ) == NULL ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( NULL );
    }
    return( (PyObject *) Py_BuildValue( "[d,d]", point->x, point->y ) );
}
/*
************************************************************
*/
static void pointwiseXY_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 ) {

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
static PyObject *pointwiseXY_C__getslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_ ) {

    int64_t index1 = (int64_t) index1_, index2 = (int64_t) index2_, overflowSize = self->ptwXY->overflowAllocatedSize;
    pointwiseXY_CPy *newPy;
    ptwXYPoints *n;
    nfu_status status_nf;

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pointwiseXY_C object had a prior error and is not usable" ) );

    pointwiseXY_C_getSliceIndices( self->ptwXY->length, &index1, &index2 );
    if( ( newPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) return( NULL );
    if( overflowSize > 10 ) overflowSize = 10;
    if( ( n = ptwXY_slice( self->ptwXY, index1, index2, overflowSize, &status_nf ) ) == NULL ) {
        Py_DECREF( newPy );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from __getslice__: %s", nfu_statusMessage( status_nf ) ) );
    }
    newPy->ptwXY = n;
    return( (PyObject *) newPy );
}
/*
************************************************************
*/
static int pointwiseXY_C__setitem__( pointwiseXY_CPy *self, Py_ssize_t index_, PyObject *value ) {
/*
*   This routine can be called either when "self[index_] = number" or "del self[index_]".
*/
    int status = 0;
    int64_t index = (int64_t) index_;
    ptwXYPoint *point;
    pointwiseXY_CPy *other = (pointwiseXY_CPy *) value;
    double x, y;

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne("pointwiseXY_C object had a prior error and is not usable") );
    if( ( index < 0 ) || ( index >= self->ptwXY->length ) ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( -1 );
    }
    if( value == NULL ) {
        ptwXY_deletePoints( self->ptwXY, index, index + 1 ); }    /* Have already check for possible errors above, no need to here. */
    else {
        switch( pyObject_NumberOrPtwXY( value ) ) {
        case pyObject_Unsupported :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid object for set method" );
            break;
        case pyObject_Number :
            point = ptwXY_getPointAtIndex( self->ptwXY, index );
            pointwiseXY_C_PyNumberToFloat( value, &y );
            point->y = y;
            break;
        case pyObject_ptwXY :
            if( other->ptwXY->length != 1 ) {
                status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid length = %ld for set object; must be 1", other->ptwXY->length ); }
            else {
                point = ptwXY_getPointAtIndex( other->ptwXY, 0 );
                if( ptwXY_setXYPairAtIndex( self->ptwXY, index, point->x, point->y ) != nfu_Okay )
                    status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "x value of set object not between index = %ld's neighbors' x values", index );
            }
            break;
        default :       /* Some kind of sequence object. */
            status = pointwiseXY_C_PythonXYPairToCPair( value, &x, &y, 0 );
            if( status == 0 ) {
                 if( ptwXY_setXYPairAtIndex( self->ptwXY, index, x, y ) != nfu_Okay ) {
                    status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "x = %e value of set object not between index = %ld neighbors' x values", 
                        x, index ); } }
            else {
                status = -1;
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C__setslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *value ) {
/*
*   This routine can be called either when "self[index1_:index2_] = number" or "del self[index1_:index2_]".
*/
    int status = 0;
    int64_t i, index1 = (int64_t) index1_, index2 = (int64_t) index2_, length;
    double *xys, x, xMin, xMax;
    nfu_status status_nf;
    ptwXYPoints *n = NULL;

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "pointwiseXY_C object had a prior error and is not usable" ) );

    pointwiseXY_C_getSliceIndices( self->ptwXY->length, &index1, &index2 );

    if( value == NULL ) {
        ptwXY_deletePoints( self->ptwXY, index1, index2 ); }    /* Have already check for possible errors above, no need to here. */
    else {
        if( ( length = pointwiseXY_C_PythonXYListToCList( value, &xys ) ) < (int64_t) 0 ) return( -1 );
        if( ( n = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) 
            return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from __setslice__: %s", nfu_statusMessage( status_nf ) ) );
        if( ( self->ptwXY->length > 0 ) && ( length > 0 ) ) {
            x = xMin = *xys;
            xMax = xys[2 * ( length - 1 )];
            for( i = 1; i < length; i++ ) { 
                if( xys[2 * i] <= x ) {
                    ptwXY_free( n );
                    free( xys );
                    return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from __setslice__: x[i] = %g >= x[i+1] = %g", x, xys[2 * i] ) );
                }
                x = xys[2 * i];
            }
            ptwXY_deletePoints( n, index1, index2 );            /* Have already check for possible errors above, no need to here. */
            if( index1 > 0 ) {                                  /* Logic here requires that ptwXY_deletePoints coalesced points. */
                i = index1 - 1;
                if( index1 >= n->length ) i = n->length - 1;
                if( xMin <= n->points[i].x ) {
                    x = n->points[i].x;
                    ptwXY_free( n );
                    free( xys );
                    return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from __setslice__:xMin: x[i] = %g >= x[i+1] = %g", x, xMin ) );
                }
            }
            if( index1 < n->length ) {
                if( xMax >= n->points[index1].x ) {
                    x = n->points[index1].x;
                    ptwXY_free( n );
                    free( xys );
                    return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from __setslice__:xMax: x[i] = %g <= x[i+1] = %g, i = %d", 
                            xMax, x, index1 ) );
                }
            }
        }
        for( i = 0; i < length; i++ ) {
            if( ( status_nf = ptwXY_setValueAtX( n, xys[2 * i], xys[2 * i + 1] ) ) != nfu_Okay ) {
                ptwXY_free( n );
                free( xys );
                return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from __setslice__: %s", nfu_statusMessage( status_nf ) ) );
            }
        }
        ptwXY_free( self->ptwXY );
        self->ptwXY = n;
    }
    return( status );
}
/*
************************************************************
*/
static PySequenceMethods pointwiseXY_CPy_sequence = {
    (lenfunc) pointwiseXY_C__len__,               /* sq_length */
    0,                                          /* sq_concat */
    0,                                          /* sq_repeat */
    (ssizeargfunc) pointwiseXY_C__getitem__,        /* sq_item */
    (ssizessizeargfunc) pointwiseXY_C__getslice__,    /* sq_slice */
    (ssizeobjargproc) pointwiseXY_C__setitem__,     /* sq_ass_item */
    (ssizessizeobjargproc) pointwiseXY_C__setslice__, /* sq_ass_slice */
    0,                                          /* sq_contains */
    0,                                          /* sq_inplace_concat */
    0                                           /* sq_inplace_repeat */
};
/*
************************************************************
*/
static void pointwiseXY_C_Get_pointwiseXY_CAsSelf( PyObject *self, PyObject *other, PyObject **s, PyObject **o ) {

    *s = self;
    *o = other;
    if( !PyObject_TypeCheck( self, &pointwiseXY_CPyType ) ) {
        *s = other;
        *o = self;
    }
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__add__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s, *o;
    nfu_status status_nf;

    pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s, &o );
    switch( pyObject_NumberOrPtwXY( o ) ) {
    case pyObject_Number :
        n = pointwiseXY_C_add_sub_mul_div_number( (pointwiseXY_CPy *) s, o, ptwXY_add_double, "__add__", &status_nf  );
        break;
    case pyObject_ptwXY :
        n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( (pointwiseXY_CPy *) s, (pointwiseXY_CPy *) o, ptwXY_add_ptwXY, "__add__", &status_nf );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__iadd__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_add_double, "__iadd__" );
        break;
    case pyObject_ptwXY :
        self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_add_ptwXY, "__iadd__" );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__sub__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s, *o;
    nfu_status status_nf;

    pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s, &o );
    switch( pyObject_NumberOrPtwXY( o ) ) {
    case pyObject_Number :
        if( s == self ) {
            n = pointwiseXY_C_add_sub_mul_div_number( (pointwiseXY_CPy *) s, o, ptwXY_sub_doubleFrom, "__sub__", &status_nf ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_number( (pointwiseXY_CPy *) s, o, ptwXY_sub_fromDouble, "__rsub__", &status_nf );
        }
        break;
    case pyObject_ptwXY :
        n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( (pointwiseXY_CPy *) s, (pointwiseXY_CPy *) o, ptwXY_sub_ptwXY, "__sub__", &status_nf );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__isub__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_sub_doubleFrom, "__isub__" );
        break;
    case pyObject_ptwXY :
        self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_sub_ptwXY, "__isub__" );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__mul__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *sPy, *oPy;
    pointwiseXY_CPy *s, *o;
    nfu_status status_nf;

    pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &sPy, &oPy );
    s = (pointwiseXY_CPy *) sPy;
    switch( pyObject_NumberOrPtwXY( oPy ) ) {
    case pyObject_Number :
        n = pointwiseXY_C_add_sub_mul_div_number( s, oPy, ptwXY_mul_double, "__mul__", &status_nf );
        break;
    case pyObject_ptwXY :
        o = (pointwiseXY_CPy *) oPy;
        if( s->infill || o->infill ) {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s, o, ptwXY_mul2_ptwXY, "__mul__", &status_nf ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s, o, ptwXY_mul_ptwXY, "__mul__", &status_nf );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__imul__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_mul_double, "__mul__" );
        break;
    case pyObject_ptwXY :
        if( s1->infill || o1->infill ) {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_mul2_ptwXY, "__imul__" ); }
        else {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_mul_ptwXY, "__imul__" );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__div__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s, *o;
    nfu_status status_nf = nfu_Okay;
    ptwXYPoints *numerator, *ptwXYresult;
    pointwiseXY_CPy *sXY;

    pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s, &o );
    sXY = (pointwiseXY_CPy *) s;
    switch( pyObject_NumberOrPtwXY( o ) ) {
    case pyObject_Number :
        if( s == self ) {
            n = pointwiseXY_C_add_sub_mul_div_number( sXY, o, ptwXY_div_doubleFrom, "__div__", &status_nf ); }
        else {
            double number = 0, xMin = ptwXY_getXMin( sXY->ptwXY ), xMax = ptwXY_getXMax( sXY->ptwXY );

            pointwiseXY_C_PyNumberToFloat( o, &number );
            if( ( numerator = ptwXY_valueTo_ptwXY( ptwXY_interpolationLinLin, xMin, xMax, number, &status_nf ) ) == NULL ) {
                PyErr_NoMemory( );
                return( NULL );
            }
            if( sXY->safeDivide ) {
                ptwXYresult = pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( numerator, sXY->ptwXY, &status_nf ); }
            else {
                ptwXYresult = pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( numerator, sXY->ptwXY, &status_nf );
            }
            ptwXY_free( numerator );
            if( ptwXYresult == NULL ) {
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from __div__: %s", nfu_statusMessage( status_nf ) );
                n = NULL; }
            else {
                if( ( n = (PyObject *) pointwiseXY_CNewInitialize( ((pointwiseXY_CPy *) s)->infill, ((pointwiseXY_CPy *) s)->safeDivide ) ) == NULL ) {
                    ptwXY_free( ptwXYresult );
                    return( NULL );
                }
                ((pointwiseXY_CPy *) n)->ptwXY = ptwXYresult;
            }
        }
        break;
    case pyObject_ptwXY :
        if( ((pointwiseXY_CPy *) s)->safeDivide || ((pointwiseXY_CPy *) s)->safeDivide ) {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( (pointwiseXY_CPy *) s, (pointwiseXY_CPy *) o, 
                pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide, "__div__", &status_nf ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( (pointwiseXY_CPy *) s, (pointwiseXY_CPy *) o, 
                pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide, "__div__", &status_nf );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    if( status_nf == nfu_divByZero ) PyErr_SetString( PyExc_ZeroDivisionError, "for pointwiseXY instance" );
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__idiv__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_div_doubleFrom, "__idiv__" );
        break;
    case pyObject_ptwXY :
        if( s1->safeDivide || s1->safeDivide ) {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide, "__idiv__" ); }
        else {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide, "__idiv__" );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status_nf ) {

    return( ptwXY_div_ptwXY( ptwXY1, ptwXY2, status_nf, 0 ) );
}
/*
************************************************************
*/
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status_nf ) {

    return( ptwXY_div_ptwXY( ptwXY1, ptwXY2, status_nf, 1 ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__mod__( PyObject *self, PyObject *other ) {

    nfu_status status_nf;
    double m;
    PyObject *n, *o;

    if( ( o = PyNumber_Float( other ) ) == NULL ) return( NULL );
    m = PyFloat_AsDouble( o );
    if( ( n = pointwiseXY_C_copy( (pointwiseXY_CPy *) self ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwXY_mod( ((pointwiseXY_CPy *) n)->ptwXY, m, 1 ) ) ) {
        Py_DECREF( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from mod: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_number( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func, const char *oper, 
        nfu_status *status_nf ) {

    ptwXYPoints *n = NULL;
    pointwiseXY_CPy *nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide );
    double number = 0.;

    if( nPy != NULL ) {
        pointwiseXY_C_PyNumberToFloat( other, &number );       /* Assume calling routine has checked that other is a number. */
        n = ptwXY_clone( self->ptwXY, status_nf );
        if( n == NULL ) {
            pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( *status_nf ) ); }
        else {
            if( ( *status_nf = func( n, number ) ) == nfu_Okay ) {
                nPy->ptwXY = n; }
            else {
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( *status_nf ) );
                n = ptwXY_free( n );
            }
        }
    }
    if( ( n == NULL ) && ( nPy != NULL ) ) {
        Py_DECREF( nPy );
        nPy = NULL;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_number_insitu( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func, char const *oper ) {

    ptwXYPoints *n1 = NULL;
    double number = 0.;
    nfu_status status_nf;

    pointwiseXY_C_PyNumberToFloat( other, &number );       /* Assume calling routine has checked that other is a number. */
    n1 = ptwXY_clone( self->ptwXY, &status_nf );
    if( n1 == NULL ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( status_nf ) ) ); }
    else {
        if( ( status_nf = func( n1, number ) ) == nfu_Okay ) {
            self->ptwXY = n1; }
        else {
            ptwXY_free( n1 );
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( status_nf ) ) );
        }
    }
    return( (PyObject *) self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( pointwiseXY_CPy *self, pointwiseXY_CPy *other, ptwXY_ptwXY_ptwXY_func func, const char *oper,
        nfu_status *status_nf ) {

    ptwXYPoints *n = NULL;
    pointwiseXY_CPy *nPy;

    nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide );
    if( nPy != NULL ) {
        if( ( n = func( self->ptwXY, other->ptwXY, status_nf ) )  == NULL ) {
            if( *status_nf == nfu_mallocError ) {
                PyErr_NoMemory( ); }
            else {
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( *status_nf ) );
            } }
        else {
            nPy->ptwXY = n;
        }
    }
    if( ( n == NULL ) && ( nPy != NULL ) ) {
        Py_DECREF( nPy );
        nPy = NULL;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( pointwiseXY_CPy *self, pointwiseXY_CPy *other, ptwXY_ptwXY_ptwXY_func func, 
        char const *oper ) {

    ptwXYPoints *n1 = NULL;
    nfu_status status_nf;

    if( ( n1 = func( self->ptwXY, other->ptwXY, &status_nf ) )  == NULL ) {
        self = NULL;
        if( status_nf == nfu_mallocError ) {
            PyErr_NoMemory( ); }
        else {
            pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from %s: %s", oper, nfu_statusMessage( status_nf ) );
        } }
    else {
        ptwXY_free( self->ptwXY );
        self->ptwXY = n1;
    }
    return( (PyObject *) self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__pow__( PyObject *self, PyObject *other, PyObject *dummy ) {

    nfu_status status_nf;
    double p;
    PyObject *o, *n;

    if( dummy != Py_None) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pow() 3rd argument not allowed unless all arguments are integers" ) );
    if( ( o = PyNumber_Float( other ) ) == NULL ) return( NULL );
    p = PyFloat_AsDouble( o );
    if( ( n = pointwiseXY_C_copy( (pointwiseXY_CPy *) self ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwXY_pow( ((pointwiseXY_CPy *) n)->ptwXY, p ) ) ) {
        Py_DECREF( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from pow: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__neg__( PyObject *self ) {

    nfu_status status_nf;
    PyObject *n;

    if( ( n = pointwiseXY_C_copy( (pointwiseXY_CPy *) self ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwXY_neg( ((pointwiseXY_CPy *) n)->ptwXY ) ) != nfu_Okay ) {
        Py_DECREF( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from __neg__: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__abs__( PyObject *self ) {

    nfu_status status_nf;
    PyObject *n;

    if( ( n = pointwiseXY_C_copy( (pointwiseXY_CPy *) self ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwXY_abs( ((pointwiseXY_CPy *) n)->ptwXY ) ) != nfu_Okay ) {
        Py_DECREF( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from __abs__: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( n );
}
/*
************************************************************
*/
static PyNumberMethods pointwiseXY_CPy_number = {
    (binaryfunc) pointwiseXY_C__add__,      /* nb_add */
    (binaryfunc) pointwiseXY_C__sub__,      /* nb_subtract */
    (binaryfunc) pointwiseXY_C__mul__,      /* nb_multiply */
    (binaryfunc) pointwiseXY_C__div__,      /* nb_divide */
    (binaryfunc) pointwiseXY_C__mod__,      /* nb_remainder */
    0,                                      /* binaryfunc nb_divmod */
    (ternaryfunc) pointwiseXY_C__pow__,     /* nb_power */
    (unaryfunc) pointwiseXY_C__neg__,       /* nb_negative */
    0,                                      /* unaryfunc nb_positive */
    (unaryfunc) pointwiseXY_C__abs__,       /* nb_absolute */
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
    (binaryfunc) pointwiseXY_C__iadd__,     /* nb_inplace_add */            /* Added in release 2.0 */
    (binaryfunc) pointwiseXY_C__isub__,     /* nb_inplace_subtract */
    (binaryfunc) pointwiseXY_C__imul__,     /* nb_inplace_multiply */
    (binaryfunc) pointwiseXY_C__idiv__,     /* nb_inplace_divide */
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
static PyObject *pointwiseXY_C_allocatedSize( pointwiseXY_CPy *self ) {

    return( (PyObject *) Py_BuildValue( "l", self->ptwXY->allocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_applyFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    PyObject *f_args[2];
    PyObject *f, *parameters;
    int biSectionMax = -1, checkForRoots = 0;
    double accuracy = -1.;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "f", "parameters", "accuracy", "biSectionMax", "checkForRoots", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OO|dii", kwlist, &f, &parameters, &accuracy, &biSectionMax, &checkForRoots ) ) return( NULL );
    f_args[0] = f;
    f_args[1] = parameters;

    if( biSectionMax < 0 ) biSectionMax = ptwXY_getBiSectionMax( self->ptwXY );
    if( accuracy < 0 ) accuracy = ptwXY_getAccuracy( self->ptwXY );

    if( !PyCallable_Check( f ) ) {
        PyErr_SetString( PyExc_TypeError, "First argument must be callable" );
        return( NULL );
    }

    if( ( n = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from applyFunction: %s", nfu_statusMessage( status_nf ) ) );
    ptwXY_setBiSectionMax( n, biSectionMax );
    ptwXY_setAccuracy( n, accuracy );

    Py_INCREF( f );
    status_nf = ptwXY_applyFunction( n, pointwiseXY_C_applyFunction2, (void *) f_args, checkForRoots );
    Py_DECREF( f );

    if( status_nf != nfu_Okay ) {
        ptwXY_free( n );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static nfu_status pointwiseXY_C_applyFunction2( ptwXYPoint *ptwXY, void *argList ) {

    PyObject **f_args = (PyObject **) argList, *result;
    nfu_status status_nf = nfu_Okay;

    result = PyEval_CallFunction( (PyObject *) f_args[0], "(d,O)", ptwXY->y, f_args[1] );
    if( result == NULL ) return( nfu_badInput );
    if( pointwiseXY_C_PyNumberToFloat( result, &(ptwXY->y) ) != 0 ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert returned value to float" );
        status_nf = nfu_badInput;
    }
    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    char *interpolationStr = NULL;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;
    double accuracy = self->ptwXY->accuracy, lowerEps = 0., upperEps = 0.;
    static char *kwlist[] = { "interpolation", "accuracy", "lowerEps", "upperEps", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|sddd", kwlist, &interpolationStr, &accuracy, &lowerEps, &upperEps ) ) return( NULL );

    if( pointwiseXY_C_getInterpolationFromString( interpolationStr, &interpolation, 0 ) != 0 ) return( NULL );
    return( pointwiseXY_C_changeInterpolation2( self, interpolation, accuracy, lowerEps, upperEps ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolationIfNeeded( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int status, length = 0;
    double accuracy = self->ptwXY->accuracy, lowerEps = 0., upperEps = 0.;
    static char *kwlist[] = { "interpolation", "accuracy", "lowerEps", "upperEps", NULL };
    PyObject *allowedInterpolations, *iterator, *interpolationItem = NULL;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin, firstInterpolation = ptwXY_interpolationLinLin;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|ddd", kwlist, &allowedInterpolations, &accuracy, &lowerEps, &upperEps ) ) return( NULL );

    if( ( iterator = PyObject_GetIter( allowedInterpolations ) ) == NULL ) return( NULL );
    for( interpolationItem = PyIter_Next( iterator ); interpolationItem != NULL; interpolationItem = PyIter_Next( iterator ) ) {
        status = pointwiseXY_C_getInterpolationFromTupleOfTwoStrings( interpolationItem, &interpolation );
        Py_DECREF( interpolationItem );
        if( status < 0 ) {
            Py_DECREF( iterator );
            return( NULL );
        }
        if( interpolation == self->ptwXY->interpolation ) {
            Py_DECREF( iterator );
            Py_INCREF( self );
            return( (PyObject *) self );
        }
        if( length == 0 ) firstInterpolation = interpolation;
        ++length;
    }
    Py_DECREF( iterator );
    if( length == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Allowed interpolation sequence is empty" ) );
    return( pointwiseXY_C_changeInterpolation2( self, firstInterpolation, accuracy, lowerEps, upperEps ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolation2( pointwiseXY_CPy *self, ptwXY_interpolation interpolation, double accuracy, double lowerEps, double upperEps ) {

    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( self->ptwXY->interpolation == ptwXY_interpolationFlat ) {
        if( interpolation != ptwXY_interpolationLinLin )
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: can only change 'flat' to 'linear,linear'" ) );
        if( ( n = ptwXY_flatInterpolationToLinear( self->ptwXY, lowerEps, upperEps, &status_nf ) ) == NULL ) {
            if( ( lowerEps == 0 ) && ( upperEps == 0 ) )
                return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: both lowerEps and upperEps are 0." ) );
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: %s", nfu_statusMessage( status_nf ) ) );
        } }
    else {
        if( ( n = ptwXY_toOtherInterpolation( self->ptwXY, interpolation, accuracy, &status_nf ) ) == NULL )
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_clip( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    double yMin, yMax, yMinP = 1e99, yMaxP = -1e99;
    static char *kwlist[] = { "yMin", "yMax", NULL };
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &yMinP, &yMaxP ) ) return( NULL );
    n = ptwXY_clone( self->ptwXY, &status_nf );
    if( n == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from clip: %s", nfu_statusMessage( status_nf ) ) );
    if( self->ptwXY->length != 0 ) {
        yMin = ptwXY_getYMin( self->ptwXY );
        yMax = ptwXY_getYMax( self->ptwXY );
        if( yMinP != 1e99 ) yMin = yMinP;
        if( yMaxP != -1e99 ) yMax = yMaxP;
        status_nf = ptwXY_clip( n, yMin, yMax );
        if( status_nf != nfu_Okay ) {
            ptwXY_free( n );
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from clip: %s", nfu_statusMessage( status_nf ) ) );
        }
    }

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_coalescePoints( pointwiseXY_CPy *self ) {

    nfu_status status_nf;

    if( ( status_nf = ptwXY_coalescePoints( self->ptwXY, 0, NULL, 0 ) ) != nfu_Okay ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from coalescePoints: %s", nfu_statusMessage( status_nf ) ) );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_convolute( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int mode = 0;
    nfu_status status_nf;
    ptwXYPoints *n;
    PyObject *otherPy;
    pointwiseXY_CPy *nPy, *other;
    static char *kwlist[] = { "other", "mode", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|i", kwlist, &otherPy, &mode ) ) return( NULL );
    other = (pointwiseXY_CPy *) otherPy;

    if( ( n = ptwXY_convolution( self->ptwXY, ((pointwiseXY_CPy *) otherPy)->ptwXY, &status_nf, mode ) ) == NULL ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from convolute: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copy( pointwiseXY_CPy *self ) {

    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    n = ptwXY_clone( self->ptwXY, &status_nf );
    if( n == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from copy: %s", nfu_statusMessage( status_nf ) ) );
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copyDataToXYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int64_t i;
    PyObject *l = PyList_New( 0 ), *ni;
    nfu_status status_nf;
    double xScale = 1.0, yScale = 1.0;
    static char *kwlist[] = { "xScale", "yScale", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xScale, &yScale ) ) return( NULL );
    if( l == NULL ) return( NULL );
    if( ( status_nf = ptwXY_coalescePoints( self->ptwXY, 0, NULL, 0 ) ) != nfu_Okay )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from copyDataToXYs: %s", nfu_statusMessage( status_nf ) ) );
    for( i = 0; i < self->ptwXY->length; i++ ) {
        ni = Py_BuildValue( "[d,d]", xScale * self->ptwXY->points[i].x, yScale * self->ptwXY->points[i].y );
        if( pointwiseXY_C_addedItemToPythonList( l, ni ) != 0 ) return( NULL );
    }
    return( l );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copyDataToXsAndYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int64_t i;
    PyObject *lx = PyList_New( 0 ), *ly, *item;
    nfu_status status_nf;
    double xScale = 1.0, yScale = 1.0;
    static char *kwlist[] = { "xScale", "yScale", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xScale, &yScale ) ) return( NULL );
    if( lx == NULL ) return( NULL );
    if( ( status_nf = ptwXY_coalescePoints( self->ptwXY, 0, NULL, 0 ) ) != nfu_Okay )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from copyDataToXYs: %s", nfu_statusMessage( status_nf ) ) );

    if( ( ly = PyList_New( 0 ) ) == NULL ) {
        Py_DECREF( lx );
        return( NULL );
    }

    for( i = 0; i < self->ptwXY->length; i++ ) {
        item = Py_BuildValue( "d", xScale * self->ptwXY->points[i].x );
        if( pointwiseXY_C_addedItemToPythonList( lx, item ) != 0 ) {
            Py_DECREF( ly );
            return( NULL );
        }
        item = Py_BuildValue( "d", yScale * self->ptwXY->points[i].y );
        if( pointwiseXY_C_addedItemToPythonList( ly, item ) != 0 ) {
            Py_DECREF( lx );
            return( NULL );
        }
    }
    return( Py_BuildValue( "[O,O]", lx, ly ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_dullEdges( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "lowerEps", "upperEps", "positiveXOnly", NULL };
    int positiveXOnly = 0;
    double lowerEps = 0., upperEps = 0.;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ddi", kwlist, &lowerEps, &upperEps, &positiveXOnly ) ) return( NULL );

    if( ( n = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from dullEdges: %s", nfu_statusMessage( status_nf ) ) );
    if( ( status_nf = ptwXY_dullEdges( n, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay ) {
        ptwXY_free( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from dullEdges: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_exp( pointwiseXY_CPy *self, PyObject *args ) {

    double a;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTuple( args, "d", &a ) ) return( NULL );

    if( ( n = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from exp: %s", nfu_statusMessage( status_nf ) ) );
    if( ( status_nf = ptwXY_exp( n, a ) ) != nfu_Okay ) {
        ptwXY_free( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from exp: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getAccuracy( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "d", ptwXY_getAccuracy( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getBiSectionMax( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "d", (double) ptwXY_getBiSectionMax( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getInfill( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "i", self->infill ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getInterpolation( pointwiseXY_CPy *self ) {

    char *interpolationStr;

    if( self->ptwXY->interpolation == ptwXY_interpolationLinLin ) {
        interpolationStr = "linear,linear"; }
    else if( self->ptwXY->interpolation == ptwXY_interpolationLinLog ) {
        interpolationStr = "linear,log"; }
    else if( self->ptwXY->interpolation == ptwXY_interpolationLogLin ) {
        interpolationStr = "log,linear"; }
    else if( self->ptwXY->interpolation == ptwXY_interpolationLogLog ) {
        interpolationStr = "log,log"; }
    else if( self->ptwXY->interpolation == ptwXY_interpolationFlat ) {
        interpolationStr = "flat"; }
    else {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "INTERNAL ERROR: unknown interpolation = %d", self->ptwXY->interpolation ) );
    }
    return( Py_BuildValue( "s", interpolationStr ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getSecondaryCacheSize( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "l", (long int) self->ptwXY->overflowAllocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getSafeDivide( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "i", self->safeDivide ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getValue( pointwiseXY_CPy *self, PyObject *args ) {

    double x, y;
    nfu_status status_nf;
    PyObject *value = NULL;

    if( !PyArg_ParseTuple( args, "d", &x ) ) return( NULL );

    switch( status_nf = ptwXY_getValueAtX( self->ptwXY, x, &y ) ) {
    case nfu_Okay :
        value = Py_BuildValue( "d", y );
        break;
    case nfu_XOutsideDomain :
        value = pointwiseXY_C_GetNone( );
        break;
    case nfu_invalidInterpolation :
        value = pointwiseXY_C_SetPyErrorExceptionReturnNull( "unsupported interpolation = %d", self->ptwXY->interpolation );
        break;
    default :
        value = pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from getValue: %s", nfu_statusMessage( status_nf ) );
    }

    return( value );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getUserFlag( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "i", ptwXY_getUserFlag( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrate( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "xMin", "xMax", NULL };
    double xMin, xMax, v;
    nfu_status status_nf;

    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    xMin = ptwXY_getXMin( self->ptwXY );
    xMax = ptwXY_getXMax( self->ptwXY );
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xMin, &xMax ) ) return( NULL );

    v = ptwXY_integrate( self->ptwXY, xMin, xMax, &status_nf );
    if( status_nf != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from integrate: %s", nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", v ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrateWithWeight_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "xMin", "xMax", NULL };
    double xMin, xMax, v;
    nfu_status status_nf;

    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    xMin = ptwXY_getXMin( self->ptwXY );
    xMax = ptwXY_getXMax( self->ptwXY );
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xMin, &xMax ) ) return( NULL );

    v = ptwXY_integrateWithWeight_x( self->ptwXY, xMin, xMax, &status_nf );
    if( status_nf != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from integrateWithWeight_x: %s", nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", v ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrateWithWeight_sqrt_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "xMin", "xMax", NULL };
    double xMin, xMax, v;
    nfu_status status_nf;

    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    xMin = ptwXY_getXMin( self->ptwXY );
    xMax = ptwXY_getXMax( self->ptwXY );
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xMin, &xMax ) ) return( NULL );

    v = ptwXY_integrateWithWeight_sqrt_x( self->ptwXY, xMin, xMax, &status_nf );
    if( status_nf != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from integrateWithWeight_sqrt_x: %s", nfu_statusMessage( status_nf ) ) );
    return( Py_BuildValue( "d", v ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_normalize( pointwiseXY_CPy *self ) {

    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    n = ptwXY_clone( self->ptwXY, &status_nf );
    if( n == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from normalize: %s", nfu_statusMessage( status_nf ) ) );
    if( ( status_nf = ptwXY_normalize( n ) ) != nfu_Okay ) {
        ptwXY_free( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from normalize: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_runningIntegral( pointwiseXY_CPy *self ) {

    nfu_status status_nf;
    ptwXPoints *ptwX;

    if( ( ptwX = ptwXY_runningIntegral( self->ptwXY, &status_nf ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from ptwXY_runningIntegral: %s", nfu_statusMessage( status_nf ) ) );
    return( pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwX ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupOneFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL;
    static char *kwlist[] = { "groupBoundaries", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|O", kwlist, &groupBoundariesPy, &normPy ) ) return( NULL );

    return( pointwiseXY_C_groupFunctionsCommon( self, NULL, NULL, groupBoundariesPy, normPy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupTwoFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL, *f2;
    static char *kwlist[] = { "f2", "groupBoundaries", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OO|O", kwlist, &f2, &groupBoundariesPy, &normPy ) ) return( NULL );

    return( pointwiseXY_C_groupFunctionsCommon( self, f2, NULL, groupBoundariesPy, normPy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupThreeFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL, *f2, *f3;
    static char *kwlist[] = { "f2", "f3", "groupBoundaries", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OOO|O", kwlist, &f2, &f3, &groupBoundariesPy, &normPy ) ) return( NULL );

    return( pointwiseXY_C_groupFunctionsCommon( self, f2, f3, groupBoundariesPy, normPy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupFunctionsCommon( pointwiseXY_CPy *f1, PyObject *f2, PyObject *f3, PyObject *groupBoundariesPy, PyObject *normPy ) {

    ptwXPoints *ptwXGBs, *ptwX_norm = NULL, *groups = NULL;
    PyObject *newPy = NULL;
    ptwXY_group_normType norm = ptwXY_group_normType_none;
    nfu_status status_nf;
    char *normChars;
    int status;

    if( f2 != NULL ) {
        if( ( status = PyObject_IsInstance( f2, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
        if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "f2 must be a pointwiseXY_C instance" ) );
        if( f3 != NULL ) {
            if( ( status = PyObject_IsInstance( f3, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
            if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "f3 must be a pointwiseXY_C instance" ) );
        }
    }

    if( ( normPy != NULL ) && ( normPy != Py_None ) ) {
        if( PyString_Check( normPy ) ) {
            if( ( normChars = PyString_AsString( normPy ) ) == NULL ) {
                return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid string norm" ) ); }
            else {
                if( strcmp( "dx", normChars ) != 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid string norm identifier, must be 'dx'" ) );
                norm = ptwXY_group_normType_dx;
            } }
        else {
            if( ( ptwX_norm = pointwiseXY_C_PyFloatList_to_ptwXPoints( normPy ) ) == NULL ) return( NULL );
            norm = ptwXY_group_normType_norm;
        }
    }

    if( ( ptwXGBs = pointwiseXY_C_PyFloatList_to_ptwXPoints( groupBoundariesPy ) ) == NULL ) {
        if( ptwX_norm != NULL ) ptwX_free( ptwX_norm );
        return( NULL );
    }

    if( f2 == NULL ) {
        groups = ptwXY_groupOneFunction( f1->ptwXY, ptwXGBs, norm, ptwX_norm, &status_nf ); }
    else if( f3 == NULL ) {
        groups = ptwXY_groupTwoFunctions( f1->ptwXY, ((pointwiseXY_CPy *) f2)->ptwXY, ptwXGBs, norm, ptwX_norm, &status_nf ); }
    else {
        groups = ptwXY_groupThreeFunctions( f1->ptwXY, ((pointwiseXY_CPy *) f2)->ptwXY, ((pointwiseXY_CPy *) f3)->ptwXY, ptwXGBs, norm, ptwX_norm, &status_nf );
    }

    ptwX_free( ptwXGBs );
    if( ptwX_norm != NULL ) ptwX_free( ptwX_norm );

    if( groups == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "Grouping error: %s", nfu_statusMessage( status_nf ) ); }
    else {
        newPy = pointwiseXY_C_ptwXPoints_to_PyFloatList( groups );
        ptwX_free( groups );
    }

    return( newPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_areDomainsMutual( pointwiseXY_CPy *self, PyObject *args ) {

    nfu_status status_nf;
    PyObject *other, *rPy = NULL;
    int status;

    if( !PyArg_ParseTuple( args, "O", &other ) ) return( NULL );
    if( ( status = PyObject_IsInstance( (PyObject* ) other, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "other must be a pointwiseXY_C instance" ) );

    switch( status_nf = ptwXY_areDomainsMutual( self->ptwXY, ((pointwiseXY_CPy *) other)->ptwXY ) ) {
    case nfu_Okay :
    case nfu_empty :
        if( ( rPy = Py_BuildValue( "i", 1 ) ) == NULL ) return( NULL );
        break;
    case nfu_domainsNotMutual :
    case nfu_tooFewPoints :
        if( ( rPy = Py_BuildValue( "i", 0 ) ) == NULL ) return( NULL );
        break;
    default :
        rPy = pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from ptwXY_areDomainsMutual: %s", nfu_statusMessage( status_nf ) );
    }

    return( rPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_mutualify( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int positiveXOnly1, positiveXOnly2, status;
    double lowerEps1, upperEps1, lowerEps2, upperEps2;
    static char *kwlist[] = { "lowerEps1", "upperEps1", "positiveXOnly1", "other", "lowerEps2", "upperEps2", "positiveXOnly2", NULL };
    PyObject *other, *lPy;
    ptwXYPoints *n1 = NULL, *n2 = NULL;
    pointwiseXY_CPy *n1Py = NULL, *n2Py = NULL;
    nfu_status status_nf;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "ddiOddi", kwlist, &lowerEps1, &upperEps1, &positiveXOnly1, &other, &lowerEps2, 
        &upperEps2, &positiveXOnly2 ) ) return( NULL );

    if( ( status = PyObject_IsInstance( other, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "other must be a pointwiseXY_C instance" ) );

    if( ( n1 = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) goto Err;
    if( ( n2 = ptwXY_clone( ((pointwiseXY_CPy *) other)->ptwXY, &status_nf ) ) == NULL ) goto Err;

    if( ( status_nf = ptwXY_mutualifyDomains( n1, lowerEps1, upperEps1, positiveXOnly1, n2, lowerEps2, upperEps2, positiveXOnly2 ) ) != nfu_Okay ) goto Err;

    if( ( n1Py = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) goto Err;
    if( ( n2Py = pointwiseXY_CNewInitialize( ((pointwiseXY_CPy *) other)->infill, ((pointwiseXY_CPy *) other)->safeDivide ) ) == NULL ) goto Err;
    if( ( lPy = Py_BuildValue( "(O,O)", n1Py, n2Py ) ) == NULL ) goto Err;
    n1Py->ptwXY = n1;
    n2Py->ptwXY = n2;

    return( lPy );

Err:
    if( n1 != NULL ) ptwXY_free( n1 );
    if( n2 != NULL ) ptwXY_free( n2 );
    if( n1Py != NULL ) { Py_DECREF( n1Py ); }
    if( n2Py != NULL ) { Py_DECREF( n2Py ); }
    return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from mutualify: %s", nfu_statusMessage( status_nf ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_overflowAllocatedSize( pointwiseXY_CPy *self ) {

    return( (PyObject *) Py_BuildValue( "l", self->ptwXY->overflowAllocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_overflowLength( pointwiseXY_CPy *self ) {

    return( (PyObject *) Py_BuildValue( "l", self->ptwXY->overflowLength ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_plot( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *moduleName, *GnuplotModule, *Gnuplot = NULL, *Data = NULL, *g = NULL, *data = NULL, *status = NULL;

    if( ( moduleName = PyString_FromString( "Gnuplot" ) ) == NULL ) return( NULL );
    GnuplotModule = PyImport_Import( moduleName );
    Py_DECREF( moduleName );
    if( GnuplotModule == NULL ) return( NULL );

    if( ( Gnuplot = PyObject_GetAttrString( GnuplotModule, "Gnuplot" ) ) == NULL ) goto err;
    if( ( g = PyEval_CallFunction( Gnuplot, "()" ) ) == NULL ) goto err;
    if( PyEval_CallFunction( g, "(s)", "set style data linespoints" ) == NULL ) goto err;

    if( ( Data = PyObject_GetAttrString( GnuplotModule, "Data" ) ) == NULL ) goto err;

    if( ptwXY_length( self->ptwXY ) > 0 ) {
        if( ( data = PyEval_CallFunction( Data, "(O)", self ) ) == NULL ) goto err; }
    else {
        if( ( data = PyEval_CallFunction( Data, "([[d,d]])", -.9999, -.9999 ) ) == NULL ) goto err;
        PyEval_CallFunction( g, "(s)", "set xrange [-1:1]" );
        PyEval_CallFunction( g, "(s)", "set yrange [-1:1]" );
        PyEval_CallFunction( g, "(s)", "set label 'Instance has no data' center at 0,0" );
    }
    if( PyEval_CallMethod( g, "plot", "(O)", data ) == NULL ) goto err;

    status = g;
    goto theEnd;

err:
    status = NULL;

theEnd:
    Py_DECREF( GnuplotModule );
    if( Data != NULL ) { Py_DECREF( Data ); }
    if( data != NULL ) { Py_DECREF( data ); }
    if( Gnuplot != NULL ) { Py_DECREF( Gnuplot ); }
    if( ( status == NULL ) && ( g != NULL ) ) { Py_DECREF( g ); }
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_reallocatePoints( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    int size;
    int forceSmallerResize = 1;
    static char *kwlist[] = { "size", "forceSmaller", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "i|i", kwlist, &size, &forceSmallerResize ) ) return( NULL );

    if( ( status_nf = ptwXY_reallocatePoints( self->ptwXY, size, forceSmallerResize ) ) != nfu_Okay )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from reallocatePoints: %s", nfu_statusMessage( status_nf ) ) );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_reallocateOverflowPoints( pointwiseXY_CPy *self, PyObject *args ) {

    nfu_status status_nf;
    int size;

    if( !PyArg_ParseTuple( args, "i", &size ) ) return( NULL );

    if( ( status_nf = ptwXY_reallocateOverflowPoints( self->ptwXY, size ) ) != nfu_Okay )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from reallocateOverflowPoints: %s", nfu_statusMessage( status_nf ) ) );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setData( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXYList;

    if( !PyArg_ParseTuple( args, "O", &PyXYList ) ) return( NULL );
    if( pointwiseXY_C_setData2( self, PyXYList ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setData2( pointwiseXY_CPy *self, PyObject *PyXYList ) {

    int status = 0;
    int64_t length;
    double *xys;
    nfu_status status_nf;

    length = pointwiseXY_C_PythonXYListToCList( PyXYList, &xys );
    if( length == -1 ) return( -1 );
    if( ( status_nf = ptwXY_setXYData( self->ptwXY, length, xys ) ) != nfu_Okay ) {
        status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( nfu_statusMessage( status_nf ) );
    }
    free( xys );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromPtwXY( pointwiseXY_CPy *self, pointwiseXY_CPy *otherPY ) {

    nfu_status status_nf = ptwXY_copy( self->ptwXY, otherPY->ptwXY );

    if( status_nf != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( nfu_statusMessage( status_nf ) ) );
    return( 0 );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setAccuracy( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );

    return( Py_BuildValue( "d", ptwXY_setAccuracy( self->ptwXY, accuracy ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setBiSectionMax( pointwiseXY_CPy *self, PyObject *args ) {

    double biSectionMax;

    if( !PyArg_ParseTuple( args, "d", &biSectionMax ) ) return( NULL );

    return( Py_BuildValue( "d", ptwXY_setBiSectionMax( self->ptwXY, biSectionMax ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setSecondaryCacheSize( pointwiseXY_CPy *self, PyObject *args ) {

    long int length;
    nfu_status status_nf;

    if( !PyArg_ParseTuple( args, "l", &length ) ) return( NULL );

    if( ( status_nf = ptwXY_reallocateOverflowPoints( self->ptwXY, length ) ) != nfu_Okay ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( nfu_statusMessage( status_nf ) ) );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setDataFromList( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXYs;

    if( !PyArg_ParseTuple( args, "O", &PyXYs) ) return( NULL );
    if( pointwiseXY_C_setDataFromList2( self, PyXYs ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromList2( pointwiseXY_CPy *self, PyObject *list ) {

    int status = 0;
    nfu_status status_nf;
    int64_t i, j, length;
    double *xys = NULL, x, xPrior = 0.;

    self->ptwXY->length = 0;
    if( ( length = pointwiseXY_C_pythonDoubleListToCList( list, &xys, 0 ) ) == -1 ) return( -1 );
    if( 2 * ( length / 2 ) != length ) {
        status = -1;
        free( xys );
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "length = %d is not even", length ); }
    else {
        if( ( status_nf = ptwXY_reallocatePoints( self->ptwXY, length / 2, 1 ) ) == nfu_Okay ) {
            for( i = 0, j = 0; i < length; i += 2, j++ ) {
                x = xys[i];
                if( i > 0 ) {
                    if( x <= xPrior ) {
                        pointwiseXY_C_SetPyErrorExceptionReturnNull( "data not in ascending order at index %d and %d", (int) i - 1, (int) i );
                        status = -1;
                        break;
                    }
                }
                self->ptwXY->points[j].x = x;
                self->ptwXY->points[j].y = xys[i+1];
                xPrior = x;
            }
            self->ptwXY->length = length / 2; }
        else {
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from setValue: %s", nfu_statusMessage( status_nf ) );
        }
        free( xys );
    }
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setDataFromXsAndYs( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXs, *PyYs;

    if( !PyArg_ParseTuple( args, "OO", &PyXs, &PyYs) ) return( NULL );
    if( pointwiseXY_C_setDataFromXsAndYs2( self, PyXs, PyYs ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromXsAndYs2( pointwiseXY_CPy *self, PyObject *PyXs, PyObject *PyYs ) {

    int status = 0;
    nfu_status status_nf;
    int64_t i, xLength, yLength;
    double *xs = NULL, *ys = NULL;

    self->ptwXY->length = 0;
    if( ( xLength = pointwiseXY_C_pythonDoubleListToCList( PyXs, &xs, 1 ) ) == -1 ) return( -1 );
    if( ( yLength = pointwiseXY_C_pythonDoubleListToCList( PyYs, &ys, 0 ) ) == -1 ) {
        free( xs );
        return( -1 );
    }
    if( xLength != yLength ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "x and y lengths are not the same" );
        status = -1;
    }
    if( status == 0 ) {
        if( ( status_nf = ptwXY_reallocatePoints( self->ptwXY, xLength, 1 ) ) == nfu_Okay ) {
            for( i = 0; i < xLength; i++ ) {
                self->ptwXY->points[i].x = xs[i];
                self->ptwXY->points[i].y = ys[i];
            }
            self->ptwXY->length = xLength; }
        else {
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Error from setValue: %s", nfu_statusMessage( status_nf ) );
        }
    }
    if( xs != NULL ) free( xs );
    if( ys != NULL ) free( ys );
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setInfill( pointwiseXY_CPy *self, PyObject *args ) {

    if( !PyArg_ParseTuple( args, "i", &(self->infill) ) ) return( NULL );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setSafeDivide( pointwiseXY_CPy *self, PyObject *args ) {

    if( !PyArg_ParseTuple( args, "i", &(self->safeDivide) ) ) return( NULL );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setValue( pointwiseXY_CPy *self, PyObject *args ) {

    double x, y;
    nfu_status status_nf;

    if( !PyArg_ParseTuple( args, "dd", &x, &y ) ) return( NULL );

    status_nf = ptwXY_setValueAtX( self->ptwXY, x, y );
    if( status_nf != nfu_Okay ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from setValue: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setUserFlag( pointwiseXY_CPy *self, PyObject *args ) {

    int flag;

    if( !PyArg_ParseTuple( args, "i", &flag ) ) return( NULL );

    ptwXY_setUserFlag( self->ptwXY, flag );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_showInteralStructure( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int printPointersAsNull = 0;
    static char *kwlist[] = { "printPointersAsNull", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &printPointersAsNull ) ) return( NULL );

    ptwXY_showInteralStructure( self->ptwXY, stdout, printPointersAsNull );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_thicken( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int sectionSubdivideMax = 1;
    double dxMax = 0, fxMax = 1.;
    static char *kwlist[] = { "sectionSubdivideMax", "dxMax", "fxMax", NULL };
    nfu_status status_nf;
    PyObject *n;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|idd", kwlist, &sectionSubdivideMax, &dxMax, &fxMax ) ) return( NULL );
    if( ( n = pointwiseXY_C_copy( self ) ) == NULL ) return( NULL );
    if( ( status_nf = ptwXY_thicken( ((pointwiseXY_CPy *) n)->ptwXY, sectionSubdivideMax, dxMax, fxMax ) ) ) {
        Py_DECREF( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from thicken: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_thin( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );

    if( ( n = ptwXY_thin( self->ptwXY, accuracy, &status_nf ) ) == NULL ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from thin: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_toString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int pairsPerLine = 1;
    char *format = " %16.8e %16.8e", *pairSeparator = "";
    static char *kwlist[] = { "pairsPerLine", "format", "pairSeparator", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|iss", kwlist, &pairsPerLine, &format, &pairSeparator ) ) return( NULL );
    return( pointwiseXY_C_toString2( self, pairsPerLine, format, pairSeparator ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_toString2( pointwiseXY_CPy *self, int pairsPerLine, char *format, char *pairSeparator ) {

    int lineFeed, extraSpace = 4, dummyLen;     /* '\n' plus three more for reserve. */
    ptwXYPoints *ptwXY = self->ptwXY;
    int64_t index, strLength, pairSeparatorLength = strlen( pairSeparator );
    char *s, *e, *p, dummy[1024];
    ptwXYPoint *point;
    PyObject *str;

    if( ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pointwiseXY_C object had a prior error and is not usable" ) );

    s = pointwiseXY_C_toString_isFormatForDouble( format, 0 );
    if( ( s = pointwiseXY_C_toString_isFormatForDouble( s, 0 ) ) != NULL ) {
        s = pointwiseXY_C_toString_isFormatForDouble( s, 1 );
/*
        snprintf( dummy, 0, "'%s' %% ( -3.1238972718972345e-21, -6.1432987987397648974e+33 )", format );
        if( PyRun_SimpleString( dummy ) < 0 ) s = NULL;
*/
    }
    if( s == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid format = '%s' for converting two doubles to a string", format ) );

    dummyLen = snprintf( dummy, 0, format, -3.1238972718972345e-21, -6.1432987987397648974e-33 );
    if( dummyLen < 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Invalid format string = '%s'", format ) );
    strLength = ptwXY->length * ( dummyLen + pairSeparatorLength + extraSpace );

    if( ( s = (char *) malloc( (size_t) strLength * sizeof( char ) ) ) == NULL ) {
        PyErr_NoMemory( );
        return( NULL );
    }
    s[0] = 0;
    for( index = 0, e = s, lineFeed = pairsPerLine; index < ptwXY->length; index++ ) {
        point = ptwXY_getPointAtIndex( ptwXY, index );
        if( point == NULL ) {
            free( s );
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Internal bug: ptwXY_getPointAtIndex is NULL" ) );
        }
        sprintf( e, format, point->x, point->y );
        while( *e ) e++;
        if( index < ptwXY->length - 1 ) {
            for( p = pairSeparator; *p; p++, e++ ) *e = *p;
            *e = 0;
        }
        lineFeed--;
        if( lineFeed <= 0 ) {
            *e = '\n';
            e++;
            *e = 0;
            lineFeed = pairsPerLine;
        }
    }
    if( lineFeed != pairsPerLine ) {
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
static PyObject *pointwiseXY_C_trim( pointwiseXY_CPy *self ) {

    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    n = ptwXY_clone( self->ptwXY, &status_nf );
    if( n == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from trim: %s", nfu_statusMessage( status_nf ) ) );
    if( ( status_nf = ptwXY_trim( n ) ) != nfu_Okay ) {
        ptwXY_free( n );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from trim: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static char *pointwiseXY_C_toString_isFormatForDouble( char *format, int returnFormat ) {

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
static PyObject *pointwiseXY_C_union(  pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int fillWithSelf = 0, trim = 0, unionOptions = 0;
    static char *kwlist[] = { "other", "fillWithSelf", "trim", NULL };
    nfu_status status_nf;
    pointwiseXY_CPy *other, *nPy;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|ii", kwlist, &other, &fillWithSelf, &trim ) ) return( NULL );

    if( fillWithSelf ) unionOptions |= ptwXY_union_fill;
    if( trim ) unionOptions |= ptwXY_union_trim;

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide ) ) == NULL ) return( NULL );
    if( ( nPy->ptwXY = ptwXY_union( self->ptwXY, other->ptwXY, &status_nf, unionOptions ) ) == NULL ) {
        Py_DECREF( (PyObject *) nPy );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from union: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_mergeClosePoints(  pointwiseXY_CPy *self, PyObject *args ) {

    double epsilon;
    nfu_status status_nf;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTuple( args, "d", &epsilon ) ) return( NULL );
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) return( NULL );

    if( ( nPy->ptwXY = ptwXY_clone( self->ptwXY, &status_nf ) ) == NULL ) {
        Py_DECREF( (PyObject *) nPy );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from mergeClosePoints: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( ( status_nf = ptwXY_mergeClosePoints( nPy->ptwXY, epsilon ) ) != nfu_Okay ) {
        Py_DECREF( (PyObject *) nPy );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from mergeClosePoints: %s", nfu_statusMessage( status_nf ) ) );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getDomainGrid( pointwiseXY_CPy *self, PyObject *args ) {

    int i;
    double scale = 1.0;
    ptwXPoints *ptwX;
    nfu_status status_nf;
    PyObject *nPy;

    if( !PyArg_ParseTuple( args, "|d", &scale ) ) return( NULL );
    ptwX = ptwXY_getXArray( self->ptwXY, &status_nf );
    if( status_nf != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from getDomainGrid: %s", nfu_statusMessage( status_nf ) ) );
    for( i = 0; i < ptwX_length( ptwX ); i++ ) ptwX_setPointAtIndex( ptwX, i, scale * ptwX_getPointAtIndex_Unsafely( ptwX, i ) );
    nPy = pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwX );
    ptwX_free( ptwX );
    return( nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_xMin( pointwiseXY_CPy *self ) {

    double xMin;

    if( isOkayAndHasData( self->ptwXY ) == -1 ) return( NULL );
    xMin = ptwXY_getXMin( self->ptwXY );                          /* ???? Should check for nan or infty. */
    return( Py_BuildValue( "d", xMin ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_xMax( pointwiseXY_CPy *self ) {

    double xMax;

    if( isOkayAndHasData( self->ptwXY ) == -1 ) return( NULL );
    xMax = ptwXY_getXMax( self->ptwXY );                          /* ???? Should check for nan or infty. */
    return( Py_BuildValue( "d", xMax ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_xSlice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int fill = 1;
    double xMin = 1e99, xMax = -1e99, dullEps = 0;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "xMin", "xMax", "fill", "dullEps", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ddid", kwlist, &xMin, &xMax, &fill, &dullEps) ) return( NULL );

    if( ( xMin ==  1e99 ) && ( self->ptwXY->length > 0 ) ) xMin = ptwXY_getXMin( self->ptwXY );
    if( ( xMax == -1e99 ) && ( self->ptwXY->length > 0 ) ) xMax = ptwXY_getXMax( self->ptwXY );

    if( ( n = ptwXY_xSlice( self->ptwXY, xMin, xMax, 10, fill, &status_nf ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from xSlice: %s", nfu_statusMessage( status_nf ) ) );
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n ); }
    else {
        nPy->ptwXY = n;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_yMin( pointwiseXY_CPy *self ) {

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pointwiseXY_C object had a prior error and is not usable" ) );
    return( Py_BuildValue( "d", ptwXY_getYMin( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_yMax( pointwiseXY_CPy *self ) {

    if( self->ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pointwiseXY_C object had a prior error and is not usable" ) );
    return( Py_BuildValue( "d", ptwXY_getYMax( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_defaultAccuracy( pointwiseXY_CPy *self, PyObject *args ) {

    return( Py_BuildValue( "d", _defaultAccuracy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_createFromFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    PyObject *f_args[3], *xs_Py;
    PyObject *f, *parameters;
    int checkForRoots = 0, infill = 1, safeDivide = 1;
    double accuracy, biSectionMax;
    ptwXPoints *xs = NULL;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "xs", "f", "parameters", "accuracy", "biSectionMax", "checkForRoots", "infill", "safeDivide", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OOOdd|iii", kwlist, &xs_Py, &f, &parameters, &accuracy, &biSectionMax, &checkForRoots,
        &infill, &safeDivide ) ) return( NULL );
    f_args[0] = f;
    f_args[1] = parameters;
    f_args[2] = NULL;

    if( !PyCallable_Check( f ) ) {
        PyErr_SetString( PyExc_TypeError, "Second argument must be callable" );
        return( NULL );
    }

    if( ( xs = pointwiseXY_C_PyFloatList_to_ptwXPoints( xs_Py ) ) == NULL ) return( NULL );

    Py_INCREF( f );
    n = ptwXY_createFromFunction2( xs, pointwiseXY_C_createFromFunction2, (void *) f_args, accuracy, checkForRoots, (int) biSectionMax, &status_nf );
    Py_DECREF( f );

    ptwX_free( xs );
    if( n == NULL ) {
        if( f_args[2] != 0 ) return( NULL );
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from ptwXY_createFromFunction: %s", nfu_statusMessage( status_nf ) ) );
    }

    if( ( nPy = pointwiseXY_CNewInitialize( infill, safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static nfu_status pointwiseXY_C_createFromFunction2( double x, double *y, void *argList ) {

    PyObject **f_args = (PyObject **) argList, *result;
    nfu_status status_nf = nfu_Okay;

    result = PyEval_CallFunction( (PyObject *) f_args[0], "(d,O)", x, f_args[1] );
    if( result == NULL ) {
        f_args[2] = f_args[0];
        return( nfu_badInput );
    }
    if( pointwiseXY_C_PyNumberToFloat( result, y ) != 0 ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert returned value from createFromFunction function to float" );
        status_nf = nfu_badInput;
        f_args[2] = f_args[0];
    }
    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_createFromString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    int infill = 1, safeDivide = 1;
    double accuracy, biSectionMax;
    ptwXYPoints *ptwXY;
    pointwiseXY_CPy *ptwXYPy;
    char *str, *interpolationStr = NULL, *endCharacter;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;
    static char *kwlist[] = { "str", "accuracy", "biSectionMax", "interpolation", "infill", "safeDivide", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "sdd|sii", kwlist, &str, &accuracy, &biSectionMax, &interpolationStr, &infill, &safeDivide ) ) return( NULL );
    if( pointwiseXY_C_getInterpolationFromString( interpolationStr, &interpolation, 0 ) != 0 ) return( NULL );

    if( ( ptwXY = ptwXY_fromString( str, interpolation, biSectionMax, accuracy, &endCharacter, &status_nf ) ) == NULL )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from ptwXY_fromString: %s", nfu_statusMessage( status_nf ) ) );

    if( ( ptwXYPy = pointwiseXY_CNewInitialize( infill, safeDivide ) ) == NULL ) {
        ptwXY_free( ptwXY );
        return( NULL );
    }
    ptwXYPy->ptwXY = ptwXY;
    return( Py_BuildValue( "(O,s)", ptwXYPy, endCharacter ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_ptwXY_interpolatePoint( PyObject *self, PyObject *args ) {

    double x, y, x1, y1, x2, y2;
    nfu_status status;
    ptwXY_interpolation interpolation;
    char *interpolationStr;

    PyArg_ParseTuple( args, "sddddd", &interpolationStr, &x, &x1, &y1, &x2, &y2 );
    if( pointwiseXY_C_getInterpolationFromString( interpolationStr, &interpolation, 0 ) != 0 ) return( NULL );
    status = ptwXY_interpolatePoint( interpolation, x, &y, x1, y1, x2, y2 );
    if( status != nfu_Okay ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( 
            "Interpolate failed: interpolation = '%s', x = %g, (x1,y1) = (%g,%g) and (x2,y2) = (%g,%g)", interpolationStr, x, x1, y1, x2, y2 ) );
    }
    return( Py_BuildValue( "d", y ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_gaussian( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    double accuracy, xMin, xMax, offset = 0., sigma = 1., amplitude = 1., dullEps = 0.;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "accuracy", "xMin", "xMax", "offset", "sigma", "amplitude", "dullEps", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "ddd|dddd", kwlist, &accuracy, &xMin, &xMax, &offset, &sigma, 
        &amplitude, &dullEps) ) return( NULL );

    if( ( n = ptwXY_createGaussian( accuracy, offset, sigma, amplitude, xMin, xMax, dullEps, &status_nf ) ) == NULL ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from gaussian: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_basicGaussian( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );
    if( ( n = ptwXY_createGaussianCenteredSigma1( accuracy, &status_nf ) ) == NULL ) {
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from basicGaussian: %s", nfu_statusMessage( status_nf ) ) );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_unitbaseInterpolate( pointwiseXY_CPy *self, PyObject *args ) {

    double w, lw, uw;
    nfu_status status_nf;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy, *lXY, *uXY;
    int status;

    if( !PyArg_ParseTuple( args, "ddOdO", &w, &lw, &lXY, &uw, &uXY ) ) return( NULL );

    if( ( status = PyObject_IsInstance( (PyObject* ) lXY, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "lXY instance is not an pointwiseXY_C instance" ) );

    if( ( status = PyObject_IsInstance( (PyObject* ) uXY, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "uXY instance is not an pointwiseXY_C instance" ) );

    if( ( n = ptwXY_unitbaseInterpolate( w, lw, lXY->ptwXY, uw, uXY->ptwXY, &status_nf ) ) == NULL )
        return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from unitbaseInterpolate: %s", nfu_statusMessage( status_nf ) ) );
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static int pointwiseXY_C_PythonXYPairToCPair( PyObject *XYPairPy, double *x, double *y, int64_t index ) {
/*
*   Status is 0) if all is ok.
*             1) if XYPairPy does not support iteration
*             2) if the list is empty (i.e., no x object).
*             3) if the list contains only one entry (i.e., no y object).
*             4) if the list contains more than 2 entries.
*             5) if x entry could not be converted to a float.
*             6) if y entry could not be converted to a float.
*/
    int status = 1;
    PyObject *iterator = PyObject_GetIter( XYPairPy ), *xPy = NULL, *yPy = NULL, *theEnd = NULL;

    if( iterator != NULL ) {
        status = 2;
        xPy = PyIter_Next( iterator );
    }
    if( xPy != NULL ) {
        status = 3;
        yPy = PyIter_Next( iterator );
        if( yPy != NULL ) {
            status = 4;
            theEnd = PyIter_Next( iterator );
        }
    }
    if( ( xPy != NULL ) && ( yPy != NULL ) && ( theEnd == NULL ) ) {
        status = 0;
        if( pointwiseXY_C_PyNumberToFloat( xPy, x ) != 0 ) status = 5;
        if( ( status == 0 ) && ( pointwiseXY_C_PyNumberToFloat( yPy, y ) != 0 ) ) status = 6;
    }
    if( iterator != NULL ) { Py_DECREF( iterator ); }
    if( xPy != NULL ) { Py_DECREF( xPy ); }
    if( yPy != NULL ) { Py_DECREF( yPy ); }
    if( theEnd != NULL ) { Py_DECREF( theEnd ); }
    if( status != 0 ) {
        switch( status ) {
        case 1 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld does not support iteration", index );
            break;
        case 2 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld is empty", index );
            break;
        case 3 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld only has one entry", index );
            break;
        case 4 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld has more than two entries", index );
            break;
        case 5 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "x value for element at index = %ld is not a number", index );
            break;
        case 6 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "y value for element at index = %ld is not a number", index );
            break;
        }
    }
    return( status );
}
/*
************************************************************
*/
static int64_t pointwiseXY_C_PythonXYListToCList( PyObject *XYListPy, double **xys ) {

    int status = 0;
    nfu_status status_nf;
    int64_t i = 0, length, n;
    double *d;
    PyObject *item, *iterator;
    pointwiseXY_CPy *ptwXY = (pointwiseXY_CPy *) XYListPy;

    *xys = NULL;
    if( PyObject_TypeCheck( XYListPy, &pointwiseXY_CPyType ) ) {
        length = ptwXY->ptwXY->length;
        if( length != (int64_t) 0 ) {
            *xys = (double *) malloc( 2 * (size_t) length * sizeof( double ) );
            if( *xys == NULL ) return( -1 );
            if( ( status_nf = ptwXY_copyToC_XY( ptwXY->ptwXY, 0, length, length, &n, *xys ) ) != nfu_Okay ) {
                free( *xys );
                return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "pointwiseXY_C error: %s", nfu_statusMessage( status_nf ) ) );
            }
        } }
    else {
        if( ( iterator = PyObject_GetIter( XYListPy ) ) == NULL ) return( -1 );
        length = (int64_t) PySequence_Size( XYListPy );
        if( length != (int64_t) 0 ) {
            *xys = d = (double *) malloc( 2 * (size_t) length * sizeof( double ) );
            if( d == NULL ) {
                PyErr_NoMemory( );
                Py_DECREF( iterator );
                return( -1 );
            }

            while( ( status == 0 ) && ( ( item = PyIter_Next( iterator ) ) != NULL ) ) {
                status = pointwiseXY_C_PythonXYPairToCPair( item, d, &(d[1]), i );
                Py_DECREF( item );
                i++;
                d += 2;
            }

            if( status != 0 ) {
                free( *xys );
                length = -1;
            }
        }
        Py_DECREF( iterator );
    }

    return( length );
}
/*
************************************************************
*/
static int64_t pointwiseXY_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds, int ascending ) {

    int status = 0;
    int64_t length, i;
    double *d, dPrior = 0.;
    PyObject *item, *iterator;

    *ds = NULL;
    if( ( iterator = PyObject_GetIter( PyDoubleList ) ) == NULL ) return( -1 );
    if( ( length = (int64_t) PySequence_Size( PyDoubleList ) ) != (int64_t) 0 ) {
        if( ( *ds = (double *) malloc( (size_t) length * sizeof( double ) ) ) == NULL ) {
            PyErr_NoMemory( );
            Py_DECREF( iterator );
            return( -1 );
        }
        for( i = 0, d = *ds, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i++, d++ ) {
            if( ( status = pointwiseXY_C_PyNumberToFloat( item, d ) ) != 0 )
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i );
            if( ( status == 0 ) && ascending && ( i != 0 ) ) {
                if( *d <= dPrior ) {
                    status = -1;
                    pointwiseXY_C_SetPyErrorExceptionReturnNull( "data not in ascending order at index %d and %d", (int) i - 1, (int) i );
                }
            }
            dPrior = *d;
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
static PyObject *pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwXPoints *ptwX ) {

    int64_t i;
    PyObject *l = PyList_New( 0 ), *ni;

    if( l == NULL ) return( NULL );
    for( i = 0; i < ptwX->length; i++ ) {
        ni = Py_BuildValue( "d", ptwX->points[i] );
        if( pointwiseXY_C_addedItemToPythonList( l, ni ) != 0 ) return( NULL );
    }
    return( l );
}
/*
************************************************************
*/
static ptwXPoints *pointwiseXY_C_PyFloatList_to_ptwXPoints( PyObject *PyFloatList ) {

    int status = 0;
    int64_t i, length;
    double d;
    PyObject *item, *iterator;
    nfu_status status_nf;
    ptwXPoints *ptwX;

    if( ( iterator = PyObject_GetIter( PyFloatList ) ) == NULL ) return( NULL );
    length = (int64_t) PySequence_Size( PyFloatList );
    if( ( ptwX = ptwX_new( length, &status_nf ) ) == NULL ) {
        PyErr_NoMemory( ); }
    else {
        ptwX->length = length;
        for( i = 0, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i++ ) {
            if( ( status = pointwiseXY_C_PyNumberToFloat( item, &d ) ) == 0 ) {
                ptwX->points[i] = d; }
            else {
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i );
            }
            Py_DECREF( item );
        }
    }
    Py_DECREF( iterator );
    if( status != 0 ) ptwX = ptwX_free( ptwX );

    return( ptwX );
}
/*
************************************************************
*/
static int pointwiseXY_C_PyNumberToFloat( PyObject *n, double *d ) {

    if( PyFloat_Check( n ) ) {
        *d = PyFloat_AsDouble( n ); }
    else if( PyInt_Check( n ) ) {
        *d = (double) PyInt_AsLong( n ); }
    else if( PyLong_Check( n ) ) {
        *d = (double) PyLong_AsLongLong( n ); }
    else {
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
static int pyObject_NumberOrPtwXY(  PyObject *other ) {

    if( PyInt_Check( other ) ) {
        return( pyObject_Number ); }
    else if( PyLong_Check( other ) ) {
        return( pyObject_Number ); }
    else if( PyFloat_Check( other ) ) {
        return( pyObject_Number ); }
    else if( PyList_Check( other ) ) {
        return( pyObject_List ); }
    else if( PyTuple_Check( other ) ) {
        return( pyObject_Tuple ); }
    else if( PyObject_TypeCheck( other, &pointwiseXY_CPyType ) ) {
        return( pyObject_ptwXY );
    }
    return( pyObject_Unsupported );
}
/*
************************************************************
*/
static int pointwiseXY_C_addedItemToPythonList( PyObject *list, PyObject *item ) {

    if( item == NULL ) {
        Py_DECREF( list );
        return( -1 ); }
    else {
        if( PyList_Append( list, item ) != 0 ) {
            Py_DECREF( item );          /* append created a new ref */
            Py_DECREF( list );
            return( -1 );
        }
        Py_DECREF( item );              /* append created a new ref */
    }
    return( 0 );
}
/*
************************************************************
*/
static int isOkayAndHasData( ptwXYPoints *ptwXY ) {

    if( ptwXY->status != nfu_Okay ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "pointwiseXY_C object had a prior error and is not usable" ) );
    if( ptwXY->length == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "pointwiseXY_C object has no data" ) );
    return( 0 );
}
/*
************************************************************
*/
static int pointwiseXY_C_getInterpolationFromString( char *interpolationStr, ptwXY_interpolation *interpolation, int allowOther ) {

    int xInterpolation = 0, yInterpolation = 0, n;
    char *c;

    if( interpolationStr != NULL ) {
        if( strcmp( interpolationStr, "flat" ) == 0 ) {
            *interpolation = ptwXY_interpolationFlat;
            return( 0 );
        }
        if( strcmp( interpolationStr, "other" ) == 0 ) goto other;
        if( ( c = strchr( interpolationStr, ',' ) ) == NULL )
            return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "interpolation = '%s' missing separator ','", interpolationStr ) );
        n = (int) ( c - interpolationStr );
        if( strncmp( interpolationStr, "other", n ) == 0 ) goto other;
        if( strncmp( interpolationStr, "linear", n ) != 0 ) {
            if( strncmp( interpolationStr, "log", n ) != 0 ) 
                return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid x interpolation in '%s'", interpolationStr ) );
            xInterpolation = 1;
        }
        c++;
        if( strcmp( c, "other" ) == 0 ) goto other;
        if( strcmp( c, "linear" ) != 0 ) {
            if( strcmp( c, "log" ) != 0 ) {
                if( strcmp( c, "flat" ) != 0 )
                    return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid y interpolation in '%s'", interpolationStr ) );
                *interpolation = ptwXY_interpolationFlat;
                return( 0 );
            }
            yInterpolation = 1;
        }
        if( xInterpolation == 0 ) {
            *interpolation = ptwXY_interpolationLinLin;
            if( yInterpolation != 0 ) *interpolation = ptwXY_interpolationLinLog; }
        else {
            *interpolation = ptwXY_interpolationLogLin;
            if( yInterpolation != 0 ) *interpolation = ptwXY_interpolationLogLog;
        }
    }
    return( 0 );

other:
    if( !allowOther ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "interpolation 'other' not allowed here" ) );
    *interpolation = ptwXY_interpolationOther;
    return( 0 );
}
/*
************************************************************
*/
static int pointwiseXY_C_getInterpolationFromTupleOfTwoStrings( PyObject *twoInterpolations, ptwXY_interpolation *interpolation ) {

    int status = -1;
    PyObject *iterator, *item1 = NULL, *item2 = NULL, *item3 = NULL;
    enum e_interpolationType i1 = e_interpolationTypeInvalid, i2 = e_interpolationTypeInvalid;
 
    if( ( iterator = PyObject_GetIter( twoInterpolations ) ) == NULL ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Interpolation object is not a sequence" ) );
    if( ( item1 = PyIter_Next( iterator ) ) != NULL ) {
        i1 = pointwiseXY_C_getInterpolationTypeOfObject( item1 );
        if( ( item2 = PyIter_Next( iterator ) ) != NULL ) {
            i2 = pointwiseXY_C_getInterpolationTypeOfObject( item2 );
            if( ( item3 = PyIter_Next( iterator ) ) != NULL ) status = -2;
            if( item3 != NULL ) { Py_DECREF( item3 ); }
            Py_DECREF( item2 );
        }
        Py_DECREF( item1 );
    }
    Py_DECREF( iterator );
    if( status == -2 ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "More than two elements in interpolation string list" ) );
    if( ( i1 == e_interpolationTypeInvalid ) || ( i2 == e_interpolationTypeInvalid ) ) 
        return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "More than two elements in interpolation string list" ) );
    status = 0;
    if( i2 == e_interpolationTypeFlat ) {
        *interpolation = ptwXY_interpolationFlat; }
    else {
        *interpolation = ptwXY_interpolationOther;
        if( i1 == e_interpolationTypeLinear ) {
            if( i2 == e_interpolationTypeLinear ) {
                *interpolation = ptwXY_interpolationLinLin; }
            else if( i2 == e_interpolationTypeLog ) {
                *interpolation = ptwXY_interpolationLinLog;
            } }
        else if( i1 == e_interpolationTypeLog ) {
            if( i2 == e_interpolationTypeLinear ) {
                *interpolation = ptwXY_interpolationLogLin; }
            else if( i2 == e_interpolationTypeLog ) {
                *interpolation = ptwXY_interpolationLogLog;
            }
        }
    }
    return( ( status < 0 ) ? -1 : 0 );
}
/*
************************************************************
*/
static enum e_interpolationType pointwiseXY_C_getInterpolationTypeOfObject( PyObject *o ) {

    char *p = PyString_AsString( o );

    if( p == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Interpolation object is not a string" );
        return( e_interpolationTypeInvalid );
    }
    if( strcmp( "linear", p ) == 0 ) return( e_interpolationTypeLinear );
    if( strcmp( "log", p ) == 0 ) return( e_interpolationTypeLog );
    if( strcmp( "flat", p ) == 0 ) return( e_interpolationTypeFlat );
    if( strcmp( "other", p ) == 0 ) return( e_interpolationTypeOther );
    pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "Invalid interpolation string = '%s'", p );
    return( e_interpolationTypeInvalid );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_GetNone( void ) {

    Py_INCREF( Py_None );
    return( Py_None );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

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
static int pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... ) {

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
static PyMethodDef pointwiseXY_CPyMethods[] = {

    { "allocatedSize", (PyCFunction) pointwiseXY_C_allocatedSize, METH_NOARGS, "Returns the size of memory allocated in the points region." },
    { "applyFunction", (PyCFunction) pointwiseXY_C_applyFunction, METH_VARARGS | METH_KEYWORDS, \
        "Returns a new instance which is y(x) = f(y_s(x)) where f is a function of one variable and y_s(x) is self's.\n" \
        "y-value at x. The function f must take two arguments. The first is y_s(x) and the second (not including self)\n" \
        "is the second argument to applyFunction. The returned object may contain more points then self as applyFunction\n" \
        "does infilling. That is, points are added recursively between points until accuracy or biSectionMax is reached.\n" \
        "\nArguments are:\n" \
        "   f             a function that returns the new y(x) given the old y(x),\n" \
        "   parameters    any python object. Passed as the second argument to f,\n" \
        "   accuracy      [o] The accuracy for infilling. If not present, taken from self,\n" \
        "   biSectionMax  [o] The maximum number of bi-sections for infilling. If not present, taken from self.\n" \
        "   checkForRoots [o] If true, and biSectionMax > 0, an addition point is added whenever two consecutive points cross the y-axis\n" \
        "                     at (or close to) the root (i.e., the y crossing point).\n" \
        "\nTo return f(x), set self's y-values to x (i.e., y(x) = x for self)." },
    { "changeInterpolation", (PyCFunction) pointwiseXY_C_changeInterpolation, METH_VARARGS | METH_KEYWORDS,
        "Returns a new instance that is equivalent to self, but with the interpolation changed to interpolation.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   interpolation  [o] the new interpolation (default is 'linear,linear'),\n" \
        "   accuracy       [o] the accuracy of the conversion from the old to the new interpolation (default is self's accuracy),\n" \
        "   lowerEps       [o] for flat to linear, the lower eps,\n" \
        "   upperEps       [o] the flat to linear, the upper eps. lowerEps and upperEps cannot both be 0.\n" },
    { "changeInterpolationIfNeeded", (PyCFunction) pointwiseXY_C_changeInterpolationIfNeeded, METH_VARARGS | METH_KEYWORDS,
        "Returns self if it is one of the allowed interpolations. Otherwise, a new instance is returned that is equivalent to self,\n" \
        "but with the interpolation changed to the interpolation of the first allowed interpolations.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   allowedInterpolations the list of allowed interpolations,\n" \
        "   accuracy              [o] the accuracy of the conversion from the old to the new interpolation (default is self's accuracy),\n" \
        "   lowerEps              [o] for flat to linear, the lower eps,\n" \
        "   upperEps              [o] the flat to linear, the upper eps. lowerEps and upperEps cannot both be 0.\n" },
    { "clip", (PyCFunction) pointwiseXY_C_clip, METH_VARARGS | METH_KEYWORDS,
        "Returns a new instance which is the same as self, but with the y-values clipped between yMin and yMax\n" \
        "clip may add points, to insure that the return instance has the same shape as self between yMin and yMax\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   yMin  [o]   the lower y-value for clipping,\n" \
        "   yMax  [o]   the upper y-value for clipping.\n" },
    { "coalescePoints", (PyCFunction) pointwiseXY_C_coalescePoints, METH_VARARGS, "Moves all points in overflow region to points region." },
    { "convolute", (PyCFunction) pointwiseXY_C_convolute, METH_VARARGS | METH_KEYWORDS, \
        "Returns a new instance which is the convolution of self with the first argument, that must also be a pointwiseXY_C instance." },
    { "copy", (PyCFunction) pointwiseXY_C_copy, METH_NOARGS, "Returns a copy of self." },
    { "copyDataToXYs", (PyCFunction) pointwiseXY_C_copyDataToXYs, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of self's data as [ x1, y1 ], ... [ xn, yn ] ].\n" \
        "\nArguments are:\n" \
        "   xScale [o]  a float to scale all x-values by,\n" \
        "   yScale [o]  a float to scale all y-values by." },
    { "copyDataToXsAndYs", (PyCFunction) pointwiseXY_C_copyDataToXsAndYs, METH_VARARGS | METH_KEYWORDS, 
        "Returns python list of length two. The first element is a python list of\n" \
        "self's x-values and the second element is a python list of self's y-values.\n" \
        "\nArguments are:\n" \
        "   xScale [o]  a float to scale all x-values by,\n" \
        "   yScale [o]  a float to scale all y-values by." },
    { "dullEdges", (PyCFunction) pointwiseXY_C_dullEdges, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new instance that is a copy of self, except the endpoints are guaranteed to have 0's for y-values\n" \
        "as long as lowerEps and/or upperEps are none zero (see below). The following algorithm is used.\n\n" \
        "If lowerEps is zero or the y-value at xMin is 0, then nothing is added at the lower end. Otherwise, the following is\n"\
        "done at the lower end.\n" \
        "   If lowerEps is positive, a point 'abs( xMin * lowerEps )' above xMin is added with its interpolated y-value; provided,\n" \
        "       this x value is 'abs( xMin * lowerEps )' less than the next x-value. In the prior sentence, if xMin is 0, then\n" \
        "       replace 'abs( xMin * lowerEps )' with 'abs( lowerEps )'. Independent of whether a point is added, xMin's y-value is\n" \
        "       set to 0.\n" \
        "   If lowerEps is negative, the logic for adding the point above xMin for positive lowerEps is followed. In addition, a\n" \
        "       point is added 'abs( xMin * lowerEps )' below xMin with a value of zero and the point at xMin is reset by\n" \
        "       interpolating the new surrounding values. However, if positiveXOnly is True and the point below xMin would cause\n" \
        "       a negative x-value (and xMin is not negative) then the logic for positive lowerEps is implemented instead.\n" \
        "\n" \
        "   The logic for upperEps is similar to lowerEps except, replace xMin with xMax, below with above and above with below.\n" \
        "       Also positiveXOnly is ignored.\n" \
        "\nArguments are:\n" \
        "   lowerEps       [o] a point (or two if lowerEps is negative) is (are) added a distance xMin * lowerEps from xMin,\n" \
        "   upperEps       [o] a point (or two if upperEps is negative) is (are) added a distance xMax * upperEps from xMax,\n" \
        "   positiveXOnly  [o] this only applies to lowerEps and only if an added point would be negative when xMin is non-negative." },
    { "exp", (PyCFunction) pointwiseXY_C_exp, METH_VARARGS,
        "self.exp( a )\n\n" \
        "Returns a new instance with ts y-values set to exp( a * [self's y-values] ). x-values are added to meet required accuaracy."},
    { "getAccuracy", (PyCFunction) pointwiseXY_C_getAccuracy, METH_NOARGS, "Returns self's accuracy value." },
    { "getBiSectionMax", (PyCFunction) pointwiseXY_C_getBiSectionMax, METH_NOARGS, "Returns self's biSectionMax value." },
    { "getInfill", (PyCFunction) pointwiseXY_C_getInfill, METH_NOARGS, "Returns self's infill flag." },
    { "getInterpolation", (PyCFunction) pointwiseXY_C_getInterpolation, METH_NOARGS, "Returns self's interpolation as a string." },
    { "getSecondaryCacheSize", (PyCFunction) pointwiseXY_C_getSecondaryCacheSize, METH_NOARGS, "Returns the size of self's secondary cache." },
    { "getSafeDivide", (PyCFunction) pointwiseXY_C_getSafeDivide, METH_NOARGS, "Returns self's safeDivide flag." },
    { "getValue", (PyCFunction) pointwiseXY_C_getValue, METH_VARARGS, "Gets the y value at x." },
    { "getUserFlag", (PyCFunction) pointwiseXY_C_getUserFlag, METH_NOARGS, "Gets the user flag." },
    { "groupOneFunction", (PyCFunction) pointwiseXY_C_groupOneFunction, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of self between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "groupTwoFunctions", (PyCFunction) pointwiseXY_C_groupTwoFunctions, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of the product of self and f2 between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   f2                  the second pointwiseXY_C function with the integrand being the product of self * f2\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "groupThreeFunctions", (PyCFunction) pointwiseXY_C_groupThreeFunctions, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of the product of self, f2 and f3 between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   f2                  the second pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   f3                  the third pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "integrate", (PyCFunction) pointwiseXY_C_integrate, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self from xMin to xMax.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xMin    [o] the lower limit of the integral, default is xMin of self,\n" \
        "   xMax    [o] the upper limit of the integral, default is xMax of self." },
    { "integrateWithWeight_x", (PyCFunction) pointwiseXY_C_integrateWithWeight_x, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self weighted by x from xMin to xMax (i.e. integral dx x * self(x)).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xMin    [o] the lower limit of the integral, default is xMin of self,\n" \
        "   xMax    [o] the upper limit of the integral, default is xMax of self." },
    { "integrateWithWeight_sqrt_x", (PyCFunction) pointwiseXY_C_integrateWithWeight_sqrt_x, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self weighted by sqrt( x ) from xMin to xMax (i.e. integral dx sqrt( x ) * self(x)).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xMin    [o] the lower limit of the integral, default is xMin of self,\n" \
        "   xMax    [o] the upper limit of the integral, default is xMax of self." },
    { "normalize", (PyCFunction) pointwiseXY_C_normalize, METH_NOARGS, 
        "Returns a new instance with the same x-values as self but with the y-values scaled so that the area of the curve is 1." },
    { "runningIntegral", (PyCFunction) pointwiseXY_C_runningIntegral, METH_NOARGS, 
        "Returns, as a python list, the integrals between successive x-values." },
    { "areDomainsMutual", (PyCFunction) pointwiseXY_C_areDomainsMutual, METH_VARARGS, 
        "This routine returns True if self and the first argument have a mutual domain, and false otherwise." },
    { "mutualify", (PyCFunction) pointwiseXY_C_mutualify, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list containing two pointwiseXY_C instances that are the mutualification of self and other.\n" \
        "Mutualification is the act of modifying, if needed, two pointwiseXY_C instances so that their domains are mutual.\n" \
        "Self and other are not altered.\n" \
        "\nArguments are:\n" \
        "   lowerEps1       the lowerEps applied to self if needed, see dullEdges for meaning,\n" \
        "   upperEps1       the upperEps applied to self if needed, see dullEdges for meaning,\n" \
        "   positiveXOnly1  the positiveXOnly applied to self if needed, see dullEdges for meaning,\n" \
        "   other           the second pointwiseXY_C instance to mutualify self with,\n" \
        "   lowerEps2       the lowerEps applied to other if needed, see dullEdges for meaning,\n" \
        "   upperEps2       the upperEps applied to other if needed, see dullEdges for meaning,\n" \
        "   positiveXOnly2  the positiveXOnly applied to other if needed, see dullEdges for meaning,\n" },
    { "overflowAllocatedSize", (PyCFunction) pointwiseXY_C_overflowAllocatedSize, METH_NOARGS, 
        "Returns the size of memory allocated for the overflow region." },
    { "overflowLength", (PyCFunction) pointwiseXY_C_overflowLength, METH_NOARGS, "Returns of number of points in the overflow region." },
    { "plot", (PyCFunction) pointwiseXY_C_plot, METH_NOARGS, 
        "Calls Gnuplot, if it exists, to plot self. This is a simple method, mainly for debugging.\n" \
        "\nArguments: None." },
    { "reallocatePoints", (PyCFunction) pointwiseXY_C_reallocatePoints, METH_VARARGS | METH_KEYWORDS, \
        "self.reallocatePoints( size, forceSmaller = True )\n\n" \
        "Adjusts the memory allocated for primary points to the maximum of size and the current length of self.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   size              the desired allocated size of self (actual size will be larger if length is greater than size),\n" \
        "   forceSmaller  [o] if False action is only taken if allocated is significantly greater than size (default is True).\n" },
    { "reallocateOverflowPoints", (PyCFunction) pointwiseXY_C_reallocateOverflowPoints, METH_VARARGS, \
        "self.reallocateOverflowPoints( size )\n\n" \
        "Adjusts the memory allocated for the overflow points to size (the first and only argument)." },
    { "setData", (PyCFunction) pointwiseXY_C_setData, METH_VARARGS, "Replaces the data in self with the frist argument. This argument must be a list. " \
        "Each item of the list must contain two floats, or objects that can be convert to float (e.g., [ [ 1, 2 ], [ 2, 4 ], [ 4, 0.5 ] ]" },
    { "setAccuracy", (PyCFunction) pointwiseXY_C_setAccuracy, METH_VARARGS, "Sets self's accuracy value and returns the actual value set." },
    { "setBiSectionMax", (PyCFunction) pointwiseXY_C_setBiSectionMax, METH_VARARGS, "Sets self's biSectionMax value and returns the actual value set." },
    { "setSecondaryCacheSize", (PyCFunction) pointwiseXY_C_setSecondaryCacheSize, METH_VARARGS, "Sets the size of self's secondary cache to size." },
    { "setDataFromList", (PyCFunction) pointwiseXY_C_setDataFromList, METH_VARARGS, 
        "Replaces the data in self with the python list given by the first argument.\n" \
        "\nArguments are:\n" \
        "   xys      the list of xy-values as a single python list (e.g., [ 1, 2, 4, 3, 5, 6]\n" },
    { "setDataFromXsAndYs", (PyCFunction) pointwiseXY_C_setDataFromXsAndYs, METH_VARARGS, 
        "Replaces the data in self with the 2 python lists given by the first 2 arguments.\n" \
        "\nArguments are:\n" \
        "   xs      the list of x-values\n" \
        "   ys      the list of y-values\n" },
    { "setInfill", (PyCFunction) pointwiseXY_C_setInfill, METH_VARARGS, "Sets self's infill flag to the first argument." },
    { "setSafeDivide", (PyCFunction) pointwiseXY_C_setSafeDivide, METH_VARARGS, "Sets self's safeDivide flag to the first argument." },
    { "setValue", (PyCFunction) pointwiseXY_C_setValue, METH_VARARGS, "Sets the y value at x." },
    { "setUserFlag", (PyCFunction) pointwiseXY_C_setUserFlag, METH_VARARGS, "Sets the users flag via the first argument." },
    { "showInteralStructure", (PyCFunction) pointwiseXY_C_showInteralStructure, METH_VARARGS | METH_KEYWORDS, 
        "For debbuging only. Dumps information on internal data." },
    { "thicken", (PyCFunction) pointwiseXY_C_thicken, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new instance with denser points than self filled in using self's interpolation. The number of points added are determined by the following arguments\n" \
            "sectionSubdivideMax    maximum number of points to insert between consecutive points (default 1),\n" \
            "dxMax                  minimum dx step (default 0),\n" \
            "fxMax                  minimum fractional step (default 1)." },
    { "thin", (PyCFunction) pointwiseXY_C_thin, METH_VARARGS, "Returns a new instance with points of self thinned to accuracy of argument one." },
    { "toString", (PyCFunction) pointwiseXY_C_toString, METH_VARARGS | METH_KEYWORDS, "Returns a string representation of self." \
        " This method has three keyword parameters:\npairsPerLine, format and pairSeparator which are defined as,\n" \
        "    pairsPerLine    the number of pairs to put on each line\n" \
        "    format          a valid format to convert an (x,y) pair (i.e., two floats) into a string (e.g. format = ' %.3f %12.5e')\n" \
        "    pairSeparator   a string to put between every pair (e.g, to put a comma to separate pairs use pairSeparator = ',')" },
    { "trim", (PyCFunction) pointwiseXY_C_trim, METH_NOARGS, "Returns a new instance with excessive 0. y-value points at beginning and end of self removed." },
    { "union", (PyCFunction) pointwiseXY_C_union, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new pointwiseXY_C object whose x values are the union of self and other." },
    { "mergeClosePoints", (PyCFunction) pointwiseXY_C_mergeClosePoints, METH_VARARGS, 
        "Returns a new pointwiseXY_C object whose x values are the union of self and other." },
    { "getDomainGrid", (PyCFunction) pointwiseXY_C_getDomainGrid, METH_VARARGS, "Returns a list of x-values for self." },
    { "xMin", (PyCFunction) pointwiseXY_C_xMin, METH_NOARGS, "Returns the first point's x value." },
    { "xMax", (PyCFunction) pointwiseXY_C_xMax, METH_NOARGS, "Returns the last point's x value." },
    { "xSlice", (PyCFunction) pointwiseXY_C_xSlice, METH_VARARGS | METH_KEYWORDS, "Returns a new instance with self sliced between xMin and xMax.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xMin    [o] the lower x-value of the slice, default is xMin of self,\n" \
        "   xMax    [o] the upper x-value of the slice, default is xMax of self,\n" \
        "   fill    [o] if True, points are added at xMin and xMax if they are not in self,\n" \
        "               else only existing points in the range [xMin, xMax] are included.\n" \
        "   dullEps [o] (Currently not implemented) the lower and upper points are dulled, default is 0." },
    { "yMin", (PyCFunction) pointwiseXY_C_yMin, METH_NOARGS, "Returns the minimum y-value in self or 0 if self is empty." },
    { "yMax", (PyCFunction) pointwiseXY_C_yMax, METH_NOARGS, "Returns the maximum y-value in self or 0 if self is empty." },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
static PyTypeObject pointwiseXY_CPyType = {
    PyObject_HEAD_INIT( NULL )
    0,                                      /* ob_size        */
    "pointwiseXY_C.pointwiseXY_C",          /* tp_name        */
    sizeof( pointwiseXY_CPy ),              /* tp_basicsize   */
    0,                                      /* tp_itemsize    */
    /* methods */ 
    (destructor) pointwiseXY_C_dealloc,     /* tp_dealloc     */
    0,                                      /* tp_print       */
    0,                                      /* tp_getattr     */
    0,                                      /* tp_setattr     */
    0,                                      /* tp_compare     */
    (reprfunc) pointwiseXY_C__repr__,       /* tp_repr        */
    &pointwiseXY_CPy_number,                /* tp_as_number   */
    &pointwiseXY_CPy_sequence,              /* tp_as_sequence */
    0,                                      /* tp_as_mapping  */
    0,                                      /* tp_hash        */
    0,                                      /* tp_call        */
    0,                                      /* tp_str         */
    0,                                      /* tp_getattro    */
    0,                                      /* tp_setattro    */
    0,                                      /* tp_as_buffer   */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, /* tp_flags     */
    pointwiseXY_C__doc__,                   /* tp_doc         */
    0,                                      /* tp_traverse    */
    0,                                      /* tp_clear       */
    0,                                      /* tp_richcompare */
    0,                                      /* tp_weaklistoffset */
    0,                                      /* tp_iter        */
    0,                                      /* tp_iternext    */
    pointwiseXY_CPyMethods,                 /* tp_methods     */
    0,                                      /* tp_members     */
    0,                                      /* tp_getset      */
    0,                                      /* tp_base        */
    0,                                      /* tp_dict        */
    0,                                      /* tp_descr_get   */
    0,                                      /* tp_descr_set   */
    0,                                      /* tp_dictoffset  */
    (initproc) pointwiseXY_C__init__,       /* tp_init        */
    0,                                      /* tp_alloc       */
    0                                       /* tp_new         */
};
/*
************************************************************
*/

static PyMethodDef pointwiseXY_CMiscPyMethods[] = {

    { "defaultAccuracy", (PyCFunction) pointwiseXY_C_defaultAccuracy, METH_NOARGS,
        "defaultAccuracy( )\n\nReturns the value of the default accuracy for the __init_ method." },
    { "createFromFunction", (PyCFunction) pointwiseXY_C_createFromFunction, METH_VARARGS | METH_KEYWORDS,
        "createFromFunction( xs, f, parameters, accuracy, biSectionMax, checkForRoots = False, infill = True, saveDivide = True )\n\n" \
        "Returns a pointwiseXY_C instance which represents the function f(x) where f is the second argument and must be a\n" \
        "function reference. The function f must take two arguments. The first arugment to f is x and the second is the third\n" \
        "argument to createFromFunction. The function f must return, as a float, the y-value at x. The returned pointwiseXY_C\n" \
        "instance may contain more points then xs (the first argument) as createFromFunction does infilling. That is, points\n" \
        "are added recursively between the points in xs until accuracy or biSectionMax is reached.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xs            An ascending list of x-values to use as an initial guess for the refinement of the returned pointwiseXY_C instance.\n" \
        "   f             A function that returns the y-value given an x-value.\n" \
        "   parameters    Any python object which may be needed by f. Passed as the second argument to f.\n" \
        "   accuracy      The accuracy for infilling.\n" \
        "   biSectionMax  The maximum number of bi-sections for infilling.\n" \
        "   checkForRoots [o] If true, and biSectionMax > 0, an additional point is added whenever two consecutive points cross the y-axis\n" \
        "                     at (or close to) the root (i.e., the y crossing point) (default is False).\n" \
        "   infill        [o] Infill value used for the returned pointwiseXY_C instance (default is True).\n" \
        "   safeDivide    [o] SafeDivide value used for the returned pointwiseXY_C instance (default is True).\n" \
        "\nExample: Return a pointwise representation of 'x * sin( x**2 )' in the domain [1, 10].\n\n" \
        "import math\n" \
        "import pointwiseXY_C\n" \
        "def f( x, args ) :\n" \
        "   return( x * math.sin( x * x ) )\n" \
        "\n" \
        "xSin_xx = pointwiseXY_C.createFromFunction( ( 1, 10 ), f, None, 1e-3, 12 )\n" \
        "\n\n# A better solution may be to input the known zeros of f and force them to be 0., as\n" \
        "\n" \
        "xs = [ 1 ] + [ math.sqrt( i * math.pi ) for i in xrange( 1, int( math.sqrt( 100 / math.pi ) ) ) ] + [ 10 ]\n" \
        "xSin_xx = pointwiseXY_C.createFromFunction( xs, f, None, 1e-3, 12, checkForRoots = True )" },
    { "createFromString", (PyCFunction) pointwiseXY_C_createFromString, METH_VARARGS | METH_KEYWORDS,
        "createFromString( str, accuracy, biSectionMax, interpolation = 'linear,linear', infill = True, safeDivide = True )\n\n" \
        "Returns a tuple of two elements. The first element is a pointwiseXY_C instance representing the float values\n" \
        "translated from 'str'. The second element is the portion of 'str' not translated\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   str           The containing a list of floats to be converted.\n" \
        "   accuracy      The accuracy for infilling.\n" \
        "   biSectionMax  The maximum number of bi-sections for infilling.\n" \
        "   interpolation [o] the interpolation string (default is 'linear,linear').\n" \
        "   infill        [o] Infill value used for the returned pointwiseXY_C instance (default is True).\n" \
        "   safeDivide    [o] SafeDivide value used for the returned pointwiseXY_C instance (default is True).\n" },
    { "interpolatePoint", (PyCFunction) pointwiseXY_C_ptwXY_interpolatePoint, METH_VARARGS,
        "interpolatePoint( interpolation, x, x1, y1, x2, y2 )\n\n" \
        "Returns interpolation of x for x1 <= x <= x2 given x1, y1, x2, y2 and interpolation law." \
        "\n" \
        "Arguments are:\n" \
        "   interpolation   a string representing the interpolation law (e.g., 'log,linear'; see constructor's docstring),\n" \
        "   x               x point of interpolated y-value,\n" \
        "   x1              lower x-value,\n" \
        "   y1              y(x1),\n" \
        "   x2              upper x-value,\n" \
        "   y2              y(x2)." },
    { "gaussian", (PyCFunction) pointwiseXY_C_gaussian, METH_VARARGS | METH_KEYWORDS, 
        "gaussian( accuracy, xMin, xMax, offset = 0., sigma = 1., amplitude = 1., dullEps = False )\n\n" \
        "Returns a new pointwiseXY_C instance constructed from the following equation\n\n" \
        "       amplitude * exp( ( ( x - offset ) / sigma )^2 / 2 )        for xMin <= x <= xMax\n" \
        "\n" \
        "Arguments are:  ([o] implies optional argument)\n" \
        "   accuracy        the accuracy of linear interpolation,\n" \
        "   xMin            the lower x point generated,\n" \
        "   xMax            the upper x point generated,\n" \
        "   offset     [o]  the x offset of the center of the Gaussian,\n" \
        "   sigma      [o]  width of the Gaussian,\n" \
        "   amplitude  [o]  the Gaussian's amplitude,\n" \
        "   dullEps    [o]  currently not used." },
    { "basicGaussian", (PyCFunction) pointwiseXY_C_basicGaussian, METH_VARARGS, 
        "basicGaussian( accuracy )\n\n" \
        "Returns a new pointwiseXY_C instance constructed from the following equation\n\n" \
        "       exp( ( x^2 / 2 )\n" \
        "\n" \
        "Arguments are:\n" \
        "   accuracy        the accuracy of linear interpolation,\n" },
    { "unitbaseInterpolate", (PyCFunction) pointwiseXY_C_unitbaseInterpolate, METH_VARARGS,
        "unitbaseInterpolate( x, lw, lXY, uw, uXY )\n\n"
        "Returns the unitbase interpolation of two XYs objects at w where the axes are labeled (w, x, y).\n" \
        "\n" \
        "Arguments are:\n" \
        "   w       the point between lw and uw to return the unitbase interpolation of lXY and uXY,\n" \
        "   lw      the w point where lXY resides,\n" \
        "   lXY     a pointwiseXY_C instance for a function y(x),\n" \
        "   uw      the w point where uXY resides,\n" \
        "   uXY     a pointwiseXY_C instance for a function y(x),\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
ptwXYPoints *pointwiseXY_C_factory_get_ptwXYPoints( PyObject *ptwXY_Py ) {

    if( is_pointwiseXY_CPyObject( ptwXY_Py ) ) return( ((pointwiseXY_CPy *) ptwXY_Py)->ptwXY );
    return( (ptwXYPoints *) pointwiseXY_C_SetPyErrorExceptionReturnNull( "Object is not a pointwiseXY_C instance" ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_factory_create( ptwXYPoints *ptwXY, int infill, int safeDivide ) {

    pointwiseXY_CPy *nPy = pointwiseXY_CNewInitialize( infill, safeDivide );

    if( nPy == NULL ) {
        ptwXY_free( ptwXY ); }
    else {
        nPy->ptwXY = ptwXY;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
DL_EXPORT( void ) initpointwiseXY_C( void ) {

    PyObject *m;
    static void *Py_pointwiseXY_C_API[Py_pointwiseXY_C_API_numberOfPointers];
    PyObject *c_api_object;

    pointwiseXY_CPyType.tp_new = PyType_GenericNew;
    if( PyType_Ready( &pointwiseXY_CPyType ) < 0 ) return;

    if( ( m = Py_InitModule3( "pointwiseXY_C", pointwiseXY_CMiscPyMethods, "A module that contains the class pointwiseXY_C." ) ) == NULL ) return;

    Py_INCREF( &pointwiseXY_CPyType );
    PyModule_AddObject( m, "pointwiseXY_C", (PyObject *) &pointwiseXY_CPyType );

    Py_pointwiseXY_C_API[Py_pointwiseXY_C_get_ptwXYPoints_NUM] = (void *) pointwiseXY_C_factory_get_ptwXYPoints;
    Py_pointwiseXY_C_API[Py_pointwiseXY_C_NUM] = (void *) pointwiseXY_C_factory_create;

/*
Cannot get the PyCapsule_Import (i.e. python 2.7) way to work but the PyImport_ImportModuleEx with 2.6 and 2.7 so use it for now.
#if( PY_VERSION_HEX < 0x02070000 )
*/
#if 1
    c_api_object = PyCObject_FromVoidPtr( (void *) Py_pointwiseXY_C_API, NULL) ;
#else
    c_api_object = PyCapsule_New( (void *) Py_pointwiseXY_C_API, "pointwiseXY_C._C_API", NULL );
#endif

    if( c_api_object != NULL ) PyModule_AddObject( m, "_C_API", c_api_object );
}
