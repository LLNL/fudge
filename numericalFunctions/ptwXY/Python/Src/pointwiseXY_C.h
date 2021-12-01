/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#ifndef POINTWISEXY_C_H
#define POINTWISEXY_C_H

#ifdef __cplusplus
    extern "C" {
#endif

/* Total number of C API pointers */
#define Py_pointwiseXY_C_API_numberOfPointers 2

/* C API functions */
#define Py_pointwiseXY_C_get_ptwXYPoints_NUM 0
#define Py_pointwiseXY_C_get_ptwXYPoints_RETURN ptwXYPoints *
#define Py_pointwiseXY_C_get_ptwXYPoints_PROTO ( PyObject * )

#define Py_pointwiseXY_C_NUM 1
#define Py_pointwiseXY_C_RETURN PyObject *
#define Py_pointwiseXY_C_PROTO ( ptwXYPoints *ptwXY, int infill, int safeDivide )

#ifdef POINTWISEXY_C_MODULE         /* This section is used when compiling pointwiseXY_C.c */

static Py_pointwiseXY_C_get_ptwXYPoints_RETURN pointwiseXY_C_factory_get_ptwXYPoints Py_pointwiseXY_C_get_ptwXYPoints_PROTO;
static Py_pointwiseXY_C_RETURN pointwiseXY_C_factory_create Py_pointwiseXY_C_PROTO;

#else                               /* This section is used in modules that use pointwiseXY_C's API */

static void **Py_pointwiseXY_C_API;
#define pointwiseXY_C_factory_get_ptwXYPoints (*(Py_pointwiseXY_C_get_ptwXYPoints_RETURN (*)Py_pointwiseXY_C_get_ptwXYPoints_PROTO) Py_pointwiseXY_C_API[Py_pointwiseXY_C_get_ptwXYPoints_NUM])
#define pointwiseXY_C_factory_create (*(Py_pointwiseXY_C_RETURN (*)Py_pointwiseXY_C_PROTO) Py_pointwiseXY_C_API[Py_pointwiseXY_C_NUM])

static int import_pointwiseXY_C( void ) {  /* Return -1 on error, 0 on success.  PyCapsule_Import will set an exception if there's an error.  */
/*
Cannot get the PyCapsule_Import way to work python 2.7 but the PyImport_ImportModuleEx did work, so use it for now.
*/
#if PY_MAJOR_VERSION < 3

    PyObject *c_api_object, *module1, *module2, *fromList;

    if( ( fromList = Py_BuildValue( "(s)", "pointwiseXY_C" ) ) == NULL ) return( -1 );
    module1 = PyImport_ImportModuleEx( "numericalFunctions", NULL, NULL, fromList );
    Py_DECREF( fromList );
    if( module1 == NULL ) return( -1 );
    module2 = PyObject_GetAttrString( module1, "pointwiseXY_C" );
    Py_DECREF( module1 );
    if( module2 == NULL ) return( -1 );

    c_api_object = PyObject_GetAttrString( module2, "_C_API" );
    Py_DECREF( module2 );
    if( c_api_object == NULL ) return( -1 );

    if( PyCObject_Check( c_api_object ) ) Py_pointwiseXY_C_API = (void **) PyCObject_AsVoidPtr( c_api_object );
    Py_DECREF( c_api_object );

    return( 0 );

#else

    Py_pointwiseXY_C_API = (void **) PyCapsule_Import( "numericalFunctions.pointwiseXY_C._C_API", 0 );
    return( ( Py_pointwiseXY_C_API != NULL ) ? 0 : -1 );

#endif

}

#endif

#ifdef __cplusplus
    }
#endif

#endif /* !defined( POINTWISEXY_C_H ) */
