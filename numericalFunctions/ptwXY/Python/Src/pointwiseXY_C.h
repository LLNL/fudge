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
Cannot get the PyCapsule_Import (i.e. python 2.7) way to work but the PyImport_ImportModuleEx with 2.6 and 2.7 so use it for now.
#if( PY_VERSION_HEX < 0x02070000 )
*/
#if 1

    PyObject *c_api_object, *module1, *module2, *fromList;

    if( ( fromList = Py_BuildValue( "(s)", "pointwiseXY_C" ) ) == NULL ) return( -1 );
    module1 = PyImport_ImportModuleEx( "numericalFunctions", NULL, NULL, fromList );
    Py_DECREF( fromList );
    if( module1 == NULL ) return -1;
    module2 = PyObject_GetAttrString( module1, "pointwiseXY_C" );
    if( module2 == NULL ) {
        Py_DECREF( module1 );
        return -1;
    }

    c_api_object = PyObject_GetAttrString( module2, "_C_API" );
    if( c_api_object == NULL ) {
        Py_DECREF( module1 );
        Py_DECREF( module2 );
        return( -1 );
    }
    if( PyCObject_Check( c_api_object ) ) Py_pointwiseXY_C_API = (void **) PyCObject_AsVoidPtr( c_api_object );

    Py_DECREF( c_api_object );
    Py_DECREF( module1 );
    Py_DECREF( module2 );
    return( 0 );

#else

    Py_pointwiseXY_C_API = (void **) PyCapsule_Import( "pointwiseXY_C._C_API", 0 );
    return( ( Py_pointwiseXY_C_API != NULL ) ? 0 : -1 );

#endif

}

#endif

#ifdef __cplusplus
    }
#endif

#endif /* !defined( POINTWISEXY_C_H ) */
