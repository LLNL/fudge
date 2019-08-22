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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <arrayobject.h>

#include "coulfg2.h"

static PyObject *getCoulombWavefunctions(PyObject *self, PyObject *args);

/* two methods required when setting up a python extension module: */
static PyMethodDef _getCoulombWavefunctionsMethods[] = {
    {"getCoulombWavefunctions", getCoulombWavefunctions, METH_VARARGS},
    {NULL,NULL} /* end of this structure */
};

DL_EXPORT( void ) init_getCoulombWavefunctions( void ) {
    (void) Py_InitModule("_getCoulombWavefunctions", _getCoulombWavefunctionsMethods);
    import_array(); // for numpy
}

static PyObject *getCoulombWavefunctions(PyObject *self, PyObject *args)
{
    // inputs: rho, eta and L.  outputs: F, G
    PyArrayObject *rho, *eta, *F = NULL, *G = NULL;
    int L;
    // internal variables:
    double *rho_c = NULL, *eta_c = NULL, *fc = NULL, *gc = NULL;
    int Ne, LMax, i;
    int returnStatus, m1;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!O!i",
                &PyArray_Type, &rho, &PyArray_Type, &eta, &L )) return NULL;

    Ne = (int) PyArray_DIMS( rho )[0];
    dims[0] = Ne;

    /* get pointers to input data arrays 'rho' and 'eta': */
    rho_c = (double*) PyArray_DATA( rho );
    eta_c = (double*) PyArray_DATA( eta );

    /* allocate F and G arrays: */
    if ( ( F= (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE) ) == NULL ) goto memErr;
    if ( ( G= (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE) ) == NULL ) goto memErr;

    LMax = ( L < 20 ? 20 : L );

    /* allocate working arrays */
    if ( ( fc = (double *) malloc( (LMax+1)*sizeof(double) ) ) == NULL ) goto memErr;
    if ( ( gc = (double *) malloc( (LMax+1)*sizeof(double) ) ) == NULL ) goto memErr;

    for (i=0; i<Ne; i++) {
        returnStatus = Coulomb_FG( rho_c[i], eta_c[i], 0, LMax, fc, gc, &m1);
        // handle potential errors:
        if ( returnStatus < 0 ) {
            PyErr_SetString(PyExc_Exception, "Error in function Coulomb_FG!");
            Py_DECREF(F); Py_DECREF(G); free(fc); free(gc);
            return NULL;
        }
        ((double*) PyArray_DATA( F ))[i] = fc[L];
        ((double*) PyArray_DATA( G ))[i] = gc[L];
    }

    free(fc); free(gc);
    return Py_BuildValue("(OO)", F, G);

memErr:
        if( F != NULL ) { Py_DECREF( F ); }
        if( G != NULL ) { Py_DECREF( G ); }
        if( fc != NULL ) { free( fc ); }
        if( gc != NULL ) { free( gc ); }
        PyErr_NoMemory();
        return NULL;
}

