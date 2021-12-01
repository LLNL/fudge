/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <arrayobject.h>

#include "coulfg2.h"

static PyObject *getCoulombWavefunctions(PyObject *self, PyObject *args);

static PyMethodDef _getCoulombWavefunctionsMethods[] = {
    { "getCoulombWavefunctions", getCoulombWavefunctions, METH_VARARGS },
    { NULL, NULL } /* end of this structure */
};

static char const doc[] = "A module that contains a function that returns the Coulomb f and g values.";

#if PY_MAJOR_VERSION < 3

    PyMODINIT_FUNC init_getCoulombWavefunctions( void ) {

        Py_InitModule3( "_getCoulombWavefunctions", _getCoulombWavefunctionsMethods, doc );
        import_array(); // for numpy
    }

#else

    static struct PyModuleDef getCoulombWavefunctions_CModule = {
        PyModuleDef_HEAD_INIT,
        "_getCoulombWavefunctions",                 /* name of module */
        doc,                                        /* module documentation, may be NULL */
        -1,                                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        _getCoulombWavefunctionsMethods
    };

    PyMODINIT_FUNC PyInit__getCoulombWavefunctions( void ) {

        import_array(); // for numpy
        return( PyModule_Create( &getCoulombWavefunctions_CModule ) );
    }

#endif

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
    PyObject *FG_Py;

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
        if ( returnStatus != 0 ) {
            char errString[128];
            sprintf( errString, "Error in function Coulomb_FG! Error code: %d", returnStatus );
            PyErr_SetString(PyExc_Exception, errString);
            Py_DECREF(F); Py_DECREF(G); free(fc); free(gc);
            return NULL;
        }
        ((double*) PyArray_DATA( F ))[i] = fc[L];
        ((double*) PyArray_DATA( G ))[i] = gc[L];
    }

    free(fc);
    free(gc);
    FG_Py = Py_BuildValue("OO", F, G);
    Py_DECREF( F );
    Py_DECREF( G );

    return( FG_Py );

memErr:
        if( F != NULL ) { Py_DECREF( F ); }
        if( G != NULL ) { Py_DECREF( G ); }
        if( fc != NULL ) { free( fc ); }
        if( gc != NULL ) { free( gc ); }
        PyErr_NoMemory();
        return NULL;
}
