/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <math.h>
#include <Python.h>
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
    double *fc = NULL, *gc = NULL;
    int Ne, LMax, i;
    int returnStatus, m1;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!O!i",
                &PyArray_Type, &rho, &PyArray_Type, &eta, &L )) return NULL;

    Ne = (int)rho->dimensions[0];
    dims[0] = Ne;

    /* allocate F and G arrays: */
    if ( ( F= (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE) ) == NULL ) goto memErr;
    if ( ( G= (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE) ) == NULL ) goto memErr;

    LMax = ( L < 20 ? 20 : L );

    /* allocate working arrays */
    if ( ( fc = (double *) malloc( (LMax+1)*sizeof(double) ) ) == NULL ) goto memErr;
    if ( ( gc = (double *) malloc( (LMax+1)*sizeof(double) ) ) == NULL ) goto memErr;

    for (i=0; i<Ne; i++) {
        returnStatus = Coulomb_FG(((double*)rho->data)[i], ((double*)eta->data)[i], 
                0, LMax, fc, gc, &m1);
        // handle potential errors:
        if ( returnStatus < 0 ) {
            PyErr_SetString(PyExc_Exception, "Error in function Coulomb_FG!");
            Py_DECREF(F); Py_DECREF(G); free(fc); free(gc);
            return NULL;
        }
        ((double*)F->data)[i] = fc[L];
        ((double*)G->data)[i] = gc[L];
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

