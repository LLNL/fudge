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
#include <arrayobject.h>    /* for numpy c interface */
#include <math.h>
#include <stdarg.h>

static PyObject *getBreitWignerSums(PyObject *self, PyObject *args);

/* .... Helper functions: allocate pointer to numpy memory: .... */
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
static double square(double x) {return x*x;}

/* #### Globals #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef _getBreitWignerSumsMethods[] = {
    {"getBreitWignerSums", getBreitWignerSums, METH_VARARGS},
    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

static char const doc[] = "Module defines a helper function used for reconstructing Single-level or multi-level \
    Breit-Wigner resonances.";

/* ==== Initialize the C_test functions ====================== */
#if PY_MAJOR_VERSION < 3

    // Module name must be _C_arraytest in compile and linked 
    PyMODINIT_FUNC init_getBreitWignerSums( void )  {

        Py_InitModule3( "_getBreitWignerSums", _getBreitWignerSumsMethods, doc );
        import_array();  // Need to import NumPy.  Called first after above line.
    }

#else

    static struct PyModuleDef _getBreitWignerSums_CModule = {
        PyModuleDef_HEAD_INIT,
        "_getBreitWignerSums",                      /* name of module */
        doc,                                        /* module documentation, may be NULL */
        -1,                                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        _getBreitWignerSumsMethods
    };

    PyMODINIT_FUNC PyInit__getBreitWignerSums( void ) {

        import_array();  // import NumPy
        return( PyModule_Create( &_getBreitWignerSums_CModule ) );
    }

#endif


static PyObject *getBreitWignerSums(PyObject *self, PyObject *args)
{
    PyArrayObject *E, *Eres, *captureWidth, *elasticWidth, *fissionWidth;
    PyArrayObject *shiftFactorAtRes, *penetrationFactor, *shiftFactor, *phi;
    int approx;

    PyArrayObject *capture, *elastic, *fission;
    PyObject *returnTuple;
    double *c_E, *c_Eres, *c_captureWidth, *c_elasticWidth, *c_fissionWidth;
    double *c_shiftFactorAtRes, *c_penetrationFactor, *c_shiftFactor, *c_phi;
    double *c_capture, *c_elastic, *c_fission;
    npy_intp dims[1]; // dimension of return arrays

    // internal variables
    int Ne, Nres, i,j;
    double dE, totalWidth, denominator, commonFactor, elasticTerm1, elasticTerm2, sin2ps, sinps2;

    /* Read in function arguments: */
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!i",
    	&PyArray_Type, &E, &PyArray_Type, &Eres,
    	&PyArray_Type, &captureWidth, &PyArray_Type, &elasticWidth, &PyArray_Type, &fissionWidth,
    	&PyArray_Type, &shiftFactorAtRes, &PyArray_Type, &penetrationFactor, &PyArray_Type, &shiftFactor,
    	&PyArray_Type, &phi, &approx))    return NULL;

    /* Get dimensions from input */
    Ne=(int) PyArray_DIMS( E )[0];
    Nres=(int) PyArray_DIMS( Eres )[0];

    /* allocate the return arrays */
    dims[0] = Ne;
    capture=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    if (capture==NULL) return NULL;
    elastic=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    if (elastic==NULL) {
        Py_DECREF(capture);
        return NULL;
    }
    fission=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    if (fission==NULL) {
        Py_DECREF(capture);
        Py_DECREF(elastic);
        return NULL;
    }
    /* initialize to zero */
    PyArray_FILLWBYTE(capture,0); PyArray_FILLWBYTE(elastic,0); PyArray_FILLWBYTE(fission,0);

    /* Get C-style pointer arrays that map to same memory as the numpy arrays: */
    c_E=pyvector_to_Carrayptrs(E);
    c_Eres=pyvector_to_Carrayptrs(Eres);
    c_captureWidth=pyvector_to_Carrayptrs(captureWidth);
    c_elasticWidth=pyvector_to_Carrayptrs(elasticWidth);
    c_fissionWidth=pyvector_to_Carrayptrs(fissionWidth);
    c_shiftFactorAtRes=pyvector_to_Carrayptrs(shiftFactorAtRes);
    c_penetrationFactor=pyvector_to_Carrayptrs(penetrationFactor);
    c_shiftFactor=pyvector_to_Carrayptrs(shiftFactor);
    c_phi=pyvector_to_Carrayptrs(phi);

    c_capture=pyvector_to_Carrayptrs(capture);
    c_elastic=pyvector_to_Carrayptrs(elastic);
    c_fission=pyvector_to_Carrayptrs(fission);

    /* Actual computation */
    for (i=0; i<Ne; i++) {      // over incident energies
      elasticTerm1 = 0;
      elasticTerm2 = 0;
      sin2ps = 2 * sin(c_phi[i]) * cos(c_phi[i]);
      sinps2 = 2 * square( sin(c_phi[i]) );

      for (j=0; j<Nres; j++) {  // over resonances
        dE = c_E[i] - (c_Eres[j] + c_elasticWidth[j] * (c_shiftFactorAtRes[j] - c_shiftFactor[i]));
        totalWidth = c_penetrationFactor[i] * c_elasticWidth[j] + c_captureWidth[j] + c_fissionWidth[j];
        denominator = square(dE) + square(totalWidth) / 4;
        commonFactor = c_penetrationFactor[i] * c_elasticWidth[j] / denominator;
        c_capture[i] += commonFactor * c_captureWidth[j];
        c_fission[i] += commonFactor * c_fissionWidth[j];

        if (approx==1) {    // Single-level Breit-Wigner
          c_elastic[i] += commonFactor * (c_penetrationFactor[i] * c_elasticWidth[j]
                - sinps2 * totalWidth + 2 * sin2ps * (c_E[i] - c_Eres[j]));
        }
        else {              // Multi-level Breit-Wigner
          elasticTerm1 += totalWidth/2 * commonFactor;
          elasticTerm2 += dE * commonFactor;
          if (j==Nres-1)
            c_elastic[i] += square(sinps2 - elasticTerm1) + square(sin2ps + elasticTerm2);
        }
      }
    }

    /* Free memory, close file and return */
    returnTuple = Py_BuildValue( "OOO", capture, elastic, fission );
    Py_DECREF(capture);
    Py_DECREF(elastic);
    Py_DECREF(fission);
    return returnTuple;
}


/* #### Numpy Utility functions ######################### */

/* ==== Create 1D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.             */
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
    return (double *) PyArray_DATA( arrayin );  /* pointer to arrayin data as double */
}
