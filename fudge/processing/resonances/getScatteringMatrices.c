/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

static PyObject *getScatteringMatrices(PyObject *self, PyObject *args);

static int offset( int x, int y, int z, int xSize, int ySize );
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( char *s, ... );

/* .... Helper functions: allocate or free pointers to numpy memory: .... */
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
static double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin);
static double **ptrvector(long n);
static void free_Carrayptrs(double **v);

/* #### Globals #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef _getScatteringMatricesMethods[] = {
    {"getScatteringMatrices", getScatteringMatrices, METH_VARARGS},
    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

static char const doc[] = "A module that contains a function that returns energy-dependent matrices R and S for resonance reconstruction. \
    Where R and S are the symmetric and anti-symmetric scattering matrices, which must be inverted to calculate cross sections. \
        Returns tuple of two NEW NumPy arrays (R and S) \
        interface:  getScatteringMatrices(E, Eres, captureWidth, widths, penetrabilities) \
            E is a vector of incident energies  \
            Eres and captureWidth: resonance energies and capture widths \
            widths, penetrabilities: 2-d arrays of width/penetrability for other channels \
        returns tuple of flattened 3d arrays (R,S). The python wrapper must reshape R and S \
        into 3d arrays.";

/* ==== Initialize the C_test functions ====================== */
#if PY_MAJOR_VERSION < 3

    // Module name must be _C_arraytest in compile and linked 
    PyMODINIT_FUNC init_getScatteringMatrices( void )  {

        Py_InitModule3( "_getScatteringMatrices", _getScatteringMatricesMethods, doc );
        import_array();  // Need to import NumPy.  Called first after above line.
    }

#else

    static struct PyModuleDef _getScatteringMatrices_CModule = {
        PyModuleDef_HEAD_INIT,
        "_getScatteringMatrices",                   /* name of module */
        doc,                                        /* module documentation, may be NULL */
        -1,                                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        _getScatteringMatricesMethods
    };

    PyMODINIT_FUNC PyInit__getScatteringMatrices( void ) {

        import_array();  // import NumPy
        return( PyModule_Create( &_getScatteringMatrices_CModule ) );
    }

#endif


static PyObject *getScatteringMatrices(PyObject *self, PyObject *args)
{
    PyArrayObject *E, *Eres, *captureWidth;
    PyArrayObject *widths, *penetrabilities;
    PyArrayObject *R, *S;
    PyObject *returnTuple;
    double *c_E, *c_Eres, *c_captureWidth;
    double **c_widths, **c_penetrabilities;
    double *c_R, *c_S;
    int Ne, Nres, Nch, i,j,k1,k2;
    double dE, DEN, dEoverDEN, captOverDEN, wid_ij, p_ij;
    npy_intp dims[1]; // dimension of (flattened) R and S arrays

    /* Read in function arguments: */
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!", 
    	&PyArray_Type, &E, &PyArray_Type, &Eres, &PyArray_Type, &captureWidth,
    	&PyArray_Type, &widths, &PyArray_Type, &penetrabilities))  return NULL;

    /* Get dimensions from input */
    Ne=(int) PyArray_DIMS( E )[0];
    Nres=(int) PyArray_DIMS( Eres )[0];
    Nch=(int) PyArray_DIMS( widths )[0];

    dims[0] = Ne * Nch * Nch;

    /* allocate the scattering matrices R and S: */
    R=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    if (R==NULL) return NULL;
    S=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    if (S==NULL) {
        Py_DECREF(R);
        return NULL;
    }
    /* initialize to zero */
    PyArray_FILLWBYTE(R,0); PyArray_FILLWBYTE(S,0);

    /* Get C-style pointer arrays that map to same memory as the numpy arrays: */
    c_E=pyvector_to_Carrayptrs(E);
    c_Eres=pyvector_to_Carrayptrs(Eres);
    c_captureWidth=pyvector_to_Carrayptrs(captureWidth);
    c_widths=pymatrix_to_Carrayptrs(widths);
    c_penetrabilities=pymatrix_to_Carrayptrs(penetrabilities);
    c_R=pyvector_to_Carrayptrs(R);
    c_S=pyvector_to_Carrayptrs(S);
    
    /* Fill R and S matrices: */
    for (k1=0; k1<Nch; k1++) {     // loop over channels
      for (k2=k1; k2<Nch; k2++) {
        for (i=0; i<Ne; i++) {      // over incident energies
          for (j=0; j<Nres; j++) {  // over resonances
            dE = c_Eres[j]-c_E[i];
            DEN = 1.0 / (dE*dE + c_captureWidth[j]*c_captureWidth[j]);
            dEoverDEN = dE * DEN;
            captOverDEN = c_captureWidth[j] * DEN;
            wid_ij = c_widths[k1][j] * c_widths[k2][j];
            p_ij = c_penetrabilities[k1][i] * c_penetrabilities[k2][i];
            c_R[ offset(k1,k2,i, Nch,Nch) ] += wid_ij * p_ij * captOverDEN;
            c_S[ offset(k1,k2,i, Nch,Nch) ] += wid_ij * p_ij * dEoverDEN;
          }
        }
      }
    }
    // symmetrize
    for (k1=0;k1<Nch;k1++) {
      for (k2=0;k2<Nch;k2++) {
        for (i=0;i<Ne;i++) {
          c_R[ offset(k2,k1,i, Nch,Nch) ] = c_R[ offset(k1,k2,i, Nch,Nch) ];
          c_S[ offset(k2,k1,i, Nch,Nch) ] = c_S[ offset(k1,k2,i, Nch,Nch) ];
	}
      }
    }

    /* Free memory, close file and return */
    free_Carrayptrs(c_widths);
    free_Carrayptrs(c_penetrabilities);
    returnTuple = Py_BuildValue( "OO", R, S );
    Py_DECREF(R);
    Py_DECREF(S);
    return returnTuple;
}

/* ==== map 3-d array indices to 1d array ==== */
static int offset( int x, int y, int z, int xSize, int ySize ) {
    return ( z * xSize * ySize ) + ( y * xSize ) + x;
}

/* ==== Exception (only used if we run out of memory) ==== */
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( char *s, ... ) {

    va_list args;
    char Str[1024];

    va_start( args, s );
    vsnprintf( Str, sizeof( Str ), s, args );
    Str[sizeof( Str ) - 1] = 0;
    PyErr_SetString( PyExc_Exception, Str );
    va_end( args );
    return( NULL );
}


/* #### Numpy Utility functions ######################### */

/* ==== Create 1D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.             */
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
    return (double *) PyArray_DATA( arrayin );  /* pointer to arrayin data as double */
}

/* ==== Create 2D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory (for vector of pointers) is allocated        */
static double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    int i,n,m;

    n=(int) PyArray_DIMS( arrayin )[0];
    m=(int) PyArray_DIMS( arrayin )[1];
    c=ptrvector(n);
    a=(double *) PyArray_DATA( arrayin );  /* pointer to arrayin data as double */
    for ( i=0; i<n; i++)  {
    	c[i]=a+i*m;  }
    return c;
}

/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
static double **ptrvector(long n)  {
    double **v;
    v=(double **)malloc((size_t) (n*sizeof(double)));
    if (!v)   {
        pointwiseXY_C_SetPyErrorExceptionReturnNull(
                "In **ptrvector. Allocation of memory for double array failed.");
    	}
    return v;
}

/* ==== Free a double *vector (vec of pointers) ========================== */ 
static void free_Carrayptrs(double **v)  {
    free((char*) v);
}
