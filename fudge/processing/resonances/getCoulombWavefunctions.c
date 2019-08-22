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

