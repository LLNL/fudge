/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <Python.h>
#include <pyport.h>
#include "structmember.h"
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>

#include <nf_specialFunctions.h>

static PyObject *nf_amc_wigner_3j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_wigner_6j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_wigner_9j_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_racah_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_clebsh_gordan_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_z_coefficient_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_zbar_coefficient_C( PyObject *self, PyObject *args );
static PyObject *nf_amc_reduced_matrix_element_C( PyObject *self, PyObject *args );

/*
************************************************************
*/
static PyObject *nf_amc_wigner_3j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_3j( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_wigner_6j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_6j( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_wigner_9j_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6, j7, j8, j9;

    if( !PyArg_ParseTuple( args, "iiiiiiiii", &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_wigner_9j( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_racah_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, j4, j5, j6;

    if( !PyArg_ParseTuple( args, "iiiiii", &j1, &j2, &j3, &j4, &j5, &j6 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_racah( j1, j2, j3, j4, j5, j6 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_clebsh_gordan_C( PyObject *self, PyObject *args ) {

    int j1, j2, j3, m1, m2;

    if( !PyArg_ParseTuple( args, "iiiii", &j1, &j2, &m1, &m2, &j3 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_clebsh_gordan( j1, j2, m1, m2, j3 ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_zbar_coefficient_C( PyObject *self, PyObject *args ) {

    int j1, j2, l1, l2, ll, s;

    if( !PyArg_ParseTuple( args, "iiiiii", &l1, &j1, &l2, &j2, &s, &ll ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_zbar_coefficient( l1, j1, l2, j2, s, ll ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_z_coefficient_C( PyObject *self, PyObject *args ) {

    int j1, j2, l1, l2, ll, s;

    if( !PyArg_ParseTuple( args, "iiiiii", &l1, &j1, &l2, &j2, &s, &ll ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_z_coefficient( l1, j1, l2, j2, s, ll ) ) );
}
/*
************************************************************
*/
static PyObject *nf_amc_reduced_matrix_element_C( PyObject *self, PyObject *args ) {

    int lt, st, jt, l0, j0, l1, j1;

    if( !PyArg_ParseTuple( args, "iiiiiii", &lt, &st, &jt, &l0, &j0, &l1, &j1 ) ) return( NULL );

    return( Py_BuildValue( "d", nf_amc_reduced_matrix_element( lt, st, jt, l0, j0, l1, j1 ) ) );
}
/*
************************************************************
*/
static PyMethodDef nf_angularMomentumCouplingMethods[] = {
/*
Testing issues:
    wigner_3j( 6, 6, 4, 0, 0, 0 ) / wigner_3j( 4, 6, 6, 0, 0, 0 ) - 1 = -4.392597396929432e-12
*/
    { "wigner_3j", (PyCFunction) nf_amc_wigner_3j_C, METH_VARARGS,
        "wigner_3j( j1, j2, j3, j4, j5, j6 )\n\n" \
        "Wigner's 3J symbol (similar to Clebsh-Gordan)\n\n" \
        "    = / j1 j2 j3 \\\n" \
        "      \\ j4 j5 j6 /\n\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "wigner_6j", (PyCFunction) nf_amc_wigner_6j_C, METH_VARARGS,
        "wigner_6j( j1, j2, j3, j4, j5, j6 )\n\n" \
        "Wigner's 6J symbol (similar to Racah)\n\n" \
        "        = { j1 j2 j3 }\n" \
        "          { j4 j5 j6 }\n\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "wigner_9j", (PyCFunction) nf_amc_wigner_9j_C, METH_VARARGS,
        "wigner_9j( j1, j2, j3, j4, j5, j6, j7, j8, j9 )\n\n" \
        "Wigner's 9J symbol\n\n" \
        "      / j1 j2 j3 \\\n" \
        "    = | j4 j5 j6 |\n" \
        "      \\ j7 j8 j9 /\n" \
        "\nArguments are:\n" \
        "    j1, j2, j3, j4, j5, j6     all integers, 2x angular momenta\n" },
    { "racah", (PyCFunction) nf_amc_racah_C, METH_VARARGS,
        "racah( j1, j2, l2, l1, j3, l3 )\n\n" \
        "Racah coefficient\n\n" \
        "    = W(j1, j2, l2, l1 ; j3, l3)\n\n" \
        "    = (-1)^(j1+j2+l1+l2) * { j1 j2 j3 }\n" \
        "                           { l1 l2 l3 }\n\n" \
        "\nArguments are: \n" \
        "    j1, j2, l2, l1, j3, l3     all integers, 2x angular momenta\n\n"\
        "We use the Racah coefficient definition in Edmonds (A.R. Edmonds, 'Angular Momentum in Quantum Mechanics', Princeton (1980).\n"\
        "Note that the call signature of W(...) appears jumbled, but hey, that's the convention.\n"\
        "This convention is exactly that used by Blatt-Biedenharn (Rev. Mod. Phys. 24, 258 (1952)) too."
    },
    { "clebsh_gordan", (PyCFunction) nf_amc_clebsh_gordan_C, METH_VARARGS,
        "clebsh_gordan( j1, j2, j3, m1, m2)\n\n" \
        "Clebsh-Gordan coefficient\n\n" \
        "    = <j1,j2,m1,m2|j3,m1+m2>\n\n" \
        "    = (-)^(j1-j2+m1+m2) * sqrt(2*j3+1) * / j1 j2   j3   \\\n" \
        "                                         \\ m1 m2 -m1-m2 /\n\n" \
        "Note: Last value m3 is preset to m1+m2.  Any other value will evaluate to 0.0.\n\n"\
        "\nArguments are: \n" \
        "    j1, j2, j3     integers, 2x the angular momentum (so for j=1/2 use j=1)\n" \
        "    m1, m2, m3     integers, 2x the projection of the angular momentum onto the z axis\n" },
    { "z_coefficient", (PyCFunction) nf_amc_z_coefficient_C, METH_VARARGS,
        "z_coefficient( j1, j2, l1, l2, ll, s )\n\n" \
        "Biedenharn's Z-coefficient coefficient\n\n" \
        "    =  Z(l1  j1  l2  j2 | S L )\n\n" \
        "\nArguments are: \n" \
        "    j1, j2     integers, 2x the total angular momentum (so for j=1/2 use j=1)\n" \
        "    l1, l2     integers, 2x the orbital angular momentum\n" \
        "    ll         integer, 2x the orbital angular momentum that l1 and l2 couple up to,\n" \
        "    s          integer, 2x the projection of the ll onto the z axis\n" },
    { "zbar_coefficient", (PyCFunction) nf_amc_zbar_coefficient_C, METH_VARARGS,
        "zbar_coefficient( j1, j2, l1, l2, ll, s )\n\n" \
        "Lane & Thomas's Zbar-coefficient coefficient\n\n" \
        "    =  ZBar(l1  j1  l2  j2 | S L )\n\n" \
        "    = (-i)^( -l1 + l2 + ll ) * Z(l1  j1  l2  j2 | S L )\n\n"
        "From Lane & Thomas Rev. Mod. Phys. 30, 257-353 (1958).\n"
        "Note, Lane & Thomas define this because they did not like the different phase convention in Blatt & Biedenharn's Z coefficient.  They changed it to get better time-reversal behavior.\n"
        "Froehner uses Lane & Thomas convention as does T. Kawano.\n\n"
        "\nArguments are: \n" \
        "    j1, j2     integers, 2x the total angular momentum (so for j=1/2 use j=1)\n" \
        "    l1, l2     integers, 2x the orbital angular momentum\n" \
        "    ll         integer, 2x the orbital angular momentum that l1 and l2 couple up to,\n" \
        "    s          integer, 2x the projection of the ll onto the z axis\n" },
    { "reduced_matrix_element", (PyCFunction) nf_amc_reduced_matrix_element_C, METH_VARARGS,
        "reduced_matrix_element( lt, st, jt, l0, j0, l1, j1 )\n\n" \
        "Reduced Matrix Element for Tensor Operator\n\n" \
        "    = < l1j1 || T(YL,sigma_S)J || l0j0 >\n\n" \
        "From M.B.Johnson, L.W.Owen, G.R.Satchler, Phys. Rev. 142, 748 (1966)\n\n" \
        "Note: definition differs from JOS by the factor sqrt(2j1+1)\n\n" \
        "\nArguments are: \n" \
        "    lt     integer, 2x the orbital angular momentum of the tensor operator,\n" \
        "    st     integer, 2x the projection of lt onto the z axis,\n" \
        "    jt     integer, 2x the total angular momentum of the tensor operator,\n" \
        "    l0     integer, 2x the orbital angular momentum of the ket,\n" \
        "    j0     integer, 2x the total angular momentum of the ket,\n" \
        "    l1     integer, 2x the orbital angular momentum of the bra,\n" \
        "    j1     integer, 2x the total angular momentum of the bra.\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};

/*
************************************************************
*/

static char const doc[] = "A module that contains special math functions for angular momentum coupling in physics.";

#if PY_MAJOR_VERSION < 3

    PyMODINIT_FUNC initangularMomentumCoupling( void ) {

        Py_InitModule3( "angularMomentumCoupling", nf_angularMomentumCouplingMethods, doc );
    }

#else

    static struct PyModuleDef angularMomentumCoupling_CModule = {
        PyModuleDef_HEAD_INIT,
        "angularMomentumCoupling",      /* name of module */
        doc,                    /* module documentation, may be NULL */
        -1,                     /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        nf_angularMomentumCouplingMethods
    };

    PyMODINIT_FUNC PyInit_angularMomentumCoupling( void ) {

        return( PyModule_Create( &angularMomentumCoupling_CModule ) );
    }

#endif
