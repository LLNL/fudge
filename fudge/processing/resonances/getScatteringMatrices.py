# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy
from . import _getScatteringMatrices

# #### Functions & Classes ###########################################

def getScatteringMatrices(E, Eres, captureWidth, widths, penetrabilities):
	# .... Check arguments:
	for dat in (E,Eres,captureWidth,widths,penetrabilities):
		if type(dat) != numpy.ndarray: raise TypeError("All inputs must be numpy arrays!")
		if dat.dtype != numpy.float64: raise TypeError("All inputs must be float arrays")
	# check dimensions. Row or column is fine for incident energy:
	NE = len(E)
	dim = len(widths)
	nRes = len(Eres)
	if E.shape not in ((NE,),(NE,1)):
		raise TypeError("E has wrong shape: %s, should be %s" % (E.shape, (NE,)))
	if Eres.shape != (nRes,):
		raise TypeError("Eres has wrong shape: %s, should be %s" % (Eres.shape, (nRes,)))
	if captureWidth.shape != (nRes,):
		raise TypeError("captureWidth has wrong shape: %s, should be %s" % (captureWidth.shape, (nRes,)))
	if widths.shape != (dim,nRes):
		raise TypeError("widths have wrong shape: %s, should be %s" % (widths.shape, (dim,nRes)))
	if penetrabilities.shape[:2] != (dim,NE):
		raise TypeError("penetrabilities have wrong shape: %s, should be %s" % (penetrabilities.shape, (dim,NE)))
		
	# .... Call C extension function, and reshape results into 3d arrays:
	R,S = _getScatteringMatrices.getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
	R.shape = (NE, dim, dim)
	S.shape = (NE, dim, dim)
	return (R,S)

# #### Run test(s) ##################################################
#
# Note: these tests are replicated in tests/test_getScatteringMatrices.py (DAB 8/13/2013)
if __name__ == '__main__':

    #==== Tests ==============================

    def python_getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities):
        NE = len(E)
        dim = len(widths)
        dE = Eres - E[:,numpy.newaxis]
        DEN = 1/(dE**2 + captureWidth**2)
        captOverDEN = captureWidth*DEN
        dEoverDEN = dE*DEN
        del dE, DEN

        R = numpy.zeros((NE,dim,dim))
        S = numpy.zeros((NE,dim,dim))
        for i in range(dim):
            for j in range(dim):
                width_ij = widths[i]*widths[j]
                p_ij = penetrabilities[i]*penetrabilities[j]
                R[:,i,j] = p_ij.T * numpy.sum( width_ij * captOverDEN, axis=1)
                S[:,i,j] = p_ij.T * numpy.sum( width_ij * dEoverDEN, axis=1)
                # symmetrize:
                R[:,j,i] = R[:,i,j]
                S[:,j,i] = S[:,i,j]
        return R,S

    import unittest
    
    class Test_Driver( unittest.TestCase ):
    
        def test_getScatteringMatrices( self ):
            E = numpy.array([1e-5, 0.025, 0.1, 1.0, 10.0, 20.0])
            Eres = numpy.array([0.35, 2.4, 7.8])
            captureWidth = numpy.array([0.08, 1.2, 1.6])
            widths = numpy.array([[0.01,0.3,0.7],
                                [0,0,1.2e-4]])
            penetrabilities = numpy.array([[0.1,0.2,0.3,0.4,0.5,0.55],
                                        [0.05,0.15,0.25,0.35,0.45,0.55]])

            R,S = getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            R1,S1 = python_getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            self.assertTrue( numpy.allclose( R,R1 ) and numpy.allclose( S,S1 ) )

    unittest.main()

