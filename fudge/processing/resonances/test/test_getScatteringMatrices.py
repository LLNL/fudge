# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import fudge.processing.resonances.getScatteringMatrices, unittest, numpy

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

if __name__ == '__main__':

    class Test_Driver( unittest.TestCase ):
    
        def test_getScatteringMatrices( self ):
            E = numpy.array([1e-5, 0.025, 0.1, 1.0, 10.0, 20.0])
            Eres = numpy.array([0.35, 2.4, 7.8])
            captureWidth = numpy.array([0.08, 1.2, 1.6])
            widths = numpy.array([[0.01,0.3,0.7],
                                [0,0,1.2e-4]])
            penetrabilities = numpy.array([[0.1,0.2,0.3,0.4,0.5,0.55],
                                        [0.05,0.15,0.25,0.35,0.45,0.55]])

            R,S = fudge.processing.resonances.getScatteringMatrices.getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            R1,S1 = python_getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            self.assertTrue( numpy.allclose( R,R1 ) and numpy.allclose( S,S1 ) )

    unittest.main()

