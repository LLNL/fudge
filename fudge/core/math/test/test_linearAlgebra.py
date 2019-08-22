#!/usr/bin/env python

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

from __future__ import print_function

from numpy import *
from numpy import linalg as LA
from fudge.core.math.linearAlgebra import *
import unittest



# --------- base class for all matrix tests -----------
class MatrixTests( unittest.TestCase ):

    def assertMatrixEqual( self, a, b ):
        """The numpy "==" operator is a universal function, so it operates on each element in the matrices.
        The result is another matrix whose elements are all a[i,j] == b[i,j].
        You need to call the .all() function to check that all elements evaluated to "True"
        Or, could just use array_equal()..."""
        #return self.assertTrue( ( a == b ).all() )
        return self.assertTrue( array_equal( a, b ) )
    
    def assertMatrixAlmostEqual( self, a, b, rtol=1e-05, atol=1e-08 ):
        return self.assertTrue( allclose( a, b, rtol=rtol, atol=atol ) )
    
    def assertMatrixTrue( self, a ):
        """You need to call the .all() function to check that all elements evaluated to "True" """
        return self.assertTrue( a.all() )



# --------- basic tests -----------
class BasicLinearAlgebraTests( MatrixTests ):
    """Very basic matrix tests"""
    
    def test_transpose( self ):
        """Transpose"""
        x = mat([[1, 0, 0],[1, 0, 0],[0, 0, 0]])
        self.assertMatrixEqual( x.T, mat([[1, 1, 0],[0, 0, 0],[0, 0, 0]]) )

    def test_inner_products(self):
        """Inner product test"""
        x = mat( [1,0,0])
        # * test
        self.assertEqual( x*x.T, [[1]] )
        # dot test
        self.assertEqual( dot(x,x.T), [[1]] )

    def test_outer_product(self):
        """Outer product test"""
        x = mat( [1,0,0])
        # outer product test
        if False: print( outer(x,x.T) == mat([[1, 0, 0],[0, 0, 0],[0, 0, 0]]) )
        self.assertMatrixEqual( outer(x,x.T), mat([[1, 0, 0],[0, 0, 0],[0, 0, 0]]) )

    def test_matrix_addition( self ):
        """Addition test"""
        x = mat( [[1,0,0],[0,0,0],[0,0,0]])
        y = mat( [[0,0,0],[0,1,0],[0,0,0]])
        z = zeros_like( x )
        z += x
        z += y
        self.assertMatrixEqual( x+y, mat([[1,0,0],[0,1,0],[0,0,0]]) )
        self.assertMatrixEqual( z, mat([[1,0,0],[0,1,0],[0,0,0]]) )

# --------- matrix composition tests -----------
class MatrixCompositionTests( MatrixTests ):
    """Very basic matrix tests"""
    
    def setUp( self ):
        self.a = mat([[ 1. ,  0.1,  0.1], [ 0.1,  1. ,  0.1], [ 0.1,  0.1,  1. ]])
        self.b = zeros( ( 3,6 ) )
        self.c = None
        self.d = mat([[1, 2], [3, 4]])
        self.e = mat([[5, 6], [7, 8]])
    
    def test_stackH( self ):
        """Testing horizontal stacking"""
        self.assertMatrixEqual( stackHorizontal( [ self.a, self.b, self.c ] ), 
            mat([[ 1.0,  0.1,  0.1,  0.,  0.,  0.,  0.,  0.,  0.],
                 [ 0.1,  1.0,  0.1,  0.,  0.,  0.,  0.,  0.,  0.],
                 [ 0.1,  0.1,  1.0,  0.,  0.,  0.,  0.,  0.,  0.]]) )

    def test_stackV( self ):
        """Testing vertical stacking"""
        self.assertMatrixEqual( stackVertical( [ self.d, self.e, self.c, self.d, self.e ] ), mat( [[ 1.,  2.], [ 3.,  4.], [ 5.,  6.], [ 7.,  8.], [ 1.,  2.], [ 3.,  4.], [ 5.,  6.], [ 7.,  8.]] ) )

    def test_stackD( self ):
        """Testing diagonal stacking"""
        self.assertMatrixEqual( stackDiagonal( [ self.d, self.e, self.a ] ), mat( [[ 1.,   2.,   0.,   0.,   0.,   0.,   0. ], [ 3.,   4.,   0.,   0.,   0.,   0.,   0. ], [ 0.,   0.,   5.,   6.,   0.,   0.,   0. ], [ 0.,   0.,   7.,   8.,   0.,   0.,   0. ], [ 0.,   0.,   0.,   0.,   1.,   0.1,  0.1], [ 0.,   0.,   0.,   0.,   0.1,  1.,   0.1], [ 0.,   0.,   0.,   0.,   0.1,  0.1,  1. ]]) )

# --------- CGLSQR tests -----------
class CGLSQRTests( MatrixTests ):

    def setUp( self ):
        self.answer = numpy.matrix([[ 1.34883721,-0.69767442, 0.34883721, 0.1, 42.0]])
        self.kernel = numpy.matrix( [ [ 1.0, 2.0, 3.0, 0.0, 0.0 ], [ 2.0, 3.0, 4.0, 0.0, 0.0 ], [ 4.5, 5.4, 2.0, 0.0, 0.0 ] ] ) # Note: lower two subspaces map to 0
        self.data = numpy.mat([ [ 1.1, 1.89, 3.05 ] ] )
        self.dataCov = numpy.mat([ [ 1.0, 0.1, 0.1 ], [ 0.1, 1.0, 0.1 ] , [ 0.1, 0.1, 1.0 ] ] )
        self.constraintVector = numpy.matrix( [[ 42.1 ]] )
        self.constraintMatrix = numpy.matrix( [[ 0.0, 0.0, 0.0, 1.0, 1.0 ]] )
        self.prior = numpy.matrix([[ 1.3, -0.7, 0.3, 0.11, 42.2 ]])
        self.priorCov = numpy.matrix(
       [[  0.07      ,   0.        ,   0.        ,   0.        ,   0.        ],
        [  0.        ,   0.5       ,   0.        ,   0.        ,   0.        ],
        [  0.        ,   0.        ,   0.4       ,   0.        ,   0.        ],
        [  0.        ,   0.        ,   0.        ,   0.1       ,   0.        ],
        [  0.        ,   0.        ,   0.        ,   0.        ,  20.        ]])
 
    def test_data_only( self ): 
        ans, ansCov, ansResid, ansChi2 = cglsqrSolve( self.data, dataUnc = None, dataCov = self.dataCov, kernel = self.kernel, prior = None, priorCov = None, constraintVector = None, constraintMatrix = None )
        if False:
            print('\ndata_only')
            print( 'model\n',    ans )
            print( 'data\n',     self.data )
            print( 'residual\n', self.data.T - self.kernel*ans )
            print( 'error\n',    self.answer.T - ans)
        self.assertMatrixAlmostEqual( ( self.answer - ans ).T, numpy.matrix( [[  0.68651163],[ -0.64302326],[  0.16651163],[  0.1       ],[ 42.        ]]) )

    def test_data_prior( self ): 
        ans, ansCov, ansResid, ansChi2 = cglsqrSolve( self.data, dataUnc = None, dataCov = self.dataCov, kernel = self.kernel, prior = self.prior, priorCov = self.priorCov, constraintVector = None, constraintMatrix = None )
        if False:
            print('\ndata_prior')
            print( 'model\n',    ans )
            print( 'data\n',     self.data )
            print( 'residual\n', self.data.T - self.kernel*ans )
            print( 'error\n',    self.answer.T - ans)
        self.assertMatrixAlmostEqual( ( self.answer - ans ).T, numpy.matrix([[ 0.04524624],[-0.05105358],[ 0.02455487],[-0.01      ],[-0.2       ]]) )

    def test_data_constraint( self ): 
        ans, ansCov, ansResid, ansChi2 = cglsqrSolve( self.data, dataUnc = None, dataCov = self.dataCov, kernel = self.kernel, prior = None, priorCov = None, constraintVector = self.constraintVector, constraintMatrix = self.constraintMatrix ) 
        if False:
            print('\ndata_constraint')
            print( 'model\n',       ans )
            print( 'data\n',        self.data )
            print( 'kernel\n',        self.kernel )
            #print( 'kernel*ans\n', self.kernel*ans.T )
            print( 'residual\n',    self.data.T - self.kernel*ans.T )
            print( 'error\n',       self.answer - ans)
            print( 'disobey constraint:\n',       self.constraintVector - self.constraintMatrix*ans.T)
        self.assertMatrixAlmostEqual( ( self.answer - ans ).T, numpy.matrix([[  0.68651163],[ -0.64302326],[  0.16651163],[-20.95      ],[ 20.95      ]]) )

    def test_data_constraint_prior( self ): 
        ans, ansCov, ansResid, ansChi2 = cglsqrSolve( self.data, dataUnc = None, dataCov = self.dataCov, kernel = self.kernel, prior = self.prior, priorCov = self.priorCov, constraintVector = self.constraintVector, constraintMatrix = self.constraintMatrix ) 
        if False:
            print('\ndata_constraint_prior')
            print( 'model\n',       ans )
            print( 'data\n',        self.data )
            print( 'residual\n',    self.data.T - self.kernel*ans.T )
            print( 'error\n',       self.answer - ans)
            print( 'disobey constraint:\n',       self.constraintVector - self.constraintMatrix*ans.T)
        self.assertMatrixAlmostEqual( ( self.answer - ans ).T, numpy.matrix(  [[ 0.04524624], [-0.05105358],[ 0.02455487],[-0.00895522],[ 0.00895522]] ) ) 

# ------------- eigendecomposition tests -----------------
class Eigendecomposition_base( MatrixTests ):

    def setUpMtx( self ): 
        raise Exception("Must be overridden by derived classes")

    def setUp( self ):   
        # The starting matrix, it is symmetric and real, but maybe not positive definite
        self.setUpMtx()
        self.ndim = self.A.shape[0]
        # The eigenvalue decomposition
        self.e, self.O = LA.eig( self.A )
        self.v = []
        for i in range( self.ndim ): self.v.append( self.O.T[i] )
        
class EigendecompositionTests:
    """ Define some tests that should be run for DERIVED classes only.
    This doesn't inherit from unittest.TestCase since the tests aren't meant to be run inside this class.
    
    Derived classes need to inherit from both this and from Eigendecomposition_base
    """

    def test_eigenmodes_are_really_eigenmodes( self ):
        """Check that A * v[i] = e[i] * v[i]"""
        for i in range( self.ndim ):
            self.assertMatrixAlmostEqual( dot( self.A, self.v[i].T ), self.e[i] * self.v[i].T )

    def test_eigenvector_orthogonality( self ):
        """Eigenvector orthogonality test"""
        for i in range( self.ndim ):
            for j in range( self.ndim ):
                if i == j: continue
                self.assertAlmostEqual( ( self.v[i] * self.v[j].T )[0,0], 0.0 )

    def test_eigenvector_normalization( self ):
        """Eigenvector normalization test:"""
        for i in range( self.ndim ): self.assertAlmostEqual( ( self.v[i] * self.v[i].T )[0,0], 1.0 )

    def test_eigenvector_orthonormality( self ):
        """Check that eigenvectors are orthonormal, better get back identity matrix here"""
        self.assertMatrixAlmostEqual( self.O.T * self.O, identity( self.ndim ) )
        self.assertMatrixAlmostEqual( self.O * self.O.T, identity( self.ndim ) )
        I = zeros_like( self.A )
        for i in range( self.ndim ): I += outer( self.v[i], self.v[i].T )
        self.assertMatrixAlmostEqual( I, identity( self.ndim ) )    
    
    def test_reconstruct_matrix_from_eigendecomposition( self ):
        """Try to reconstruct the matrix using the eigenvalue decomposition"""
        B = self.O * diag( self.e ) * self.O.T
        self.assertMatrixAlmostEqual( self.A, B ) 

    def test_construct_matrixinverse_from_eigendecomposition( self ):
        """Try to construct the matrix inverse using the eigenvalue decomposition"""
        B = self.O * diag( 1.0/self.e ) * self.O.T
        self.assertMatrixAlmostEqual( self.A * B, identity( self.ndim ) ) 
        self.assertMatrixAlmostEqual( B * self.A, identity( self.ndim ) ) 
 
    @unittest.expectedFailure
    def test_construct_matrixinverse_by_pruning( self ):
        """Try to construct the matrix inverse using the PRUNED eigenvalue decomposition"""
        with numpy.errstate(divide='ignore'):
            B = pruned_matrix_inverse( self.A )
        self.assertMatrixAlmostEqual( B * self.A, identity( self.ndim ), rtol=1e-05, atol=1e-08 ) 
        
    def test_construct_matrix_by_pruning( self ): 
        """Try to reconstruct the matix using the PRUNED eigenvalue decomposition"""
        B = pruned_matrix( self.A )
        self.assertMatrixAlmostEqual( self.A, B, rtol=1e-04, atol=1e-06 )
    
class EigendecompositionTests_allShouldFail( EigendecompositionTests ):
    """ For some input matrices, all tests are expected to fail. They should inherit from
    this class instead of EigendecompositionTests """

    @unittest.expectedFailure    
    def test_eigenmodes_are_really_eigenmodes( self ):
        super(self).test_eigenmodes_are_really_eigenmodes()

    @unittest.expectedFailure    
    def test_eigenvector_orthogonality( self ):
        super(self).test_eigenvector_orthogonality()

    @unittest.expectedFailure    
    def test_eigenvector_normalization( self ):
        super(self).test_eigenvector_normalization()

    @unittest.expectedFailure    
    def test_eigenvector_orthonormality( self ):
        super(self).test_eigenvector_orthonormality()
    
    @unittest.expectedFailure    
    def test_reconstruct_matrix_from_eigendecomposition( self ):
        super(self).test_reconstruct_matrix_from_eigendecomposition()

    @unittest.expectedFailure    
    def test_construct_matrixinverse_from_eigendecomposition( self ):
        super(self).test_construct_matrixinverse_from_eigendecomposition()
 
    @unittest.expectedFailure
    def test_construct_matrixinverse_by_pruning( self ):
        super(self).test_construct_matrixinverse_by_pruning() 
        
    @unittest.expectedFailure    
    def test_construct_matrix_by_pruning( self ): 
        super(self).test_construct_matrix_by_pruning() 

# --------- Specific Covariance Test Cases -----------

class BoringCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Boring, a diagonal matrix"""
    def setUpMtx( self ): self.A = mat( [[1.,0.,0.],[0.,2.,0.],[0.,0.,3.]] )

class SomeOffDiagonalCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Some off-diagonal-ness"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,0.],[1.,2.,0.],[0.,0.,3.]] )

class LotsOffDiagonalCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Lots of off-diagonal"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,1.],[1.,2.,1.],[1.,1.,3.]] )

class ZeroSubspaceCovarianceTests( Eigendecomposition_base, EigendecompositionTests_allShouldFail ):
    """Zeros? it still works!, but has negative eigenvalues"""
    def setUpMtx( self ): self.A = mat( [[0.,0.,1.],[0.,2.,1.],[1.,1.,3.]] ) 

class AllCorrelatedCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Pathological, all correlated, legal, but NOT INVERTIBLE"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,1.],[1.,1.,1.],[1.,1.,1.]] ) 

    @unittest.expectedFailure
    def test_construct_matrixinverse_from_eigendecomposition( self ):
        super(self).test_construct_matrixinverse_from_eigendecomposition()

    @unittest.expectedFailure
    def test_eigenvector_orthogonality( self ):
        super(self).test_eigenvector_orthogonality()

    @unittest.expectedFailure
    def test_eigenvector_orthonormality( self ):
        super(self).test_eigenvector_orthonormality()

class AlmostAllCorrelatedCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Pathological, all correlated, legal, but NOT INVERTIBLE"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,1.],[1.,1.,1.],[1.,1.,1.0001]] ) 

    @unittest.expectedFailure
    def test_construct_matrixinverse_from_eigendecomposition( self ):
        super(self).test_construct_matrixinverse_from_eigendecomposition()


class OnDiagonalBarelyPathologicalCovarianceTests( Eigendecomposition_base, EigendecompositionTests_allShouldFail ):
    """Pathological, all correlated, barely illegal and very much not invertible"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,1.],[1.,1.,1.],[1.,1.,0.9999]] ) 

class OffDiagonalBarelyPathologicalCovarianceTests( Eigendecomposition_base, EigendecompositionTests_allShouldFail ):
    """Pathological, all correlated, barely illegal and very much not invertible, and slightly off-diagonal"""
    def setUpMtx( self ): self.A = mat( [[1.,1.,1.],[1.,1.,1.00001],[1.,1.,0.9999]] )  

class IllegalCovarianceTest( Eigendecomposition_base, EigendecompositionTests_allShouldFail ): 
    """
    A bad matrix for testing:

        Eigensystem for A.  It is OK using numpy.linalg.eig, but not (apparently) with numpy.linalg.eigh, I checked against mathematica:
        >>> m = {{1., 2., 3.}, {1., 2., 1.}, {3., 2., 1.}}
        >>> {e, v} = Eigensystem[ m ]
        >>> {{5.23607, -2., 0.763932}, {{-0.647936, -0.400447, -0.647936}, {-0.707107, 1.08315*10^-16, 0.707107}, {0.465341, -0.752938, 0.465341}}}

        numpy.linalg.eig gives:
        >>> evals:
        >>>   [-2.14644241  0.78156156  5.36488085]
        >>> evecs:
        >>>   [[-0.64135318  0.52747554 -0.55716753]
        >>>   [-0.20230251 -0.81675359 -0.54035845]
        >>>   [ 0.74009445  0.23384422 -0.63053714]]
    
    However, note that this A is not a valid covariance as the off-diagonal elements are bigger than the on-diagonal ones.
    """
    def setUpMtx( self ): self.A = mat( [[1.,2.,3.],[1.,2.,1.],[3.,2.,1.]] )

    @unittest.expectedFailure
    def test_construct_matrix_by_pruning( self ):
        super(self).test_construct_matrix_by_pruning()

    @unittest.expectedFailure
    def test_construct_matrixinverse_from_eigendecomposition( self ):
        super(self).test_construct_matrixinverse_from_eigendecomposition()

class RealCovarianceTests( Eigendecomposition_base, EigendecompositionTests ):
    """Real covariance for testing, H(n,tot) cs cov."""
    def setUpMtx( self ):
        self.A = mat([
         [  8.77542100e-06,   1.38848800e-05,   2.79801400e-05,   3.88661200e-05,
            5.25228700e-05,   6.77747100e-05,   7.46361700e-05,   7.52854500e-05,
            7.06024400e-05,   6.12296700e-05,   4.77549800e-05,   3.07353700e-05,
            1.07019000e-05],
         [  1.38848800e-05,   2.21044000e-05,   4.46560200e-05,   6.20403900e-05,
            8.38345200e-05,   1.08155700e-04,   1.19094300e-04,   1.20119000e-04,
            1.12637300e-04,   9.76801600e-05,   7.61837900e-05,   4.90328600e-05,
            1.70763000e-05],
         [  2.79801400e-05,   4.46560200e-05,   9.03092300e-05,   1.25482100e-04,
            1.69564800e-04,   2.18748600e-04,   2.40870000e-04,   2.42936900e-04,
            2.27793700e-04,   1.97555100e-04,   1.54076600e-04,   9.91667500e-05,
            3.45412500e-05],
         [  3.88661200e-05,   6.20403900e-05,   1.25482100e-04,   1.74359700e-04,
            2.35628000e-04,   3.03993100e-04,   3.34742600e-04,   3.37618300e-04,
            3.16593500e-04,   2.74561300e-04,   2.14138000e-04,   1.37822000e-04,
            4.80018200e-05],
         [  5.25228700e-05,   8.38345200e-05,   1.69564800e-04,   2.35628000e-04,
            3.18449800e-04,   4.10888600e-04,   4.52472000e-04,   4.56392600e-04,
            4.27988300e-04,   3.71158600e-04,   2.89474100e-04,   1.86312900e-04,
            6.48883500e-05],
         [  6.77747100e-05,   1.08155700e-04,   2.18748600e-04,   3.03993100e-04,
            4.10888600e-04,   5.30240600e-04,   5.83960200e-04,   5.89030300e-04,
            5.52377700e-04,   4.79052800e-04,   3.73630300e-04,   2.40474700e-04,
            8.37559700e-05],
         [  7.46361700e-05,   1.19094300e-04,   2.40870000e-04,   3.34742600e-04,
            4.52472000e-04,   5.83960200e-04,   6.43158400e-04,   6.48730300e-04,
            6.08404400e-04,   5.27657300e-04,   4.11540800e-04,   2.64883200e-04,
            9.22628000e-05],
         [  7.52854500e-05,   1.20119000e-04,   2.42936900e-04,   3.37618300e-04,
            4.56392600e-04,   5.89030300e-04,   6.48730300e-04,   6.54432800e-04,
            6.13741100e-04,   5.32291300e-04,   4.15170400e-04,   2.67219300e-04,
            9.30924600e-05],
         [  7.06024400e-05,   1.12637300e-04,   2.27793700e-04,   3.16593500e-04,
            4.27988300e-04,   5.52377700e-04,   6.08404400e-04,   6.13741100e-04,
            5.75582300e-04,   4.99212100e-04,   3.89369700e-04,   2.50637600e-04,
            8.73319300e-05],
         [  6.12296700e-05,   9.76801600e-05,   1.97555100e-04,   2.74561300e-04,
            3.71158600e-04,   4.79052800e-04,   5.27657300e-04,   5.32291300e-04,
            4.99212100e-04,   4.32971200e-04,   3.37730200e-04,   2.17401500e-04,
            7.57803500e-05],
         [  4.77549800e-05,   7.61837900e-05,   1.54076600e-04,   2.14138000e-04,
            2.89474100e-04,   3.73630300e-04,   4.11540800e-04,   4.15170400e-04,
            3.89369700e-04,   3.37730200e-04,   2.63446700e-04,   1.69609800e-04,
            5.91586900e-05],
         [  3.07353700e-05,   4.90328600e-05,   9.91667500e-05,   1.37822000e-04,
            1.86312900e-04,   2.40474700e-04,   2.64883200e-04,   2.67219300e-04,
            2.50637600e-04,   2.17401500e-04,   1.69609800e-04,   1.09231900e-04,
            3.81553500e-05],
         [  1.07019000e-05,   1.70763000e-05,   3.45412500e-05,   4.80018200e-05,
            6.48883500e-05,   8.37559700e-05,   9.22628000e-05,   9.30924600e-05,
            8.73319300e-05,   7.57803500e-05,   5.91586900e-05,   3.81553500e-05,
            1.34232900e-05]] )

if __name__ == "__main__":
    unittest.main(verbosity=0)
    #In numpy, have around( obj[, decimals[, out] ), _round() too
    #numpy.testing.assert_approx_equal(actual, desired, significant=7, err_msg='', verbose=True)
