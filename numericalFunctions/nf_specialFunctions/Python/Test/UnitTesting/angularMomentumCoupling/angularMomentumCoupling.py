# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest
import math

import os
import glob

import numericalFunctions

# Deployment via Makefiles places shared libraries in .../fudge/numericalFunctions/lib
if len(glob.glob(os.path.join(numericalFunctions.__path__[0], 'lib', '*angularMomentumCoupling*'))) > 0:
    from numericalFunctions.lib import angularMomentumCoupling as nf_amc

# Deployment via `pip install` places shared libraries in .../site-packages/numericalFunctions
else:
    from numericalFunctions import angularMomentumCoupling as nf_amc

class Test_ClebschGordanCoefficient(unittest.TestCase):

    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 0, 0, 0, 0, 0 ), 1.0 )
    
    def test_odd_sum_of_js(self):
        '''J = j1 + j2 + j3 == odd, then is zero'''
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 0, 2, 0, 0, 0 ), 0.0 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 6, 0, 0, 0, 0 ), 0.0 )

    def test_wikipedia_values(self):
        '''Tougher test value, from Wikipedia: http://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients'''
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 4, 1, 4, -1, 3 ), 0.8944271909999159 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 4, 1, 0, 1, 5 ), 0.7745966692414834 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 4, 4, 0, 0, 0 ), 0.4472135954999579 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 4, 0, 0, 0, 4 ), 1.0 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 8, 8, 0, 0, 0 ), 0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 8, 0, 0, 0, 8 ), 1.0 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 8, 8, 4, -4, 0 ), 0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 8, 8, 6, -6, 0 ), -0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 3, 2, 3, 0, 5 ), 0.6324555320336759 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 3, 2, 3, 0, 3 ), 0.7745966692414834 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 3, 2, 1, 2, 5 ), 0.7745966692414834 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 3, 2, 1, 2, 3 ), -0.6324555320336759 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 4, 6, 0, 0, 6 ), -0.516398, 6 )
        self.assertAlmostEqual( nf_amc.clebsh_gordan( 6, 4, 0, 0, 6 ), -0.516398, 6 )



class Test_Wigner3jSymbol(unittest.TestCase):

    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.wigner_3j( 0, 0, 0, 0, 0, 0 ), 1.0 )
    
    def test_odd_sum_of_js(self):
        '''J = j1 + j2 + j3 == odd, then is zero'''
        self.assertAlmostEqual( nf_amc.wigner_3j( 0, 2, 0, 0, 0, 0 ), 0.0 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 6, 0, 0, 0, 0, 0 ), 0.0 )

    def test_permutations(self): 
        '''
        Tougher test values, from Edmonds, pp. 50-51
        permutation test
        '''
        self.assertAlmostEqual( nf_amc.wigner_3j( 6, 4, 6, 0, 0, 0 ), 0.19518001458970666 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 6, 6, 4, 0, 0, 0 ), 0.19518001458970666 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 4, 6, 6, 0, 0, 0 ), 0.19518001458970666 )
        
    def test_other_nontrivial_values(self): 
        '''
        Tougher test values, from Edmonds, pp. 50-51
        '''
        self.assertAlmostEqual( nf_amc.wigner_3j( 8, 4, 8, 0, 0, 0 ), -0.16988239714587516 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 4, 4, 0, 0, 0, 0 ), 0.4472135954999579 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 4, 0, 4, 0, 0, 0 ), 0.4472135954999579 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 8, 8, 0, 0, 0, 0 ), 0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 8, 0, 8, 0, 0, 0 ), 0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 8, 8, 0, 4, -4, 0 ), 0.33333333333333333 )
        self.assertAlmostEqual( nf_amc.wigner_3j( 8, 8, 0, 6, -6, 0 ), -0.33333333333333333 )



class Test_Wigner6jSymbol(unittest.TestCase):

    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.wigner_6j( 0, 0, 0, 0, 0, 0 ), 1.0 )

    def test_nontrivial_values(self): 
        '''
        Values computed using "Anthony Stone's Wigner coefficient calculator"
        http://www-stone.ch.cam.ac.uk/wigner.shtml
        '''
        self.assertAlmostEqual( nf_amc.wigner_6j( 8, 4, 4, 4, 4, 8 ), 0.11167656571008167, 6 ) # (1/21)*sqrt(11/2)
        self.assertAlmostEqual( nf_amc.wigner_6j( 4, 2, 2, 2, 2, 4 ), -0.22360679774997896, 6 ) #-(1/2)*sqrt(1/5)
        self.assertAlmostEqual( nf_amc.wigner_6j( 4, 4, 2, 2, 2, 4 ), 0.07453559924999298 ) # (1/6)*sqrt(1/5)

    def test_nontrivial_values_sequence(self): 
        '''
        Values computed using "Anthony Stone's Wigner coefficient calculator"
        http://www-stone.ch.cam.ac.uk/wigner.shtml
        
        Values checked against A.Simon, J.H. VanderSluis, L.C. Biedenharn "Tables of the Racah Coefficients", ORNL-1679, March 26, 1954
        '''
        self.assertAlmostEqual( nf_amc.wigner_6j(  5,  8, 5,  8, 11, 6 ),  0.0161413111 ) #  (1/18) * sqrt(13/154)
        self.assertAlmostEqual( nf_amc.wigner_6j(  7,  8, 5,  8, 11, 6 ), -0.0435382494 ) # -(17/66)* sqrt(1/35)
        self.assertAlmostEqual( nf_amc.wigner_6j(  9,  8, 5,  8, 11, 6 ),  0.0623609564 ) #  (1/30) * sqrt(7/2)
        self.assertAlmostEqual( nf_amc.wigner_6j( 11,  8, 5,  8, 11, 6 ), -0.0341078118 ) # -(7/18) * sqrt(1/130)
        self.assertAlmostEqual( nf_amc.wigner_6j( 13,  8, 5,  8, 11, 6 ), -0.0544676990 ) # -(1/11) * sqrt(14/39)


class Test_Wigner9jSymbol(unittest.TestCase):
    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.wigner_9j( 0, 0, 0, 0, 0, 0, 0, 0, 0 ), 1.0 )

    def test_nontrivial_values(self): 
        self.assertAlmostEqual( nf_amc.wigner_9j( 6, 2, 4, 2, 2, 2, 4, 2, 2 ), 0.03333333333333333 )



class Test_RacahCoefficient(unittest.TestCase):
    '''Note: Racah Coeffs are related to Wigner 6j symbols.  Some of these tests are redundant (by design)'''

    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.racah( 0, 0, 0, 0, 0, 0 ), 1.0 )

    def test_nontrivial_values_from_6j_tests(self): 
        '''
        Values computed using "Anthony Stone's Wigner coefficient calculator"
        http://www-stone.ch.cam.ac.uk/wigner.shtml
        '''
        self.assertAlmostEqual( nf_amc.racah( 8, 4, 4, 4, 4, 8 ), 0.11167656571008167, 6 ) # (1/21)*sqrt(11/2)
        self.assertAlmostEqual( nf_amc.racah( 4, 2, 2, 2, 2, 4 ), 0.22360679774997896, 6 ) # (1/2)*sqrt(1/5)
        self.assertAlmostEqual( nf_amc.racah( 4, 4, 2, 2, 2, 4 ), 0.07453559924999298, 6 ) # (1/6)*sqrt(1/5)
        
    def test_nontrivial_values_sequence(self): 
        '''
        Values computed using "Anthony Stone's Wigner coefficient calculator"
        http://www-stone.ch.cam.ac.uk/wigner.shtml
        
        Values checked against A.Simon, J.H. Vander Sluis, L.C. Biedenharn "Tables of the Racah Coefficients", ORNL-1679, March 26, 1954
        '''
        self.assertAlmostEqual( nf_amc.racah(  5,  8, 11,  8, 5, 6 ),  0.0161413111 ) #  (1/18) * sqrt(13/154)
        self.assertAlmostEqual( nf_amc.racah(  7,  8, 11,  8, 5, 6 ),  0.0435382494 ) #  (17/66)* sqrt(1/35)
        self.assertAlmostEqual( nf_amc.racah(  9,  8, 11,  8, 5, 6 ),  0.0623609564 ) #  (1/30) * sqrt(7/2)
        self.assertAlmostEqual( nf_amc.racah( 11,  8, 11,  8, 5, 6 ),  0.0341078118 ) #  (7/18) * sqrt(1/130)
        self.assertAlmostEqual( nf_amc.racah( 13,  8, 11,  8, 5, 6 ), -0.0544676990 ) # -(1/11) * sqrt(14/39)



class Test_ZCoefficient(unittest.TestCase):

    def test_easy(self):
        '''
        For L=0, have analytic form from Froehner p. 49, eq. (149) (F. Froehner "Evaluation and Analysis of Nuclear Resonance Data", JEFF Report 18, OECD (2000)):
            ..math::
                \bar{Z}( \ell_1, J_1, \ell_2, J_2, s, 0 ) = (-)^{ J_1 + s }\sqrt{ 2 J_1 + 1 }\delta_{J_1 J_2}\delta_{\ell_1 \ell_2}
        Because we use spin 1/2 particles, only get non-zero results with 1/2-integer J1        
        '''
        for J1 in [0.5]:
            self.assertAlmostEqual( nf_amc.z_coefficient( 0, int(2.*J1), 0, int(2.*J1), 1, 0 ), math.sqrt(2.*J1+1.) )
         
    def test_something(self): 
        l1=1
        l2=1
        J1=5.5
        J2=5.5
        s=4.5 
        for L in range( 0, 11 ):
            cgc = nf_amc.clebsh_gordan( int(2*l1), int(2*l2), 0, 0, int(2*L) )
            racah = nf_amc.racah( int(2*l1), int(2*J1), int(2*l2), int(2*J2), int(2*s), int(2*L ) )
            if ( L - l1 + l2 ) % 8 == 0 : sign = 1.0
            else: sign = -1.0
            if L == 1:
                self.assertAlmostEqual( racah, -0.12811768579269484 ) 
                self.assertAlmostEqual( cgc, 0.0 )
            self.assertAlmostEqual( nf_amc.z_coefficient( int(2.*l1), int(2.*J1), int(2.*l2), int(2.*J2), int(2.*s), int(2.*L) ), sign*math.sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*J1+1.)*(2.*J2+1.))*cgc*racah )

    def test_ZCoefficient_Swave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual( math.sqrt(2.0),   nf_amc.z_coefficient(  0, 1,  0, 1, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(2.0),   nf_amc.z_coefficient(  2, 1,  2, 1, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(10.),   nf_amc.z_coefficient(  8, 9,  8, 9, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(10.),   nf_amc.z_coefficient( 10, 9, 10, 9, 1, 0 ) )
    
    def test_ZCoefficient_Pwave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual( -math.sqrt(40./3.),  nf_amc.z_coefficient(  8, 9, 6, 7, 1, 2 ) )
        self.assertAlmostEqual( -math.sqrt(40./3.),  nf_amc.z_coefficient( 10, 9, 8, 7, 1, 2 ) )


    def test_ZCoefficient_Dwave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual(   math.sqrt(8.0/7.0),  nf_amc.z_coefficient( 8, 7, 4, 5, 1, 4 ) )
        self.assertAlmostEqual(  -math.sqrt(8.0/7.0),  nf_amc.z_coefficient( 6, 7, 6, 5, 1, 4 ) )


class Test_ZBarCoefficient(unittest.TestCase):
    def test_ZBarCoefficient_Swave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual( math.sqrt(2.0),   nf_amc.zbar_coefficient(  0, 1,  0, 1, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(2.0),   nf_amc.zbar_coefficient(  2, 1,  2, 1, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(10.),   nf_amc.zbar_coefficient(  8, 9,  8, 9, 1, 0 ) )
        self.assertAlmostEqual( math.sqrt(10.),   nf_amc.zbar_coefficient( 10, 9, 10, 9, 1, 0 ) )
    
    def test_ZBarCoefficient_Pwave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual( -math.sqrt(40./3.),  nf_amc.zbar_coefficient(  8, 9, 6, 7, 1, 2 ) )
        self.assertAlmostEqual( -math.sqrt(40./3.),  nf_amc.zbar_coefficient( 10, 9, 8, 7, 1, 2 ) )


    def test_ZBarCoefficient_Dwave(self):
        """
        Values checked against L.C. Biedenharn, "Revised Z Tables of the Racah Coefficients", ORNL-1501, May 28, 1953.
        """
        self.assertAlmostEqual(   math.sqrt(8.0/7.0),  nf_amc.zbar_coefficient( 8, 7, 4, 5, 1, 4 ) )
        self.assertAlmostEqual(   math.sqrt(8.0/7.0),  nf_amc.zbar_coefficient( 6, 7, 6, 5, 1, 4 ) )


class Test_ReducedMatrixElement(unittest.TestCase):

    def test_something(self): 
        pass



if __name__=="__main__":
    unittest.main()
