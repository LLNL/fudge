import unittest, math
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
        self.assertAlmostEqual( nf_amc.wigner_6j( 8, 4, 4, 4, 4, 8 ), 0.111677, 6 )
        self.assertAlmostEqual( nf_amc.wigner_6j( 4, 2, 2, 2, 2, 4 ), -0.223607, 6 )
        self.assertAlmostEqual( nf_amc.wigner_6j( 4, 4, 2, 2, 2, 4 ), 0.0745356 )



class Test_Wigner9jSymbol(unittest.TestCase):
    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.wigner_9j( 0, 0, 0, 0, 0, 0, 0, 0, 0 ), 1.0 )

    def test_nontrivial_values(self): 
        self.assertAlmostEqual( nf_amc.wigner_9j( 6, 2, 4, 2, 2, 2, 4, 2, 2 ), 0.03333333333333333 )



class Test_RacahCoefficient(unittest.TestCase):
    '''Note: Racah Coeffs are related to Wigner 6j symbols.  Some of these tests are redundant (by design)'''

    def test_all_zeroes(self): 
        self.assertAlmostEqual( nf_amc.racah( 0, 0, 0, 0, 0, 0 ), 1.0 )

    def test_nontrivial_values(self): 
        self.assertAlmostEqual( nf_amc.racah( 8, 4, 4, 4, 4, 8 ), 0.111677, 6 )
        self.assertAlmostEqual( nf_amc.racah( 4, 2, 2, 2, 2, 4 ), 0.223607, 6 )
        self.assertAlmostEqual( nf_amc.racah( 4, 4, 2, 2, 2, 4 ), 0.0745356 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 0 ), 0.16666666666550894 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 1 ), 0.03553345272565833 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 2 ), -0.12811768579269484 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 3 ), -0.040514369565350455 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 4 ), 0.06779350703073678 )
        self.assertAlmostEqual( nf_amc.racah( 2, 11, 2, 11, 9, 5 ), 0.022597835676191603 )



class Test_ZCoefficient(unittest.TestCase):

    def test_easy(self):
        '''
        For L=0, have analytic form from Froehner:
            ..math::
                \bar{Z}( \ell_1, J_1, \ell_2, J_2, s, 0 ) = (-)^{ J_1 + s }\sqrt{ 2 J_1 + 1 }\delta_{J_1 J_2}\delta_{\ell_1 \ell_2}
            '''
        for J1 in [0]:
            #print J1, nf_amc.z_coefficient( 0, int(2.*J1), 0, int(2.*J1), 1, 0 ), math.sqrt(2.*J1+1.)
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
            if L == 1:
                self.assertAlmostEqual( racah, -0.12811768579269484 ) 
                self.assertAlmostEqual( cgc, 0.0 )
            self.assertAlmostEqual( nf_amc.z_coefficient( int(2.*l1), int(2.*J1), int(2.*l2), int(2.*J2), int(2.*s), int(2.*L) ), math.sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*J1+1.)*(2.*J2+1.))*cgc*racah )



class Test_ReducedMatrixElement(unittest.TestCase):

    def test_something(self): 
        pass



if __name__=="__main__":
    unittest.main()
