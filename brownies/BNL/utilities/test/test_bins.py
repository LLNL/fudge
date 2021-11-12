import unittest
from brownies.BNL.utilities.bins import *


class Test_Utilities(unittest.TestCase):

    def test_equal_lethargy_bins(self):
        self.maxDiff = None
        answer = [1.0e-05, 0.00037606030930863936, 0.014142135623730951, 0.53182958969449889, 20.0]
        calc = equal_lethargy_bins(5, domainMax=20.0)
        for i in range(len(calc)):
            self.assertAlmostEqual(answer[i], calc[i], 7)


if __name__=="__main__":
    unittest.main()
