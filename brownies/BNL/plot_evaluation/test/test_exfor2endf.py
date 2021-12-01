import unittest
from brownies.BNL.plot_evaluation.exfor2endf import *

# ------------------------------------------------
# unit tests
# ------------------------------------------------


class TestQueries(unittest.TestCase):

    def test_query(self):
        s = session()
        out = s.query(EXFOR2ENDF).filter_by(exfor_reaction="(N,INL)PAR,DA")[0]
        self.assertEqual(out.projectile_za, 1)
        self.assertEqual(out.endf_mt, 51)
        self.assertEqual(out.endf_mf, 4)

    def test_get_exfor_reactions(self):
        out = get_exfor_reactions(1, 3, 16)
        self.assertEqual(out, ['(N,2N),SIG'])
        out = get_exfor_reactions(1, mt=16)
        self.assertEqual(out, ['(N,2N),SIG', '(N,2N),DA/DE'])
        out = get_exfor_reactions(1, mf=3, mt=18)
        self.assertEqual(out, ['(N,F),SIG,,RTE', '(N,F),SIG', 'NF', '(N,F),SIG,,AV', '(N,F),SIG,,MXW'])

    def test_get_endf_mf_mt(self):
        out = get_endf_mf_mt("(P,EL),DA,,RTH")
        self.assertEqual(out, [(4, 2)])


# ------------------------------------------------
# Main !!
# ------------------------------------------------

if __name__ == "__main__":
    #    try:
    #        import xmlrunner
    #        unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-results'))
    #    except ImportError:
    unittest.main()
#        print( )
