from brownies.BNL.plot_evaluation.plotstyles import *

if __name__ == "__main__":
    import unittest


    class TestFunctions(unittest.TestCase):

        def setUp(self):
            self.userdefs = {
                "plotStyles": {
                    "crossSection": {
                        "title": "Craptastic plot for testing",
                        "xAxis": {
                            "min": 5.0,
                            "log": False},
                        "yAxis": {
                            "min": 0.0,
                            "max": 5.0,
                            "log": False}}},
                "evaluation": {
                    "junkFile": {
                        "filePattern": "junk.endf",
                        "legend": "This is crap!",
                        "lineColor": "Brown",
                        "errorColor": "Brown"}}}
            self.observable = 'crossSection'
            self.setkind = 'evaluation'
            self.setname = 'junk.endf'

        def test_getPlotSymbol(self):
            self.assertEqual(getPlotSymbol(6), 'p')

        def test_getPlotColor(self):
            self.assertEqual(getPlotColor(6), 'gray')

        def test_getThisPlotStyle(self):
            self.assertEqual(
                getThisPlotStyle(self.userdefs, self.observable),
                {
                    "title": "Craptastic plot for testing",
                    "showLegend": True,
                    "legendX": 0.05,
                    "legendY": 0.95,
                    "referenceFrame": None,
                    "xAxis": {
                        "unit": "MeV",
                        "min": 5.0,
                        "max": 20.0,
                        "log": False,
                        "label": "Incident energy"},
                    "yAxis": {
                        "unit": "b",
                        "min": 0.0,
                        "max": 5.0,
                        "log": False,
                        "label": "Cross section"}})

        def test_getThisSetStyle(self):
            self.assertEqual(
                getThisSetStyle(self.userdefs, self.setkind, self.setname),
                {
                    "filePattern": "junk.endf",
                    "legend": "This is crap!",
                    "lineStyle": "-",
                    "lineColor": "Brown",
                    "lineWidth": 3,
                    "errorColor": "Brown"})

        def test_readUserStyles(self):
            import tempfile
            styleFileHandle, styleFileName = tempfile.mkstemp()
            with open(styleFileName, mode='w') as jsonfile:
                jsonfile.write(json.dumps(self.userdefs, indent=2))
            self.assertEqual(readUserStyles(styleFileName), self.userdefs)
            os.remove(styleFileName)


    unittest.main()

