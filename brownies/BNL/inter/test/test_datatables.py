import unittest
from brownies.BNL.inter.datatables import *


class Test_DataTable(unittest.TestCase):

    def setUp(self):
        self.a = DataTable(
            data=[[1, 2, 3], [4, 5, 6], [11, 12, 13], [14, 15, 16], [21, 22.2342314, 23], [24, 25, 26],
                  [31, 32, 33], [34, 35, 36]],
            rows=['A', 'B', '1A', '1B', '2A', '2B', 'kinkajou', '3B'],
            columns=[
                xData.table.columnHeader(name='fred', unit='eV', index=0),
                xData.table.columnHeader(name='harry', unit='eV', index=1),
                xData.table.columnHeader(name='wow', unit='b', index=2)])

    def test_init(self):
        self.assertEqual(self.a[0][1], 2)
        self.assertEqual(self.a.columns[0].unit, 'eV')
        self.assertEqual(self.a.columns[0].name, 'fred')
        self.assertEqual(self.a.rows[0], 'A')
        self.assertEqual(self.a.rows[1], 'B')

    def test_str(self):
        expectedResult = """<table rows="8" columns="3">
  <columnHeaders>
    <column index="0" name="fred" unit="eV"/>
    <column index="1" name="harry" unit="eV"/>
    <column index="2" name="wow" unit="b"/></columnHeaders>
  <data>
   <!--    name | fred |      harry | wow  -->
   <!--         |   eV |         eV |   b  -->
               A      1            2     3
               B      4            5     6
              1A     11           12    13
              1B     14           15    16
              2A     21   22.2342314    23
              2B     24           25    26
        kinkajou     31           32    33
              3B     34           35    36</data></table>"""
        expectedOutlineResult = """<!--    name | fred |      harry | wow  -->
<!--         |   eV |         eV |   b  -->
           A      1            2     3
           B      4            5     6
          1A     11           12    13
 ...
          2B     24           25    26
    kinkajou     31           32    33
          3B     34           35    36"""
        self.maxDiff=None
        self.assertEqual(str(self.a), expectedResult)
        self.assertEqual('\n'.join(self.a.toStringList(outline=True)), expectedOutlineResult)

    def test_getElement(self):
        self.assertEqual(self.a.getElement('B', 'fred'), 4)


if __name__ == "__main__":

    unittest.main()
