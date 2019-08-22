import pandas, numpy
from pqu.PQU import PhysicalQuantityWithUncertainty as PQU

'''
A wrapper around pandas DataFrame class needs factory function creators.  Here are the ones in the pandas module that
we should target:

    pandas.read_clipboard  pandas.read_fwf        pandas.read_html       pandas.read_pickle     pandas.read_sql_table
    pandas.read_csv        pandas.read_gbq        pandas.read_json       pandas.read_sql        pandas.read_stata
    pandas.read_excel      pandas.read_hdf        pandas.read_msgpack    pandas.read_sql_query  pandas.read_table
'''

class DataTable( pandas.DataFrame ):

    """
    Simple wrapper around pandas.DataFrame that allows units on a column.

    Add these functions:

        a.convert_to_unit( column, newUnit )
        a.to_XML
    """

    def __init__( self, data=[], index=[], columns=[], units=None, **kw ):
        """
        Two-dimensional size-mutable, potentially heterogeneous tabular data structure with labeled axes (rows and columns).
        Arithmetic operations align on both row and column labels.
        Can be thought of as a dict-like container for Series objects.

        The primary pandas data structure

        Our version adds the ability to associate units with a column

        :param data: numpy ndarray (structured or homogeneous), dict, or DataFrame;
                     Dict can contain Series, arrays, constants, or list-like objects
        :param index: Index or array-like
                      Index to use for resulting frame. Will default to np.arange(n) if no indexing information part of input data and no index provided
        :param columns: Index or array-like
                        Column labels to use for resulting frame. Will default to np.arange(n) if no column labels are provided
        :param units: Index or array-like
                      Units to associate with the column in the resulting frame.  Will default to None if not provided.
        :param dtype : dtype, default None
                       Data type to force, otherwise infer
        :param copy : boolean, default False
                      Copy data from inputs. Only affects DataFrame / 2d ndarray input
        :return:
        """
        # Check that units match up with column labels
        if units is not None:
            if len(units) != len(columns):
                raise ValueError( "Size of units list (%i) != number of columns (%i)" % (len(units),len(columns)) )
        else:
            units = [ None for x in columns ]
        # Initialize the pandas.DataFrame instance
        pandas.DataFrame.__init__(
            self,
            data=data,
            index=index,
            columns=pandas.MultiIndex.from_tuples( list(zip(columns,units)), names=['name','unit'] ),
            **kw)

    def __getitem__(self,key):
        levelOne = pandas.DataFrame.__getitem__(self,key)
        return levelOne[levelOne.keys()[0]]

    def get_unit(self,key):
        return pandas.DataFrame.__getitem__(self,key).keys()[0]

    def to_XML(self):
        raise NotImplementedError("write me")

    def convert_to_unit(self, column, unit):
        """
        Converts a column from current unit to a new unit
        :param column: column name
        :param unit: new column unit
        :return: None
        """
        # Which column are we modifying?
        oldUnit = self.get_unit(column)
        iCol = self.columns.get_loc(column).start

        # Convert the column
        convertUnits=numpy.frompyfunc(lambda x: float(PQU(x,oldUnit).inUnitsOf(unit)),1,1) # make a verctorized unit conversion widget
        self[column]=convertUnits(self[column]) # converts data in a column

        # Rebuild the column index
        newUnits=[]
        newCols=[]
        for i in range(self.columns.size):
            newCols.append(self.columns[i][0])
            if i == iCol: newUnits.append(unit)
            else: newUnits.append(self.columns[i][1])
        self.columns=pandas.MultiIndex.from_tuples( list(zip(newCols,newUnits)), names=['name','unit'] )

if __name__ == "__main__":
    import unittest

    class Test_DataTable( unittest.TestCase ):

        def setUp(self):
            self.a = DataTable(data=[[1,2,3],[4,5,6]], index=['A','B'], columns=['fred','harry','wow'], units=['eV','eV','b'])

        def test_init(self):
            self.assertEqual( self.a['wow']['B'], 6 )
            self.assertEqual( self.a.columns[0][1], 'eV' )
            self.assertEqual( self.a.columns[0][0], 'fred' )

        def test_str(self):
            self.assertEqual(str(self.a),'name fred harry wow\nunit   eV    eV   b\nA       1     2   3\nB       4     5   6')

        def test_to_json(self):
            self.assertEqual(self.a.to_json(), '{"["fred","eV"]":{"A":1,"B":4},"["harry","eV"]":{"A":2,"B":5},"["wow","b"]":{"A":3,"B":6}}')

        def test_to_csv(self):
            self.assertEqual(self.a.to_csv(), 'name,fred,harry,wow\nunit,eV,eV,b\n,,,\nA,1,2,3\nB,4,5,6\n')

        def test_get_unit(self):
            self.assertEqual(self.a.get_unit('fred'),'eV')
            self.assertEqual(self.a.get_unit('wow'),'b')

        def test_column_features(self):
            self.assertEqual(str(self.a.fred),'A    1\nB    4\nName: eV, dtype: int64')

        def test_add_remove_columns(self):
            self.a.insert(self.a.columns.size,('frodo','amu'),value=[7.0,5.0])
            self.assertEqual(str(self.a),'name fred harry wow frodo\nunit   eV    eV   b   amu\nA       1     2   3     7\nB       4     5   6     5')
            self.a.pop('wow')
            self.assertEqual(str(self.a),'name fred harry frodo\nunit   eV    eV   amu\nA       1     2     7\nB       4     5     5')

        def test_convert_to_unit(self):
            self.a.convert_to_unit('fred','keV')
            self.assertEqual(str(self.a),'name   fred harry wow\nunit    keV    eV   b\nA     0.001     2   3\nB     0.004     5   6')

    unittest.main()