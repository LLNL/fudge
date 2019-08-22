import xData.table
import pqu.PQU as PQUModule

class DataTable( xData.table.table ):

    def __init__( self, columns = None, rows = None, data = None ):

        xData.table.table.__init__(self,columns=columns,data=data)
        self.rows=rows
        if self.rows and len(self.rows)!=self.nRows:
            raise ValueError("Data is the wrong shape for a table with %i rows!" % self.nRows)

    def __str__(self): return '\n'.join(self.toXMLList())

    def hasMatchingRow(self,key,value):
        try:
            icol=[x.name for x in self.columns].index(key)
        except ValueError: return False
        for row in self.data:
            if row[icol]==value:return True
        return False

    def getMatchingRow(self,key,value):
        icol=[x.name for x in self.columns].index(key)
        for i,row in enumerate(self.data):
            if row[icol] == value: return i,self.rows[i],row
        return (None,None,None)

    def toStringList(self, indent='',**kwargs):

        # If not row labels, revert to base class behavior
        if not self.rows: return xData.table.table.toStringList(self,indent=indent,kwargs=kwargs)

        addHeader = kwargs.get( 'addHeader', True )
        addHeaderUnits = kwargs.get( 'addHeaderUnits', True )
        outline = kwargs.get( 'outline', False )
        columnWidths = [0] * (self.nColumns+1)
        for col in range(self.nColumns):
            columnDat = [row[col] for row in self.data if not isinstance(row[col], xData.table.blank)]
            asStrings = map( PQUModule.toShortestString, columnDat )
            columnWidths[col+1] = max( map(len, asStrings) )
        columnWidths[0]=max(map(len,self.rows))

        if addHeader:
            """ put column labels at the top of the table """
            names = ['name']+[col.name for col in self.columns]
            lengths = [len(name) for name in names]
            if any( [' width' in name for name in names] ):
                # special treatment for RML widths: split up onto two lines
                nameL = [[],[]]
                for name in names:
                    if ' width' in name: l,r = name.split(' width'); r = ' width'+r+' '
                    else: l,r = name, ''
                    nameL[0].append(l); nameL[1].append(r)
                names = nameL
                lengths = [len(name) for name in names[0]]
            else: names = [names]
            columnWidths = [max(columnWidths[i],lengths[i]) for i in range(len(columnWidths))]
            header = ['%s<!--' % (' '*(len(indent)-1)) + ' | '.join([('%%%is'%columnWidths[i]) % nameList[i] for i in range(len(columnWidths))]) + '  -->' for nameList in names]
            xml = header

        if addHeaderUnits:
            """ put column units at the top of the table """
            units = ['']+[str(col.unit).replace('None','') for col in self.columns]
            header = ['%s<!--' % (' '*(len(indent)-1)) + ' | '.join([('%s'%u).rjust(columnWidths[i]) for i,u in enumerate(units) ] ) + '  -->']
            xml += header

        template = ['%s' % (indent + ' ')] + ['%%%is' % l for l in columnWidths]

        def toString(val):
            if isinstance(val, xData.table.blank): return str(val)
            if isinstance(val, str): return val
            return PQUModule.toShortestString(val)

        if outline:
            xml += [('   '.join(template) % tuple( map(toString, [self.rows[irow]]+dataRow))).rstrip() for irow,dataRow in enumerate(self.data[:3])]
            xml += ['%s ...' % indent]
            xml += [('   '.join(template) % tuple( map(toString, [self.rows[-3+irow]]+dataRow))).rstrip() for irow,dataRow in enumerate(self.data[-3:])]
        else:
            xml += [('   '.join(template) % tuple( map(toString, [self.rows[irow]]+dataRow))).rstrip() for irow,dataRow in enumerate(self.data)]
        return xml

    def getElement(self,rowName,colName):
        column = [a for a in self.columns if a.name==colName]
        if not column: return None
        if len(column) > 1: raise ValueError("Column named '%s' is not unique!" % colName)
        icol= self.columns.index( column[0] )
        row = [a for a in self.rows if a==rowName]
        if not row: return None
        if len(row) > 1: raise ValueError("Row named '%s' is not unique!" % rowName)
        irow = self.rows.index( row[0] )
        return self.data[irow][icol]


if __name__ == "__main__":
    import unittest

    class Test_DataTable( unittest.TestCase ):

        def setUp(self):
            self.a = DataTable(data=[[1,2,3],[4,5,6],[11,12,13],[14,15,16],[21,22.2342314,23],[24,25,26],[31,32,33],[34,35,36]],
                               rows=['A','B','1A','1B','2A','2B','kinkajou','3B'],
                               columns=[
                                    xData.table.columnHeader(name='fred',unit='eV',index=0),
                                    xData.table.columnHeader(name='harry',unit='eV',index=1),
                                    xData.table.columnHeader(name='wow',unit='b',index=2)])

        def test_init(self):
            self.assertEqual( self.a[0][1], 2 )
            self.assertEqual( self.a.columns[0].unit, 'eV' )
            self.assertEqual( self.a.columns[0].name, 'fred' )
            self.assertEqual( self.a.rows[0], 'A' )
            self.assertEqual( self.a.rows[1], 'B' )

        def test_str(self):
            expectedResult="""<table rows="8" columns="3">
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
            self.assertEqual(str(self.a),expectedResult)
            self.assertEqual('\n'.join(self.a.toStringList(outline=True)), expectedOutlineResult)

        def test_getElement(self):
            self.assertEqual(self.a.getElement('B','fred'),4)

#        def test_to_json(self):
#            self.assertEqual(self.a.to_json(), '{"["fred","eV"]":{"A":1,"B":4},"["harry","eV"]":{"A":2,"B":5},"["wow","b"]":{"A":3,"B":6}}')
#
#        def test_to_csv(self):
#            self.assertEqual(self.a.to_csv(), 'name,fred,harry,wow\nunit,eV,eV,b\n,,,\nA,1,2,3\nB,4,5,6\n')


    unittest.main()