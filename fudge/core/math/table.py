# <<BEGIN-copyright>>
# <<END-copyright>>

__metaclass__ = type

from fudge.core.ancestry import ancestry

tableToken = 'table'

class table(ancestry):
    def __init__( self, columns = None, data = None ):
        ancestry.__init__(self, tableToken, None)
        self.columns = columns or []
        self.data = data or []
        self.nColumns, self.nRows = len(self.columns), len(self.data)
        if not all( [len(d)==self.nColumns for d in self.data] ):
            raise Exception, ("Data is the wrong shape for a table with %i columns!" % self.nColumns)

    def __len__( self ):
        return len(self.data)

    def __getitem__( self, indices ):
        if type(indices) is int: return self.data[ indices ]
        if len(indices)==2:
            i,j = indices
            if type(i) is slice:
                return [d[j] for d in self.data[i]]
            return self.data[i][j]
        raise IndexError("invalid index")

    def addRow( self, dataRow ):
        if not len(dataRow) == self.nColumns:
            raise Exception, ("New row has %i columns, should have %i!" % (len(dataRow),self.nColumns))
        self.data.append( dataRow )
        self.nRows += 1

    def addColumn( self, columnHeader, index=None ):
        """ add another column, either at 'index' or at the end of the table """
        if index:
            self.columns.insert(index, columnHeader)
            [row.insert(index, 0) for row in self.data]
        else:
            self.columns.append( columnHeader )
            [row.append(0) for row in self.data]
        self.nColumns += 1
        for idx, col in enumerate(self.columns): col.index = idx

    def getColumn( self, columnName, units=None ):
        """ get data from one column, identified by the column 'name' attribute.
        Convert results to unit if units are specified """
        column = [a for a in self.columns if a.name==columnName]
        if not column: return None
        if len(column) > 1: raise Exception("Column named '%s' is not unique!" % columnName)
        index = self.columns.index( column[0] )
        if units:
            from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
            cf = ( PQU(1, column[0].units) / PQU(1, units) ).getValue()
            return [cf * v for v in self[:,index]]
        return self[:,index]

    def removeColumn( self, columnName ):
        """ remove one column from the table """
        column = [a for a in self.columns if a.name==columnName]
        if not column: raise ValueError("column '%s' isn't present in the table!" % columnName)
        if len(column) > 1: raise Exception("Column named '%s' is not unique!" % columnName)
        index = self.columns.index( column[0] )
        self.columns.pop( index )
        [row.pop(index) for row in self.data]
        self.nColumns -= 1
        for idx, col in enumerate(self.columns): col.index = idx

    def toXMLList( self, indent='', prettyPrint=True, addHeader=True ):
        xml = ['%s<table rows="%i" columns="%i">' % (indent,self.nRows,self.nColumns),
                '%s<columnHeaders>' % (indent+'  ')]
        for column in self.columns: xml += column.toXMLList( indent+'    ' )
        xml[-1] += '</columnHeaders>'
        if not self.data:
            xml.append( '%s<data/></table>' % (indent+'  ') )
            return xml

        xml.append( '%s<data>' % (indent+'  ') )
        columnWidths = [0] * self.nColumns
        for col in range(self.nColumns):
            columnDat = [row[col] for row in self.data if not isinstance(row[col], blank)]
            asStrings = ['%s' % val for val in columnDat]
            if map(float, asStrings) == columnDat:
                columnWidths[col] = max( [len(val) for val in asStrings] )
            else:
                raise Exception("need higher precision when writing to table!")
                # '%s' doesn't preserve enough precision. This doesn't come up when translating ENDF,
                # but will have to be dealt with eventually. Not sure if the following works correctly:
                while precision>1:
                    format = ('%%-.%ie ' % precision) * len(columnDat)
                    string = format % tuple(columnDat)
                    if columnDat == map(float, string.split()): precision -= 1
                    else:
                        precision += 1
                        break
                template.append( '%% -.%ie' % precision )
        if addHeader:
            """ put column labels at the top of the table """
            names = [col.name for col in self.columns]
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
            columnWidths = [max(columnWidths[i],lengths[i]) for i in range(self.nColumns)]

            header = ['%s<!--' % (' '*(len(indent)-1)) + ' | '.join([('%%%is'%columnWidths[i]) % nameList[i]
                for i in range(self.nColumns)]) + '  -->'   for nameList in names]
            xml += header

        template = ['%s' % (indent + ' ')] + ['%%%is' % l for l in columnWidths]
        xml += [('   '.join(template) % tuple(dataRow)).rstrip() for dataRow in self.data]
        xml[-1] += '</data></table>'
        return xml

class columnHeader:
    """ defines one column in a table """
    def __init__( self, index, name, units=None, **kwargs ):
        self.index = index
        self.name = name
        self.units = units
        self.attributes = kwargs

    def __getitem__(self, key): return self.attributes[key]

    def __setitem__(self, key, value): self.attributes[key] = value

    def toXMLList( self, indent='' ):
        xmlStr = '%s<column index="%i" name="%s"' % (indent,self.index,self.name)
        if self.units: xmlStr += ' units="%s"' % self.units
        for key in sorted(self.attributes):
            if self.attributes[key] is not None: xmlStr += ' %s="%s"' % (key,self.attributes[key])
        xmlStr += '/>'
        return [xmlStr]

class blank:
    """ Blank table entry, to indicate missing data """
    def __init__( self ): pass
    def __str__( self ): return '_'
    def __add__( self, other ): return self
    def __radd__( self, other ): return self
    def __mul__( self, other ): return self
    def __rmul__( self, other ): return self

def parseXMLNode(element, conversionTable={}):
    """Read a table element from xml into python. To convert a column or attribute from string to some other type,
    enter the new type in the conversionTable: {'index':int, 'scatteringRadius':PhysicalQuantityWithUncertainty, etc}. """
    def fixAttributes( items ):
        attrs = dict(items)
        for key in attrs:
            if key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
        return attrs

    def floatOrBlank( val ):
        if val=='_': return blank()
        return float( val )

    nRows,nColumns = int( element.get('rows') ), int( element.get('columns') )
    columns, data = element[:]
    columns = [columnHeader( **fixAttributes(column.items()) ) for column in columns]

    if data.text: data = map(floatOrBlank, data.text.split())
    else: data = []
    for i in range(len(columns)):
        if columns[i].name in conversionTable:
            data[i::nColumns] = map( conversionTable[columns[i].name], data[i::nColumns] )
    assert len(data) == nRows * nColumns
    data = [data[i*nColumns:(i+1)*nColumns] for i in range(nRows)]
    return table( columns, data )


if __name__ == '__main__':
    """Sample uses for the table class: """
    tt = table( columns=[
        columnHeader( 0, 'energy', 'eV' ),
        columnHeader( 1, 'spin', '' ),
        columnHeader( 2, 'neutronWidth', 'eV' ),
        columnHeader( 3, 'captureWidth', 'eV' ), ] )
    for dat in ([-1.7, blank(), 3.2, 0.089],
            [3.4, 0.5, 4.2, 0.072],
            [5.6, 1.5, 2.76, blank()]):
        tt.addRow( dat )
    xmlstring = '\n'.join( tt.toXMLList() )
    print xmlstring

    from xml.etree import cElementTree as parser
    element = parser.fromstring( xmlstring )
    tt2 = parseXMLNode( element, {'index':int} )
    xmlstring2 = '\n'.join( tt2.toXMLList() )
    assert xmlstring == xmlstring2

    # extract a column, converting eV -> MeV
    print tt2.getColumn('captureWidth','MeV')
