# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

__metaclass__ = type

import ancestry as ancestryModule
import pqu.PQU as PQUModule

class table( ancestryModule.ancestry ):

    moniker = 'table'

    def __init__( self, columns = None, data = None ):
        ancestryModule.ancestry.__init__( self )
        self.columns = columns or []
        self.data = data or []
        if not all( [len(d)==self.nColumns for d in self.data] ):
            raise Exception, ("Data is the wrong shape for a table with %i columns!" % self.nColumns)

    @property
    def nColumns(self): return len(self.columns)

    @property
    def nRows(self): return len(self.data)

    def __len__( self ):
        return self.nRows

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

    def addColumn( self, columnHeader, index=None ):
        """ add another column, either at 'index' or at the end of the table """
        if index:
            self.columns.insert(index, columnHeader)
            [row.insert(index, 0) for row in self.data]
        else:
            self.columns.append( columnHeader )
            [row.append(0) for row in self.data]
        for idx, col in enumerate(self.columns): col.index = idx

    def getColumn( self, columnName, units=None ):
        """ get data from one column, identified by the column 'name' attribute.
        Convert results to unit if units are specified """
        column = [a for a in self.columns if a.name==columnName]
        if not column: return None
        if len(column) > 1: raise Exception("Column named '%s' is not unique!" % columnName)
        index = self.columns.index( column[0] )
        if units:
            cf = ( PQUModule.PQU(1, column[0].units) / PQUModule.PQU(1, units) ).getValue()
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
        for idx, col in enumerate(self.columns): col.index = idx

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )
        addHeader = kwargs.get( 'addHeader', True )
        outline = kwargs.get( 'outline', False )
        if len(self.data) < 10: outline = False

        xml = ['%s<%s rows="%i" columns="%i">' % (indent,self.moniker,self.nRows,self.nColumns)]
        xml.append( '%s<columnHeaders>' % (indent2) )
        for column in self.columns: xml += column.toXMLList( indent3 )
        xml[-1] += '</columnHeaders>'

        if not self.data:
            xml.append( '%s<data/></%s>' % (indent2, self.moniker) )
            return xml

        xml.append( '%s<data>' % (indent2) )
        columnWidths = [0] * self.nColumns
        for col in range(self.nColumns):
            columnDat = [row[col] for row in self.data if not isinstance(row[col], blank)]
            asStrings = map( PQUModule.toShortestString, columnDat )
            columnWidths[col] = max( map(len, asStrings) )

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

        def toString(val):
            if isinstance(val, blank): return str(val)
            return PQUModule.toShortestString(val)

        if outline:
            xml += [('   '.join(template) % tuple( map(toString, dataRow))).rstrip() for dataRow in self.data[:3]]
            xml += ['%s ...' % indent]
            xml += [('   '.join(template) % tuple( map(toString, dataRow))).rstrip() for dataRow in self.data[-3:]]
        else:
            xml += [('   '.join(template) % tuple( map(toString, dataRow))).rstrip() for dataRow in self.data]
        xml[-1] += '</data></%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):
        """Read a table element from xml into python. To convert a column or attribute from string to some other type,
        enter the new type in the conversionTable: {'index':int, 'scatteringRadius':PhysicalQuantityWithUncertainty, etc}. """

        xPath.append( element.tag )
        def fixAttributes( items ):
            attrs = dict(items)
            for key in attrs:
                if key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
            return attrs

        def floatOrBlank( val ):
            if val=='_': return blank()
            return float( val )

        conversionTable = linkData.get('conversionTable',{})
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
        Table = cls( columns, data )
        xPath.pop()
        return Table



class columnHeader:
    """ defines one column in a table """
    def __init__( self, index, name, units=None, **kwargs ):
        self.index = index
        self.name = name
        self.units = units
        self.attributes = kwargs

    def __getitem__(self, key): return self.attributes[key]

    def __setitem__(self, key, value): self.attributes[key] = value

    def toXMLList( self, indent = '', **kwargs ) :

        xmlStr = '%s<column index="%d" name="%s"' % ( indent, self.index, self.name )
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
    tt2 = parseXMLNode( element, xPath = [], linkData = {'conversionTable' : {'index':int} } )
    xmlstring2 = '\n'.join( tt2.toXMLList() )
    assert xmlstring == xmlstring2

    # extract a column, converting eV -> MeV
    print tt2.getColumn('captureWidth','MeV')
