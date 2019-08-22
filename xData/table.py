# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

from __future__ import print_function

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
            raise ValueError("Data is the wrong shape for a table with %i columns!" % self.nColumns)

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
            raise ValueError("New row has %i columns, should have %i!" % (len(dataRow),self.nColumns))
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

    def convertUnits( self, unitMap ):
        """
        unitMap is a dictionary of the form { 'eV' : 'MeV', 'b' : 'mb' }.
        Converts all columns whose units appear as keys in unitMap
        """

        for idx, column in enumerate(self.columns):
            if column.unit in unitMap:
                factor = PQUModule.PQU(1, column.unit).getValueAs(unitMap[column.unit])
                column.unit = unitMap[column.unit]
                for row in self.data:
                    row[idx] *= factor

    def getColumn( self, columnNameOrIndex, unit=None ):
        """ get data from one column, identified by the column 'name' attribute.
        Convert results to desired unit if specified """
        if isinstance(columnNameOrIndex, int):
            index = columnNameOrIndex
            column = self.columns[index]
        else:
            column = [a for a in self.columns if a.name==columnNameOrIndex]
            if not column: return None
            if len(column) > 1: raise ValueError("Column named '%s' is not unique!" % columnNameOrIndex)
            column = column[0]
            index = self.columns.index( column )
        if unit:
            cf = PQUModule.PQU(1, column.unit).convertToUnit( unit ).getValue()
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

    def toStringList(self, indent='',**kwargs):
        addHeader = kwargs.get( 'addHeader', True )
        addHeaderUnit = kwargs.get( 'addHeaderUnit', False )
        outline = kwargs.get( 'outline', False )
        columnWidths = [0] * self.nColumns
        xml = []
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

        if addHeaderUnit:
            """ put column unit at the top of the table """
            unit = [str(col.unit) for col in self.columns]
            lengths = [len(u) for u in unit]
            columnWidths = [max(columnWidths[i],lengths[i]) for i in range(self.nColumns)]
            header = ['%s<!--' % (' ' * (len(indent) - 1)) + ' | '.join( [('%s'%u).rjust(columnWidths[i]) for i,u in enumerate(unit) ] ) + '  -->']
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
        return xml

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )
        if len(self.data) < 10: outline = False

        xml = ['%s<%s rows="%i" columns="%i">' % (indent,self.moniker,self.nRows,self.nColumns)]
        xml.append( '%s<columnHeaders>' % (indent2) )
        for column in self.columns: xml += column.toXMLList( indent3 )
        xml[-1] += '</columnHeaders>'

        if not self.data:
            xml.append( '%s<data/></%s>' % (indent2, self.moniker) )
            return xml

        xml.append( '%s<data>' % (indent2) )
        xml+=self.toStringList(indent=indent3,kwargs=kwargs)
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
    def __init__( self, index, name, unit ):
        self.index = index
        self.name = name
        self.unit = unit

    def __str__(self): return '%s (%s)'%(self.name,self.unit)

    def toXMLList( self, indent = '', **kwargs ) :

        xmlStr = '%s<column index="%d" name="%s" unit="%s"/>' % ( indent, self.index, self.name, self.unit )
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
    print(xmlstring)

    from xml.etree import cElementTree as parser
    element = parser.fromstring( xmlstring )
    tt2 = table.parseXMLNode( element, xPath = [], linkData = {'conversionTable' : {'index':int} } )
    xmlstring2 = '\n'.join( tt2.toXMLList() )
    assert xmlstring == xmlstring2

    # extract a column, converting eV -> MeV
    print(tt2.getColumn('captureWidth', 'MeV'))
