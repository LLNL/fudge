import pqu.PQU as PQUModule
import json
import csv

from xData import table as tableModule

class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, complex):
            return [obj.real, obj.imag]
        if isinstance(obj, set):
            return tuple(obj)
        if isinstance(obj, PQUModule.PQU):
            return str(obj)
        if isinstance(obj, PQUModule.PhysicalUnit):
            return str(obj)
        if isinstance(obj, DataTable):
            return '\n'.join(obj.toStringList())
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class DataTable(tableModule.Table):

    def __init__(self, columns=None, rows=None, data=None):

        tableModule.Table.__init__(self, columns=columns, data=data)
        self.rows = rows
        if self.rows and len(self.rows) != self.nRows:
            raise ValueError("Data is the wrong shape for a table with %i rows!" % self.nRows)

    def __str__(self):
        return '\n'.join(self.toXML_strList())

    def hasMatchingRow(self, key, value):
        try:
            icol = [x.name for x in self.columns].index(key)
        except ValueError:
            return False
        for row in self.data:
            if row[icol] == value:
                return True
        return False

    def getMatchingRow(self, key, value):
        icol = [x.name for x in self.columns].index(key)
        for i, row in enumerate(self.data):
            if row[icol] == value:
                return i, self.rows[i], row
        return None, None, None

    def toStringList(self, indent='', **kwargs):

        # If not row labels, revert to base class behavior
        if not self.rows:
            return tableModule.Table.toStringList(self, indent=indent, kwargs=kwargs)

        addHeader = kwargs.get('addHeader', True)
        addHeaderUnits = kwargs.get('addHeaderUnits', True)
        outline = kwargs.get('outline', False)
        columnWidths = [0] * (self.nColumns + 1)
        for col in range(self.nColumns):
            columnDat = [row[col] for row in self.data if not isinstance(row[col], tableModule.Blank)]
            asStrings = list(map(PQUModule.toShortestString, columnDat))
            columnWidths[col + 1] = max(list(map(len, asStrings)))
        columnWidths[0] = max(list(map(len, self.rows)))

        if addHeader:
            """ put column labels at the top of the table """
            names = ['name'] + [col.name for col in self.columns]
            lengths = [len(name) for name in names]
            if any([' width' in name for name in names]):
                # special treatment for RML widths: split up onto two lines
                nameL = [[], []]
                for name in names:
                    if ' width' in name:
                        l, r = name.split(' width')
                        r = ' width' + r + ' '
                    else:
                        l, r = name, ''
                    nameL[0].append(l)
                    nameL[1].append(r)
                names = nameL
                lengths = [len(name) for name in names[0]]
            else:
                names = [names]
            columnWidths = [max(columnWidths[i], lengths[i]) for i in range(len(columnWidths))]
            header = ['%s<!--' % (' ' * (len(indent) - 1)) + ' | '.join(
                [('%%%is' % columnWidths[i]) % nameList[i] for i in range(len(columnWidths))]) + '  -->' for nameList in
                      names]
            xml = header

        if addHeaderUnits:
            """ put column units at the top of the table """
            units = [''] + [str(col.unit).replace('None', '') for col in self.columns]
            header = ['%s<!--' % (' ' * (len(indent) - 1)) + ' | '.join(
                [('%s' % u).rjust(columnWidths[i]) for i, u in enumerate(units)]) + '  -->']
            xml += header

        template = ['%s' % (indent + ' ')] + ['%%%is' % l for l in columnWidths]

        def toString(val):
            if isinstance(val, tableModule.Blank):
                return str(val)
            if isinstance(val, str):
                return val
            return PQUModule.toShortestString(val)

        if outline:
            xml += [('   '.join(template) % tuple(map(toString, [self.rows[irow]] + dataRow))).rstrip() for
                    irow, dataRow in enumerate(self.data[:3])]
            xml += ['%s ...' % indent]
            xml += [('   '.join(template) % tuple(map(toString, [self.rows[-3 + irow]] + dataRow))).rstrip() for
                    irow, dataRow in enumerate(self.data[-3:])]
        else:
            xml += [('   '.join(template) % tuple(map(toString, [self.rows[irow]] + dataRow))).rstrip() for
                    irow, dataRow in enumerate(self.data)]
        return xml

    def getElement(self, rowName, colName):
        column = [a for a in self.columns if a.name == colName]
        if not column:
            return None
        if len(column) > 1:
            raise ValueError("Column named '%s' is not unique!" % colName)
        icol = self.columns.index(column[0])
        row = [a for a in self.rows if a == rowName]
        if not row:
            return None
        if len(row) > 1:
            print("WARNING: Row named '%s' is not unique! %s" % (rowName, str(row)))
        #            raise ValueError("Row named '%s' is not unique!" % rowName)
        irow = self.rows.index(row[0])
        return self.data[irow][icol]

    def as_dict(self):
        return {
            'units':{self.columns[i].name:self.columns[i].unit for i in range(len(self.columns))},
            "datatable":{
                row:{self.columns[i].name:self.data[irow][i] for i in range(len(self.columns))} for irow, row in enumerate(self.rows)}}

    def json_report(self):
        return json.dumps(self.as_dict(), cls=ComplexEncoder)

    def csv_report(self):
        raise Exception("FIXME: Write me")
