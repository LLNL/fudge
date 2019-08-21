# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry

matrixToken = 'matrix'
diagonalFormToken = 'diagonal'
symmetricFormToken = 'symmetric'
sparse_symmetricFormToken = 'sparse_symmetric'
asymmetricFormToken = 'asymmetric'
sparseFormToken = 'sparse_asymmetric'

class matrix( ancestry ):
    """
    Class for storing a matrix or two-dimensional array.
    This class supports diagonal, symmetric, asymmetric, and sparse matrices.
    Matrix data is stored as a list of lists, which is easily converted into a numpy array.
    """

    def __init__(self, data, form=None, precision=None):
        """
        @data: list-of-lists containing matrix data
        @form: diagonal, symmetric, sparse, etc
        @precision: number of digits when printing. Necessary only if 'pretty' matrix output is desired.
        """
        ancestry.__init__( self, matrixToken, None )
        self.form = form
        self.nrows, self.ncols = len(data),len(data[0])
        self.data = data
        self.precision = precision
    
    def __getitem__(self, idx):
        return self.data[idx]

    def check( self, info ): 
        """ 
        Run some simple checks on covariance matrices: 
            * look for negative eigenvals,
            * check ratio of smallest/largest eigenvals 
        """
        from fudge.gnd import warning
        warnings = []

        if self.form in (symmetricFormToken, diagonalFormToken, sparse_symmetricFormToken):
            import numpy
            A = numpy.matrix( self.data )
            # Check that the matrix is positive definite
            vals, vecs = numpy.linalg.eigh( A )
            if min(vals) < info['negativeEigenTolerance']:
                warnings.append( warning.negativeEigenvalues( len(vals[vals<0]), min(vals), self ))
            minpos, maxpos = min(vals[vals>=0]),max(vals[vals>=0])
            # Check that the condition number of the matrix is reasonable
            if minpos/maxpos < info['eigenvalueRatioTolerance']: 
                warnings.append( warning.badEigenvalueRatio( minpos/maxpos, self ) )
        return warnings

    def fix( self, **kw ):
        """
        Fix common eigenvalue problems: 
            * negative eigenvals,
            * bad ratio of smallest/largest eigenvals 
        """
        warnings = []
        if self.form in (symmetricFormToken,diagonalFormToken,sparse_symmetricFormToken) and ( kw['removeSmallEVs'] or kw['removeNegativeEVs'] ) :
            import numpy
            A = numpy.matrix( self.data )
            # Decompose A for later fixing
            e, O = numpy.linalg.eigh( A ) # stronger and more general than eigh
            v = [ O.T[i] for i in range( A.shape[0] ) ] # the eigenvectors are the columns of the orthogonal matrix O
            B = numpy.zeros_like( A )
            minval, maxval = min(e[e>=0]), max(e[e>=0]) # min and max eigenvalues
            ratio = minval/maxval
            evTol = max( kw['eigenvalueRatioTolerance'] * maxval, kw['eigenvalueAbsoluteTolerance'] )
            # Now construct B using the good eigenspaces of A
            for i in range( A.shape[0] ):
                if e[i] < 0.0 and ( abs( e[i] ) > max( evTol, abs( kw['negativeEigenTolerance'] ) ) ): #BAD BAD BAD    
                    warnings.append("WARNING: very bad bad bad eigenvalue = %e  for MF/MT=%s\n" % (e[i], kw['rowENDF_MFMT']))
                    if kw['removeSmallEVs'] and abs( e[i] ) < evTol: continue
                    if kw['removeNegativeEVs'] and e[i] < 0.0 : continue
                try:
                    B += numpy.outer( e[i]*v[i], v[i].T )
                except TypeError as err:
                    print B, e[i], v[i]
                    raise err
            self.data = B.tolist()
        return warnings
        
    def plot( self ):
        from numpy import mat
        from matplotlib.pyplot import matshow, show
        matshow( mat( self.data ) )
        show()


    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<matrix rows="%s" columns="%s" form="%s"' 
                % (self.nrows, self.ncols, self.form)]
        if self.precision: xmlString[0]+=' precision="%i"' % self.precision
        xmlString[0]+= '>'
        
        if self.precision: format = '%% -.%ie' % self.precision
        elif flags: format = flags.get('covarianceFormat') or '%s'
        else: format = '%s'

        # for some matrix forms we can save lots of space on xml output:
        if self.form == diagonalFormToken:
            for row in range(self.nrows):
                xmlString[-1] += format % self.data[row][row] + ' '
        elif self.form == symmetricFormToken:
            for row in range(self.nrows):
                nCols = row+1
                template = indent+'  '+' '.join([format]*nCols)
                xmlString.append( template % tuple(self.data[row][:nCols]) )
        elif self.form == asymmetricFormToken:
            for row in self.data:
                template = indent+'  '+' '.join([format]*len(row))
                xmlString.append( template % tuple(row) )
        elif self.form in (sparse_symmetricFormToken, sparseFormToken):
            # sparse matrices: give non-zero elements by listing
            #  row column length ...data...
            # For example, a diagonal 2x2 matrix would be:
            #
            # <matrix rows="2" columns="2" form="sparse_symmetric">
            #   0 0 1 0.234
            #   1 1 1 0.173
            # </matrix>

            # the buffer (buf) determines maximum distance between non-zero elements
            # before a new index must be given
            buf = 6
            rows,cols = self.nrows, self.ncols
            for row in range(rows):
                rowString = ''
                if self.form == sparse_symmetricFormToken: end = row+1
                else: end = self.ncols
                dat = self.data[row][:end]
                col, stop = 0, 0
                while col<end:
                    if dat[col]!=0:
                        stop = col+1
                        while stop < end:
                            if dat[stop]!=0: stop += 1
                            elif any(dat[stop:stop+buf]):
                                for i in range(buf):
                                    if stop+i<end and dat[stop+i]!=0: stop += i
                            else: break
                        template = ' '.join([format]*(stop-col))
                        rowString += '  %i %i %i ' % (row,col,stop-col) + template % tuple(dat[col:stop])
                        col = stop
                    col += 1
                if rowString:
                    xmlString.append(indent+rowString)
        else:
            raise TypeError("Don't know how to write matrix of form %s" % self.form)
        xmlString[-1] += '</matrix>'
        return xmlString

def parseXMLNode( mat, linkData={} ):
    import numpy
    if mat.text: matrixDat = map(float, mat.text.split())
    else: matrixDat = []
    nRows, nCols = int( mat.get('rows') ), int( mat.get('columns') )
    form = mat.get('form')

    if form==diagonalFormToken:
        if len(matrixDat) != nRows:
            raise TypeError("Wrong number of elements for diagonal matrix!")
        tmp = numpy.eye( nRows ) * numpy.array( matrixDat )
    elif form==symmetricFormToken:
        if len(matrixDat) != nRows * (nCols+1)/2:
            raise TypeError("Wrong number of elements for symmetric matrix!")
        start = 0
        tmp = numpy.zeros((nRows,nCols))
        for i in range(nRows):
            tmp[i,:i+1] = matrixDat[start:start+i+1]
            tmp[:i,i] = matrixDat[start:start+i]    # upper-diagonal portion
            start += i+1
    elif form==asymmetricFormToken:
        if len(matrixDat) != nRows * nCols:
            raise TypeError("Wrong number of elements in matrix!")
        tmp = numpy.array([matrixDat[j*nCols:(j+1)*nCols] for j in range(nRows)])
    elif form in (sparseFormToken, sparse_symmetricFormToken):
        tmp = numpy.zeros((nRows,nCols))
        idx = 0
        while idx<len(matrixDat):
            row, col, length = map(int, matrixDat[idx:idx+3])
            tmp[row, col:col+length] = matrixDat[idx+3: idx+3+length]
            idx += 3+length
        if form == sparse_symmetricFormToken:
            # symmetrize
            tmp = tmp + tmp.T - numpy.diag( tmp.diagonal() )
    else:
        raise Exception("Don't know how to read matrix of form %s!" % form)

    if 'precision' in mat.keys():
        precision=int(mat.get("precision"))
    else:
        precision=None

    return matrix( tmp, form, precision=precision )
