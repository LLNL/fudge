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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""Base classes for covariances: matrix, axes."""
from fudge.core.ancestry import ancestry
from . import tokens
from fudge.core.math import matrix as gndMatrix
from fudge.core.math.xData import XYs, axes
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

class covarianceMatrix( ancestry ):
    """
    Base class for all covariances. covarianceMatrix contains axes, a list of energy boundaries,
    and the matrix data. Some details :

        * Matrix data is stored in an :py:class:`xData.matrix class`. May be diagonal, symmetric, sparse, etc
        * Symmetric matrices only require one set of energy bounds, but asymmetric
          matrices require bounds for both axes.
          
    """

    moniker = tokens.covarianceFormToken

    def __init__(self, index=None, type=tokens.absoluteToken, axes=None, matrix=None, energyBounds=None,
            ENDFconversionFlag=None):
        ancestry.__init__( self, 'covarianceMatrix', None, attribute = 'index' )
        self.index = index  #: an int, required inside a 'mixed' section
        self.type = type #: 'relative' or 'absolute'
        self.axes = axes or [] #: list of covarianceAxis instances.  There should be 3: the rows, the columns and the value of the covariance itself
        self.matrix = matrix or [] #: a :py:class:`fudge.core.math.matrix` instance containing the actual covariance matrix
        self.energyBounds = energyBounds #: duh, the energy bounds
        self.ENDFconversionFlag = ENDFconversionFlag #: yes, this is a crutch to help when converting back to ENDF

    def convertAxesToUnits( self, units ):
        '''
        Converts all the axes' units.
        The parameter ``units`` should be a list of units with the same length as self.axes
        '''
        if not type( units ) in [ list, tuple ]: raise TypeError()
        if len( units ) != len( self.axes ): raise ValueError()
        for i,a in enumerate( self.axes ): a.convertToUnit( units[i] )

    def toCovarianceMatrix( self ): 
        import copy
        return copy.copy( self )
        
    def toCorrelationMatrix( self ):
        '''
        Returns the correlation matrix generated from self's covariance matrix.  This is
        essentially a copy of self, but renormalized by the uncertainty:
        
            correlation[i,j] = covariance[i,j]/sqrt(covariance[i,i])/sqrt(covariance[j,j])
        
        We reuse the covariance matrix class so that we can do plotting, etc.  If you have
        a correlation matrix, you can safely recover it provided you have the uncertainty 
        vector.
        
        Currently only works for a square covariance matrix and not a off-diagonal part of 
        another covariance.
        '''
        # Check if is full, square covariance matrix
        if not self.matrix.form in (gndMatrix.symmetricFormToken, 
            gndMatrix.diagonalFormToken,gndMatrix.sparse_symmetricFormToken): 
                raise TypeError( "Can only extract correlation matrices from symmetric covariance matrices" )
                
        # Copy to result matrix
        import copy
        correlation = copy.copy( self )
        
        # Rescale result matrix
        import numpy
        theCorrelationMatrix = numpy.matrix( self.matrix.data )
        theUncertainty = numpy.diag( theCorrelationMatrix )
        theUncertainty[ theUncertainty < 0.0 ] = 0.0
        theUncertainty = numpy.sqrt( theUncertainty )
        for i in range( theCorrelationMatrix.shape[0] ):
            for j in range( theCorrelationMatrix.shape[1] ):
                theCorrelationMatrix[i,j] /= ( theUncertainty[i] * theUncertainty[j] )

        # Return the result
        correlation.axes[2].unit = ''
        correlation.matrix.data = theCorrelationMatrix.tolist()
        return correlation
        
    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param XYs rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param XYs colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note:    If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        import copy
        result = copy.copy( self )
        if self.type==tokens.absoluteToken: return result

        # Make sure we have usable data to rescale with
        if not isinstance( rowData, XYs.XYs ): raise TypeError( 'rowData must be of type XYs, found '+str(type(rowData)) )
        gRowData = rowData.group( self.axes[0].data, self.axes[0].unit )
        if isinstance( colData, XYs.XYs ): gColData = colData.group( self.axes[1].data, self.axes[1].unit )
        else: gColData = gRowData
        
        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.data):
            newRow = []
            for j, cell in enumerate(row):
                newRow.append( gRowData[i][1]*gColData[j][1]*cell )
            newData.append( newRow )
        result.matrix.data = newData
        result.type=tokens.absoluteToken
        return result
        
    def toRelative( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a absolute covariance) using XYs rowData and colData
        to convert self into a relative covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        import copy
        result = copy.copy( self )
        if self.type==tokens.relativeToken: return result

        # Make sure we have usable data to rescale with
        if not isinstance( rowData, XYs.XYs ): raise TypeError( 'rowData must be of type XYs, found '+str(type(rowData)) )
        gRowData = rowData.group( self.axes[0].data, self.axes[0].unit )
        if not self.axes[1].mirrorOtherAxis:
            if not isinstance( colData, XYs.XYs ): raise TypeError( 'colData must be of type XYs, found '+str(type(colData)) )
            gColData = colData.group( self.axes[1].data, self.axes[1].unit )
        else: gColData = gRowData

        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.data):
            newRow = []
            for j, cell in enumerate(row):
                newRow.append( cell/gRowData[i][1]/gColData[j][1] )
            newData.append( newRow )
        result.matrix.data = newData
        result.type=tokens.relativeToken
        return result
        
    def check( self, info ): 
        """Check if uncertainty in the bounds passed into the checker.  
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        from fudge.gnd import warning
        warnings = []

        if self.matrix.form in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken) and info['checkUncLimits']: 
            import numpy
            A = numpy.matrix( self.matrix.data )
            relative = self.type == 'relative'
            if relative:
                varMin = info['minRelUnc']*info['minRelUnc']
                varMax = info['maxRelUnc']*info['maxRelUnc']
            for i in range( A.shape[0] ):
                if not relative:
                    if info['theData'] != None:
                        uncMin = info['minRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[0].data[i]+self.axes[0].data[i+1]) )
                        uncMax = info['maxRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[0].data[i]+self.axes[0].data[i+1]) )
                        varMin = uncMin*uncMin
                        varMax = uncMax*uncMax
                    else:
                        #warnings.append( "WARNING: can't check absolute uncertainties without data to compare to!\n" )
                        break
                if varMin <= A[i,i] and varMax >= A[i,i]: pass # unc is where is should be
                elif varMin >= A[i,i]:                         # uncertainty too small
                    warnings.append( warning.varianceTooSmall( i, A[i,i], self ) )
                else:                                          # uncertainty too big
                    warnings.append( warning.varianceTooLarge( i, A[i,i], self ) )
        return warnings + self.matrix.check( info )
    
    def fix( self, **kw ): 
        """Fix uncertainty using the bounds passed into the fixer.
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        warnings = []
        if self.matrix.form in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken) and kw['fixUncLimits']: 
            import numpy
            A = numpy.matrix( self.matrix.data )
            relative = self.type == 'relative'
            for i in range( A.shape[0] ):
                eMin = self.axes[0].data[i]
                eMax = self.axes[0].data[i+1]
                if relative:
                    uncMin = kw['minRelUnc']
                    uncMax = kw['maxRelUnc']
                else:
                    uncMin = kw['theData'].getValue( 0.5*(eMax+eMin) )
                    uncMax = kw['theData'].getValue( 0.5*(eMax+eMin) )
                uncMin2 = uncMin*uncMin
                uncMax2 = uncMax*uncMax
                #eThresh = threshold.getValueAs( component.axes[0].units )
                if uncMin2 <= A[i,i] and uncMax2 >= A[i,i]: pass   # unc is where is should be
                elif uncMin2 >= A[i,i]: 
                    if i+1 < A.shape[0] and A[i+1,i+1] >= uncMin2: A[i,i] = A[i+1,i+1]
                    else: A[i,i] = uncMin2
                # else:                                            # above threshold and uncertainty out of bounds
                    # continue #skip this fix for now
                    # if uncMin2 > A[i,i]: # uncertainty too small
                        # print '    WARNING: bin', i, 'uncertainty is too small!!!', '(', uncMin2, '>', A[i,i], ')'
                        # if A[i, i] == 0.0:
                            # for j in range( i+1, A.shape[0] ): 
                                # A[i,j] = 0.0
                                # A[j,i] = 0.0
                            # A[i,i] == uncMin2
                        # else:
                            # for j in range( i+1, A.shape[0] ): 
                                # A[i,j] *= uncMin / math.sqrt( A[i,i] )
                                # A[j,i] = A[i,j]
                            # A[i,i] == uncMin2
                    # else:                # uncertainty too big
                        # print '    WARNING: bin', i, 'uncertainty is too big!!!', '(', A[i,i], '>', uncMax2, ')'
                        # for j in range( i+1, A.shape[0] ): 
                            # A[i,j] *= uncMax / math.sqrt( A[i,i] )
                            # A[j,i] = A[i,j]
                        # A[i,i] = uncMax2
            self.data = B.tolist()
        return warnings + self.matrix.fix( **kw )
        
    def group( self, groupBoundaries = ( None, None ), groupUnit = ( None, None ) ):
        '''
        Group the matrix in self
        
        :param original: the original covarianceMatrix instance we intend to regroup
        :param groupBoundaries: a 2 element list containing the group boundaries for the rows 
                                and columns (respectively) of the covariance to be regrouped
        :param groupUnit: a 2 element list containing the units in which group boundaries are 
                          specified for the rows and columns (respectively) of the covariance 
                          to be regrouped
            
        .. note::  We still need to do flux weighting
            
            
        .. rubric:: Regrouping Theory
        
        Given a function :math:`f(E)`, we write the grouped data using fudge's ``flat`` interpolation 
        scheme.  We note that we could write this scheme as an expansion over basis functions:
        
        .. math::    
            f(E) = \sum_{i=0}^{N+1} w_i(E) * f_i
        
        where the weight functions :math:`w_i(E)` are
        
        .. math::
            w_i(E) = 1  \;\\text{for}\; E_i <= E <= E_{i+1}; \;\; 0 \;\\textrm{otherwise}
            
        These weights are an orthogonal (bot not orthonormal) basis, with 
        
        .. math::
            (E_{i+1}-E_i) \delta_{ij} = \int dE w_i(E) * w_j(E)
        
        So, to transform from basis :math:`w_i(E)` to :math:`v_i(E)` (which has group boundaries 
        :math:`[ E'_0, ... ]`), do: 
        
        .. math::
            f'_j = \sum_i m_{ji} f_i
            
        where :math:`f'` is the regrouped function coefficients and :math:`m_{ji}` is the matrix
        
        .. math::
            m_{ij} = (E'_{i+1}-E'_i)^{-1} \int dE v_i(E) w_j(E) 

            
        .. rubric:: Applying regrouping theory to covariance matrices   
        
        When we are given a covariance matrix :math:`c_{ij}` in ENDF, it is meant to be interpreted
        as a grouped covariance in both the direction of the matrix rows and the matrix 
        columns.  Therefore, we must regroup in both the row direction and the column 
        direction.  The ENDF format gives both the group boundaries for the rows and columns.
        In other words, ENDF gives us the following rule for evaluating the continuous row-
        column covariance:
        
        .. math::
            c( E1, E2 ) = \sum_{ij} w_i(E1) w_j(E2) c_{ij}
            
        Computing :math:`m_{ij}` as before, 
            
        .. math::
            cc_{ij} = \sum_{i',j'} m_{ii'} c_{i'j'} m_{j'j}
            
        It is straightforward to generalize to the case where the row and column bases are 
        different.
        
        In the routine below, we abuse :py:class:`fudge.core.math.xData.XYs` to specify the functions 
        :math:`w_i(E)` and use the :py:func:`XYs.groupOneFunction()` method to perform the integrals to get
        the regrouping matrix.  We do this separately for the rows and the columns.
        The matrix multiplication that converts a covariance from one pair of bases (group 
        structures) to another is accomplished using numpy.

        
        .. rubric:: An explanation of fudge's 'flat' interpolation
        
        Suppose we have a function :math:`f(E)` specified using fudge's `'flat'` interpolation.  
        Then we have :math:`N` entries :math:`[f_0, f_1, ..., f_{N-1}]` and a set of group 
        boundaries :math:`[E_0, E_1, ..., E_N]` and the following rule for interpolation:
        
            * Below :math:`E_0`, :math:`f(E)` evaluates to :math:`0.0`
            * From :math:`E_0 \\rightarrow E_1`, :math:`f(E)` evaluates to :math:`f_0`
            * From :math:`E_1 \\rightarrow E_2`, :math:`f(E)` evaluates to :math:`f_1`
            * ...
            * From :math:`E_{i} \\rightarrow E_{i+1}`, :math:`f(E)` evaluates to :math:`f_i`
            * ...
            * From :math:`E_{N-1} \\rightarrow E_N`, :math:`f(E)` evaluates to :math:`f_{N-1}`
            * Above :math:`E_N`, :math:`f(E)` evaluates to :math:`0.0`
        '''
        import numpy, copy

        # determine where to get the settings for the potentially mirrored second axis
        if self.axes[1].mirrorOtherAxis:    axis1index = 0
        else:                               axis1index = 1
        
        # setup the old axes in a form we can (ab)use in the XYs class
        axes0_ = axes.defaultAxes( 
            labelsUnits={ 
                0:( self.axes[0].label, self.axes[0].unit ),
                1:( 'dummy', '' )},
            dependentInterpolation='flat' )        
        axes1_ = axes.defaultAxes( 
            labelsUnits={ 
                0:( self.axes[axis1index].label, self.axes[axis1index].unit ),
                1:( 'dummy', '' )},
            dependentInterpolation='flat' )        
        
        # define basis functions for the rows and columns
        basis0 = XYs.XYs( axes0_, [ ( x, 0.0 ) for x in self.axes[0].data ], 0.0001 )
        basis1 = XYs.XYs( axes1_, [ ( x, 0.0 ) for x in self.axes[axis1index].data ], 0.0001 )
        basis0 = basis0.convertAxisToUnit( 0, groupUnit[0] )
        basis1 = basis1.convertAxisToUnit( 0, groupUnit[1] )
    
        # build the regrouping matrices for the two bases
        w0 = []
        for i in range( self.matrix.nrows ):
            basis0[i] = ( basis0[i][0], 1.0 )
            w0.append( basis0.groupOneFunction( groupBoundaries[0], norm = 'dx' ) )
            basis0[i] = ( basis0[i][0], 0.0 )
        w0 = numpy.mat( w0 )
        w1 = []
        for j in range( self.matrix.ncols ):
            basis1[j] = ( basis1[j][0], 1.0 )
            w1.append( basis1.groupOneFunction( groupBoundaries[1], norm = 'dx' ) )
            basis1[j] = ( basis1[j][0], 0.0 )
        w1 = numpy.mat( w1 )
                
        # set up the regrouped covariance matrix
        grouped = copy.copy( self )
        grouped.axes[0].data = groupBoundaries[0]
        grouped.axes[1].data = groupBoundaries[1]
        grouped.axes[0].unit = groupUnit[0]
        grouped.axes[1].unit = groupUnit[1]
        odata = numpy.mat( self.matrix.data )
        gdata = w0.T * odata * w1
        grouped.matrix.data = gdata.tolist()
        grouped.matrix.nrows = len(grouped.matrix.data)
        grouped.matrix.ncols = len(grouped.matrix.data[0]) 
        return grouped

    def removeExtraZeros(self):
        import numpy
        matrix = numpy.array( self.matrix.data )
        rowStart, colStart = 0,0
        rowEnd, colEnd = matrix.shape
        while numpy.all(matrix[rowStart,:]==0):
            rowStart += 1
        while numpy.all(matrix[:,colStart]==0):
            colStart += 1
        while numpy.all(matrix[rowEnd-1,:]==0):
            rowEnd -= 1
        while numpy.all(matrix[:,colEnd-1]==0):
            colEnd -= 1

        matrix = matrix[rowStart:rowEnd, colStart:colEnd]
        if len(self.axes)==1:
            assert (rowStart,rowEnd) == (colStart,colEnd)
            self.axes[0].data = self.axes[0].data[rowStart:rowEnd]
        else:
            self.axes[0].data = self.axes[0].data[rowStart:rowEnd]
            self.axes[1].data = self.axes[1].data[colStart:colEnd]
        self.matrix.data = matrix.tolist()
        self.matrix.nrows, self.matrix.ncols = matrix.shape

    def getUncertaintyVector( self, theData=None, relative=True ):
        """ 
        Get an XYs object containing uncertainty for this matrix.
        Convert relative/absolute if requested (if so, must also pass central values as theData)

        Examples:

            - if the covariance matrix is relative and we want relative uncertainty vector, just do:
            
                >>> matrix.getUncertaintyVector()
                
            - if we want the absolute matrix instead:
            
                >>> matrix.getUncertaintyVector( theData=<XYs instance>, relative=False ) 
                
        """

        if self.matrix.form not in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken):
            raise ValueError("getUncertaintyVector only applies to symmetric matrices!")
        import numpy
        energies = self.axes[0].data
        diag = list( numpy.sqrt( numpy.diagonal( numpy.array( self.matrix.data ) ) ) )
        diag.append( diag[-1] )                             # point corresponding to final energy bin
        yunit = self.axes[-1].unit
        if yunit != '': # get square root of the unit
            yunit = PhysicalQuantityWithUncertainty(1,yunit).sqrt().getUnitSymbol()
        axes_ = axes.defaultAxes( labelsUnits={0:('energy_in',self.axes[0].unit),1:('uncertainty',yunit)},
                dependentInterpolation='flat' )
        uncert = XYs.XYs( axes_, zip(energies,diag), 0.0001 )    # what should accuracy be set to?
        uncert = uncert.changeInterpolation('linear','linear',None,1e-8,1e-8)

        # do we need to convert absolute->relative or vice versa?
        if (relative and self.type==tokens.absoluteToken) or (not relative and self.type==tokens.relativeToken):
            if theData==None:
                raise ValueError("Need central values to convert relative<->absolute uncertainties!")
            try:
                theData = theData.toPointwise_withLinearXYs(1e-8,1e-8)
                uncert, theData = uncert.mutualify(1e-8, 1e-8, False, theData, 1e-8, 1e-8, False)
                if relative: #convert absolute to relative
                    uncert /= theData
                else: #relative to absolute
                    uncert *= theData
            except Exception as err:
                print len( uncert ), uncert.copyDataToXYs()[0], uncert.copyDataToXYs()[-1]
                print len( theData ), theData.copyDataToXYs()[0], theData.copyDataToXYs()[-1]
                raise Exception( err.message )
        return uncert
        
    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import QuadMesh
        
        x = self.axes[0].data
        if self.axes[1].mirrorOtherAxis:    y = self.axes[0].data 
        else:                               y = self.axes[1].data 

        X, Y = np.meshgrid( x, y )        
        XY = np.hstack((X.ravel()[:,np.newaxis], Y.ravel()[:,np.newaxis]))
        Z = (np.array(self.matrix.data)).ravel()
        
        ax = plt.subplot(1,1,1)
        if title == None: title = str( self.toXLink() )
        plt.suptitle(title)

        qc = QuadMesh(
            meshWidth=len(x)-1, 
            meshHeight=len(y)-1, 
            coordinates=XY, 
#            showedges=True, 
            antialiased=True, 
            shading='flat',
            transOffset=ax.transData)
            
        qc.set_array(Z) 
        ax.add_collection(qc,autolim=True)

        if xlim == None:    ax.set_xlim( x[0], x[-1] )
        else:               ax.set_xlim( xlim[0], xlim[1] )
        if ylim == None:    ax.set_ylim( y[0], y[-1] )
        else:               ax.set_ylim( ylim[0], ylim[1] )
        if xlog: ax.set_xscale( 'log' )
        if ylog: ax.set_yscale( 'log' )

        xlabel = self.axes[0].label + ' (' + self.axes[0].unit +')'
        if self.axes[1].mirrorOtherAxis:    ylabel = xlabel 
        else:                               ylabel = self.axes[1].label + ' (' + self.axes[1].unit +')'

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        cbar = plt.colorbar(qc)
        if scalelabel != None: cbar.set_label(scalelabel)
        else: cbar.set_label(str(self.type)+' covariance ('+str(self.axes[2].unit)+')')
        plt.show()

    def toXMLList(self, flags=None, indent=''):
        indent2 = indent+'  '; indent3 = indent2+'  '
        xmlString = [indent+'<%s' % self.moniker]
        if self.index!=None: xmlString[0] += ' index="%s"' % self.index
        xmlString[0] += ' type="%s"' % self.type
        if self.energyBounds: xmlString[0] += (
                ' lowerBound="%s" upperBound="%s"' %
                self.energyBounds )
        if self.ENDFconversionFlag: xmlString[0] += (
                ' ENDFconversionFlag="%s"' % self.ENDFconversionFlag )
        xmlString[0] += '>'
        xmlString.append( '%s<axes>' % indent2 )
        for axis in self.axes: xmlString.append( axis.toXML(indent=indent3) )
        xmlString[-1] += '</axes>'
        xmlString += self.matrix.toXMLList(flags, indent2)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode(element, xPath=[], linkData={}):
        """Translate <covarianceMatrix> element from xml into python class."""

        xPath.append( element.tag )
        ax, mat = element[0], element[1]
        axes_ = covarianceAxis.parseXMLNode( ax )
        matrix_ = gndMatrix.parseXMLNode( mat )
        CM = covarianceMatrix( type=element.get('type'), axes=axes_, matrix=matrix_,
                ENDFconversionFlag=element.get("ENDFconversionFlag") )
        xPath.pop()
        return CM
    
    def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
        from fudge.legacy.converting import endfFormats
        endf = []
        rowdat, coldat = targetInfo['dataPointer']
        MF,MT1 = map(int, rowdat['ENDF_MFMT'].split(','))
        if not inCovarianceGroup:
            # print header for this subsection (contains one NL sub-subsection)
            MAT1 = targetInfo['MAT1']
            XMF1,XLFS1,NC,NI = 0,0,0,1
            if coldat:
                MF1, MT1 = map(int, coldat['ENDF_MFMT'].split(','))
            if MF in (31,33):
                endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        # header for matrix:
        rows,cols = self.matrix.nrows, self.matrix.ncols
        if self.matrix.form==gndMatrix.diagonalFormToken:
            LS = 0; LB = 1; NP = len(self.axes[0]); NT = 2*NP
            if self.type=='absolute': LB = 0
            if self.ENDFconversionFlag:
                LB = int( self.ENDFconversionFlag.split('=')[1] )
            matrixData = []
            for p in range(NP-1):
                matrixData.extend( [self.axes[0][p], self.matrix[p][p] ] )
            matrixData.extend( [self.axes[0][NP-1], 0] )
        elif self.matrix.form==gndMatrix.symmetricFormToken:
            LS = 1; LB = 5; NT = (rows+1) + rows*(rows+1)/2; NP = rows+1
            matrixData = self.axes[0].data + [
                    i for j in range(rows) for i in self.matrix.data[j][j:]]
        elif self.axes[1].mirrorOtherAxis:
            LS = 0; LB = 5; NT = (rows+1) + rows*cols; NP = rows+1
            matrixData = self.axes[0].data + [
                    i for sublist in self.matrix.data for i in sublist]
        else:
            LS = 0; LB = 6; NT = (rows+1) + (cols+1) + rows*cols; NP = rows+1
            matrixData = self.axes[0].data + self.axes[1].data + [
                    i for sublist in self.matrix.data for i in sublist]
        if MF==35:  # header for fission spectra is different:
            E1,E2 = [a.getValueAs('eV') for a in self.energyBounds]
            if LS: LB = 7
            else:
                raise Exception ("Unknown spectrum (MF35) covariance format")
            endf.append( endfFormats.endfHeadLine( E1,E2,LS,LB,NT,NP ) )
        else:
            endf.append( endfFormats.endfHeadLine( 0,0,LS,LB,NT,NP ) )
        endf += endfFormats.endfDataList( matrixData )
        return endf

class covarianceAxis:
    """
    Store energy group boundaries (similar to xData.axes.axis).
    Covariances in ENDF are often symmetric, NxN matrices. The
    axis then must have N+1 points defining upper and lower bounds for each
    bin in the matrix. 
    """
    def __init__(self, index=0, label=None, unit=None, interpolation=None,
            data=None, mirrorOtherAxis=False):
        self.index = index #: an int
        self.label = label #: a str to label the axis
        self.unit = unit #: any acceptable unit for a PhysicalQuantityWithUncertainty
        self.data = data or [] #: a list of floats in unit "unit" denoting the nodes (in ENDF, they are all group boundaries) for interpolating the rows and or columns of a covariance matrix
        self.interpolation = interpolation #: a string describing how the the covariance is meant to be interpreted.  Currently all ENDF covariance is 'linear,flat'
        self.mirrorOtherAxis = mirrorOtherAxis #: if True, then self is a copy of the accompanying row/column axis

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]
        
    def convertToUnit( self, unit ):
        '''
        Converts the units of a covarianceAxis, but only if the covarianceAxis isn't a mirror of another one
        '''
        if self.mirrorOtherAxis: return
        self.data = [ PhysicalQuantityWithUncertainty( x, self.unit ).inUnitsOf( unit ).getValue() for x in self.data ]
        self.unit = unit

    def toXML(self, indent = ''):
        from pqu.physicalQuantityWithUncertainty import toShortestString
        interpolationFlag, axisFlag, lengthFlag = '','',''
        if self.interpolation: interpolationFlag=' interpolation="%s"' % self.interpolation
        if self.mirrorOtherAxis: axisFlag=' mirror_row_energy_bounds="true"'
        if len(self.data)!=0: lengthFlag = ' length="%i"' % len(self.data)
        xmlString = indent+'<axis index="%i" label="%s" unit="%s"%s%s%s>' % (
                self.index, self.label, self.unit, interpolationFlag, axisFlag, lengthFlag)
        # self.data is usually list of floats (energy boundaries). Could also be list of input parameters
        if len(self.data)!=0:
            xmlString += ''.join([' %s' % toShortestString(val) for val in self.data])
            xmlString += '</axis>'
        else: xmlString = xmlString[:-1] + '/>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate covariance <axes> element from xml."""

        xPath.append( element.tag )
        axes_ = []
        for axis in element:
            data = None
            if axis.text is not None: data = map(float, axis.text.split())
            mirrorOtherAxis = False
            if axis.get('mirror_row_energy_bounds') in ('true','True'): mirrorOtherAxis=True
            interp = axis.get('interpolation')
            if interp is not None: interp = axes.interpolationXY( *interp.split(',') )
            axes_.append( covarianceAxis( index=int(axis.get('index')), label=axis.get('label'),
                unit=axis.get('unit'), data=data, interpolation=interp, mirrorOtherAxis=mirrorOtherAxis ) )
        xPath.pop()
        return axes_

