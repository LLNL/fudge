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

"""Base classes for covariances: matrix, axes."""
import copy, numpy
from . import tokens

from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import gridded as griddedModule
from xData import array as arrayModule
from xData import XYs as XYsModule
from xData import link as linkModule
from xData import values as valuesModule

from pqu import PQU

__metaclass__ = type

class covarianceMatrix( ancestryModule.ancestry ):
    """
    Base class for all covariances. covarianceMatrix contains a label, a 'gridded' container with
    the matrix data + axes and a list of energy boundaries, and a covariance 'type' which might be 'absolute',
    'relative' or 'correlation'. Some details :

        * Matrix data is stored in an :py:class:`xData.gridded.gridded` class. May be diagonal, symmetric, sparse, etc
        * Symmetric matrices only require one set of energy bounds, but asymmetric
          matrices require bounds for both axes.
    """

    moniker = 'covarianceMatrix'

    def __init__(self, label, type=tokens.absoluteToken, matrix=None, energyBounds=None, ENDFconversionFlag=None):
        """

        :param label: refers either to a style (i.e., 'eval'), or to the index inside a 'mixed' covariance
        :param type: 'relative' or 'absolute'
        :type type: str
        :param energyBounds: stores the incident energy range of validity for a PFNS covariance matrix
        :type energyBounds: tuple( :py:class:`PQU.PQU`, :py:class:`PQU.PQU` )
        :param matrix: the covariance matrix itself
        :type matrix: :py:class:`xData.gridded.gridded`
        :param ENDFconversionFlag:  yes, this is a crutch to help when converting back to ENDF
        :type ENDFconversionFlag: bool
        :return:
        :rtype: covarianceMatrix
        """

        ancestryModule.ancestry.__init__( self )
        self.__label = label
        self.type = type #: 'relative' or 'absolute'
        self.energyBounds = energyBounds
        self.matrix = matrix #: a :py:class:`xData.gridded.gridded` instance containing the matrix
        self.matrix.ancestor = self
        self.ENDFconversionFlag = ENDFconversionFlag #: yes, this is a crutch to help when converting back to ENDF

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value

    def getValue( self, x, y ):
        ix = self.matrix.axes[2].getIndexOfValue(x)
        if self.matrix.axes[1].style=='link': iy=self.matrix.axes[2].getIndexOfValue(y)
        else: iy=self.matrix.axes[1].getIndexOfValue(y)
        return self.matrix.array.constructArray()[ix,iy]

    def getRowBounds(self,unit=None):
        '''Get the bounds of the row.  If unit is specified, return the bounds in that unit.'''
        return (
            PQU.PQU(self.matrix.axes[-1].values[0],  self.matrix.axes[-1].unit),
            PQU.PQU(self.matrix.axes[-1].values[-1], self.matrix.axes[-1].unit))

    def getColumnBounds(self,unit=None):
        '''Get the bounds of the column.  If unit is specified, return the bounds in that unit.'''
        if self.matrix.axes[-2].style=='link': return self.getRowBounds(unit)
        return (
            PQU.PQU(self.matrix.axes[-2].values[0],  self.matrix.axes[-2].unit),
            PQU.PQU(self.matrix.axes[-2].values[-1], self.matrix.axes[-2].unit))

    def isSymmetric(self):
        """
        Simple test to determin if an underlying matrix is symmetric
        :return:
        """

        if ( self.matrix.array.compression == arrayModule.diagonal.moniker or
                self.matrix.array.symmetry in (arrayModule.symmetryUpperToken,arrayModule.symmetryLowerToken) ):
            return True
        # could still be symmetric even if it doesn't use compression
        arr = self.matrix.array.constructArray()
        return numpy.all( arr==arr.T )

    def convertAxesToUnits( self, units ):
        '''
        Converts all the axes' units.
        The parameter ``units`` should be a list of units with the same length as self.axes
        '''
        if not type( units ) in [ list, tuple ]: raise TypeError()
        if len( units ) != len( self.matrix.axes ): raise ValueError()
        for i,a in enumerate( self.matrix.axes ): a.convertToUnit( units[i] )

    def toCovarianceMatrix( self ):
        """
        Return self, converted to a `base` instance.  Since it already is one, we just return a copy.
        :return:
        """
        return copy.copy( self ) # copy.deepcopy( self )
        
    def getCorrelationMatrix( self ):
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
        if not self.isSymmetric():
            raise TypeError( "Can only extract correlation matrices from symmetric covariance matrices" )
                
        # Rescale result matrix
        theCorrelationMatrix = self.matrix.array.constructArray()
        theUncertainty = copy.copy( numpy.diag( theCorrelationMatrix ) )
        theUncertainty[ theUncertainty < 0.0 ] = 0.0
        theUncertainty = numpy.sqrt( theUncertainty )
        for i1 in range( theCorrelationMatrix.shape[0] ):
            for i2 in range( theCorrelationMatrix.shape[1] ):
                theCorrelationMatrix[i1,i2] /= ( theUncertainty[i1] * theUncertainty[i2] )

        # Return the result
        tridata = theCorrelationMatrix[numpy.tri(theCorrelationMatrix.shape[0])==1.0].tolist()
        correlation = griddedModule.gridded(
                axes=self.matrix.axes.copy(),
                array=arrayModule.full(shape=theCorrelationMatrix.shape, data=tridata, symmetry=arrayModule.symmetryLowerToken) )
        correlation.axes[0].unit = ''
        return correlation
        
    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs1d rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param XYs1d rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite
        :param XYs1d colData: an XYs1d instance containing data to rescale covariance in the "col direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite

        .. note:    If the column axis is a link, only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        result = copy.copy( self ) #copy.deepcopy( self )
        if self.type==tokens.absoluteToken: return result

        # Make sure we have usable row data to rescale with
        if rowData is None: 
            rowData = self.ancestor.rowData.link.toPointwise_withLinearXYs()
        if not isinstance( rowData, XYsModule.XYs1d ):
            raise TypeError( 'rowData must be of type XYs1d, found %s' % type(rowData) )
        gRowData = rowData.group( self.matrix.axes[2].values )

        # Only generate the column rescaling if we need to
        if not self.matrix.axes[1].style=='link':
            if colData is None: 
                colData = self.ancestor.columnData.link.toPointwise_withLinearXYs()
            if not isinstance( colData, XYsModule.XYs1d ):
                raise TypeError( 'colData must be of type XYs1d, found %s' % type(colData) )
            gColData = colData.group( self.matrix.axes[1].values )
        else: gColData = gRowData
        
        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.array.constructArray()):
            newRow = []
            for j, cell in enumerate(row):
                newRow.append( gRowData[i][1]*gColData[j][1]*cell )
            newData.append( newRow )
        result.matrix.data = newData
        result.type=tokens.absoluteToken
        return result
        
    def toRelative( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a absolute covariance) using XYs1d rowData and colData
        to convert self into a relative covariance matrix.
        
        :param XYs1d rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite
        :param XYs1d colData: an XYs1d instance containing data to rescale covariance in the "col direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite

        .. note::   If the column axis is a link, only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        result = copy.copy( self ) #copy.deepcopy( self )
        if self.type==tokens.relativeToken: return result

        # Make sure we have usable row data to rescale with
        if rowData is None: 
            rowData = self.ancestor.rowData.link.toPointwise_withLinearXYs(1e-8,1e-8)
        if not isinstance( rowData, XYsModule.XYs1d ):
            raise TypeError( 'rowData must be of type XYs1d, found %s' % type(rowData) )
        gRowData = rowData.group( self.matrix.axes[2].values )

        # Only generate the column rescaling if we need to
        if not self.matrix.axes[1].style=='link':
            if colData is None: 
                colData = self.ancestor.columnData.link.toPointwise_withLinearXYs(1e-8,1e-8)
            if not isinstance( colData, XYsModule.XYs1d ):
                raise TypeError( 'colData must be of type XYs1d, found %s' % type(colData) )
            gColData = colData.group( self.matrix.axes[1].values )
        else: gColData = gRowData

        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.array.constructArray()):
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

        if self.isSymmetric() and info['checkUncLimits']:
            A = self.matrix.array.constructArray()
            relative = self.type == 'relative'
            if relative:
                varMin = info['minRelUnc']*info['minRelUnc']
                varMax = info['maxRelUnc']*info['maxRelUnc']
            for idx in range( A.shape[0] ):
                if not relative:
                    if info['theData'] is not None:
                        uncMin = info['minRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[2].grid[idx]+self.axes[2].grid[idx+1]) )
                        uncMax = info['maxRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[2].grid[idx]+self.axes[2].grid[idx+1]) )
                        varMin = uncMin*uncMin
                        varMax = uncMax*uncMax
                    else:
                        #warnings.append( "WARNING: can't check absolute uncertainties without data to compare to!\n" )
                        break
                if varMin <= A[idx,idx] and varMax >= A[idx,idx]: pass # unc is where is should be
                elif varMin >= A[idx,idx]:                         # uncertainty too small
                    warnings.append( warning.varianceTooSmall( idx, A[idx,idx], self ) )
                else:                                          # uncertainty too big
                    warnings.append( warning.varianceTooLarge( idx, A[idx,idx], self ) )

            # FIXME: is this the right place for eigenvalue checks? They used to live in fudge.core.math.matrix,
            # but that no longer exists
            vals, vecs = numpy.linalg.eigh( A )
            if min(vals) < info['negativeEigenTolerance']:
                warnings.append( warning.negativeEigenvalues( len(vals[vals<0]), min(vals), self ))
            minpos, maxpos = min(vals[vals>=0]),max(vals[vals>=0])
            # Check that the condition number of the matrix is reasonable
            if minpos/maxpos < info['eigenvalueRatioTolerance'] and A.size != 1:
                warnings.append( warning.badEigenvalueRatio( minpos/maxpos, self ) )

        return warnings
    
    def fix( self, **kw ): 
        """Fix uncertainty using the bounds passed into the fixer.
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        warnings = []
        if self.isSymmetric() and kw['fixUncLimits']:
            A = numpy.matrix( self.matrix.data )
            relative = self.type == 'relative'
            for idx in range( A.shape[0] ):
                eMin = self.axes[2].grid[idx]
                eMax = self.axes[2].grid[idx+1]
                if relative:
                    uncMin = kw['minRelUnc']
                    uncMax = kw['maxRelUnc']
                else:
                    uncMin = kw['theData'].getValue( 0.5*(eMax+eMin) )
                    uncMax = kw['theData'].getValue( 0.5*(eMax+eMin) )
                uncMin2 = uncMin*uncMin
                uncMax2 = uncMax*uncMax
                #eThresh = threshold.getValueAs( component.axes[0].units )
                if uncMin2 <= A[idx,idx] <= uncMax2: pass   # unc is where is should be
                elif uncMin2 >= A[idx,idx]:
                    if idx+1 < A.shape[0] and A[idx+1,idx+1] >= uncMin2: A[idx,idx] = A[idx+1,idx+1]
                    else: A[idx,idx] = uncMin2
                # else:                                            # above threshold and uncertainty out of bounds
                    # continue #skip this fix for now
                    # if uncMin2 > A[idx,idx]: # uncertainty too small
                        # print '    WARNING: bin', idx, 'uncertainty is too small!!!', '(', uncMin2, '>', A[idx,idx], ')'
                        # if A[idx, idx] == 0.0:
                            # for jdx in range( idx+1, A.shape[0] ):
                                # A[idx,jdx] = 0.0
                                # A[jdx,idx] = 0.0
                            # A[idx,idx] == uncMin2
                        # else:
                            # for jdx in range( idx+1, A.shape[0] ):
                                # A[idx,jdx] *= uncMin / math.sqrt( A[idx,idx] )
                                # A[jdx,idx] = A[idx,jdx]
                            # A[idx,idx] == uncMin2
                    # else:                # uncertainty too big
                        # print '    WARNING: bin', idx, 'uncertainty is too big!!!', '(', A[idx,idx], '>', uncMax2, ')'
                        # for jdx in range( idx+1, A.shape[0] ):
                            # A[idx,jdx] *= uncMax / math.sqrt( A[idx,idx] )
                            # A[jdx,idx] = A[idx,jdx]
                        # A[idx,idx] = uncMax2
            self.data = A.tolist()
        return warnings + self.matrix.fix( **kw )
        
    def group( self, groupBoundaries = ( None, None ), groupUnit = ( None, None ) ):
        '''
        Group the matrix in self

        :param groupBoundaries: a 2 element list containing the group boundaries for the rows 
                                and columns (respectively) of the covariance to be regrouped
                                rows go in the first element, columns in the second
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
        
        In the routine below, we abuse :py:class:`xData.XYs1d` to specify the functions 
        :math:`w_i(E)` and use the :py:func:`XYs1d.groupOneFunction()` method to perform the integrals to get
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
        # determine where to get the settings for the potentially mirrored second axis
        if self.matrix.axes[1].style == 'link': axis1index = 2
        else:                                   axis1index = 1
        
        # setup the old axes in a form we can (ab)use in the XYs1d class
        axes2_ = axesModule.axes( labelsUnits={1:( self.matrix.axes[2].label, self.matrix.axes[2].unit ),0:( 'dummy', '' )} )
        axes1_ = axesModule.axes( labelsUnits={1:( self.matrix.axes[axis1index].label, self.matrix.axes[axis1index].unit ),0:( 'dummy', '' )} )
        
        # define basis functions for the rows and columns
        basis2 = XYsModule.XYs1d( axes=axes2_, data=[ ( x, 0.0 ) for x in self.matrix.axes[2].values ], interpolation='flat' )
        basis1 = XYsModule.XYs1d( axes=axes1_, data=[ ( x, 0.0 ) for x in self.matrix.axes[axis1index].values ], interpolation='flat' )
        basis2 = basis2.convertAxisToUnit( 1, groupUnit[0] )
        basis1 = basis1.convertAxisToUnit( 1, groupUnit[1] )
    
        # build the regrouping matrices for the two bases
        w0 = []
        for idx in range( self.matrix.array.shape[0] ):
            basis2[idx] = ( basis2[idx][0], 1.0 )
            w0.append( basis2.groupOneFunction( groupBoundaries[0], norm = 'dx' ) )
            basis2[idx] = ( basis2[idx][0], 0.0 )
        w0 = numpy.mat( w0 )
        w1 = []
        for j in range( self.matrix.array.shape[1] ):
            basis1[j] = ( basis1[j][0], 1.0 )
            w1.append( basis1.groupOneFunction( groupBoundaries[1], norm = 'dx' ) )
            basis1[j] = ( basis1[j][0], 0.0 )
        w1 = numpy.mat( w1 )
                
        # set up the regrouped covariance matrix
        grouped = copy.copy( self )
        grouped.matrix.axes[2].data = groupBoundaries[0]
        grouped.matrix.axes[1].data = groupBoundaries[1]
        grouped.matrix.axes[2].unit = groupUnit[0]
        grouped.matrix.axes[1].unit = groupUnit[1]
        odata = numpy.mat( self.matrix.array.constructArray() )
        gdata = w0.T * odata * w1
        trigdata = gdata[numpy.tri(gdata.shape[0])==1.0].tolist()[0]
        grouped.matrix.array = arrayModule.full(shape=gdata.shape,data=trigdata,symmetry=arrayModule.symmetryLowerToken)
        return grouped

    def removeExtraZeros(self):
        """
        Remove all extra zeros from the underlying matrix
        :return:
        """
        matrix = self.matrix.array.constructArray()
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
        if self.matrix.axes[-2].style=='link':
            assert (rowStart,rowEnd) == (colStart,colEnd)
            self.matrix.axes[-1].__values = self.matrix.axes[-1].values[rowStart:rowEnd]
        else:
            self.matrix.axes[-1].__values = self.matrix.axes[-1].values[rowStart:rowEnd]
            self.matrix.axes[-2].__values = self.matrix.axes[-2].values[colStart:colEnd]
        self.matrix.array = arrayModule.full(shape=matrix.shape,data=matrix[numpy.tri(matrix.shape[0])==1.0].tolist(),symmetry=arrayModule.symmetryLowerToken)

    def getUncertaintyVector( self, theData=None, relative=True ):
        """ 
        Get an XYs1d object containing uncertainty for this matrix.
        Convert relative/absolute if requested (if so, must also pass central values as theData)

        Examples:

            - if the covariance matrix is relative and we want relative uncertainty vector, just do:
            
                > matrix.getUncertaintyVector()
                
            - if we want the absolute matrix instead:
            
                > matrix.getUncertaintyVector( theData=<XYs1d instance>, relative=False )
                
        """
        if not self.isSymmetric():
            raise ValueError("getUncertaintyVector only applies to symmetric matrices!")
        energies = list( self.matrix.axes[-1].values )
        diag = numpy.diagonal( self.matrix.array.constructArray() ).copy()
        diag[ diag < 0.0 ] = 0.0
        diag = list( numpy.sqrt( diag ) )
        diag.append( diag[-1] )                             # point corresponding to final energy bin
        yunit = self.matrix.axes[0].unit
        if yunit != '': # get square root of the unit
            yunit = PQU.PQU(1,yunit).sqrt().getUnitSymbol()
        axes_ = axesModule.axes(
                                 labelsUnits={1:('energy_in',self.matrix.axes[2].unit),
                                              0:('uncertainty',yunit)} )
        uncert = XYsModule.XYs1d( zip(energies,copy.deepcopy(diag)), axes=axes_, interpolation='flat',
                                accuracy=0.0001, )    # what should accuracy be set to?
        uncert = uncert.changeInterpolation('lin-lin',None,1e-8,1e-8)

        # do we need to convert absolute->relative or vice versa?
        if (relative and self.type==tokens.absoluteToken) or (not relative and self.type==tokens.relativeToken):
            if theData is None:
                theData = self.ancestor.rowData.link.toPointwise_withLinearXYs(1e-8,1e-8)
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
        """

        :param title:
        :param scalelabel:
        :param xlim:
        :param ylim:
        :param xlog:
        :param ylog:
        :return:
        """
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
        if title is None: title = str( self.toXLink() )
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

        if xlim is None:    ax.set_xlim( x[0], x[-1] )
        else:               ax.set_xlim( xlim[0], xlim[1] )
        if ylim is None:    ax.set_ylim( y[0], y[-1] )
        else:               ax.set_ylim( ylim[0], ylim[1] )
        if xlog: ax.set_xscale( 'log' )
        if ylog: ax.set_yscale( 'log' )

        xlabel = self.axes[0].label + ' (' + self.axes[0].unit +')'
        if self.axes[1].mirrorOtherAxis:    ylabel = xlabel 
        else:                               ylabel = self.axes[1].label + ' (' + self.axes[1].unit +')'

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        cbar = plt.colorbar(qc)
        if scalelabel is not None: cbar.set_label(scalelabel)
        else: cbar.set_label(str(self.type)+' covariance ('+str(self.axes[2].unit)+')')
        plt.show()

    def toXMLList( self, indent = '', **kwargs ) :
        """

        :param indent:
        :param kwargs:
        :return:
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s' % ( indent, self.moniker ) ]
        if self.label is not None: xmlString[0] += ' label="%s"' % self.label
        xmlString[0] += ' type="%s"' % self.type
        if( self.energyBounds ) :
            lowerBound = self.energyBounds[0].toString( keepPeriod = False )
            upperBound = self.energyBounds[1].toString( keepPeriod = False )
            xmlString[0] += ( ' lowerBound="%s" upperBound="%s"' % ( lowerBound, upperBound ) )
        if self.ENDFconversionFlag: xmlString[0] += (
                ' ENDFconversionFlag="%s"' % self.ENDFconversionFlag )
        xmlString[0] += '>'
        xmlString += self.matrix.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode(element, xPath, linkData):
        """Translate <covarianceMatrix> element from xml into python class."""

        xPath.append( element.tag )
        energyBounds = None
        if element.get('energyBounds'):
            raise NotImplementedError("need to add support for incident energy bounds")
        matrix_ = griddedModule.gridded.parseXMLNode( element[0], xPath, linkData )
        CM = covarianceMatrix( label = element.get( 'label' ), type=element.get('type'), matrix=matrix_,
                energyBounds=energyBounds, ENDFconversionFlag=element.get("ENDFconversionFlag") )
        xPath.pop()
        return CM
