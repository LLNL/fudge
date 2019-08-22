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

from xData.ancestry import ancestry
from xData import gridded as griddedModule
from xData import axes as axesModule
from xData import link as linkModule
import copy

__metaclass__ = type

class mixedForm( ancestry ):
    """
    Covariance for a single quantity, stored as several separate matrices that must be summed together.
    In general, the energy bounds for these matrices can overlap (unlike regions1d cross section data). 
    """

    moniker = 'mixed'

    def __init__(self, label = None, components=None):
        ancestry.__init__( self )
        self.__label = label
        self.components = components or [] #: a Python list containing instances of ``mixedForm``, ``summedCovariance``, and ``covarianceMatrix``
        for c in self.components:
            c.setAncestor(self)

    def __getitem__(self, idx):
        return self.components[idx]

    def __len__(self):
        return len(self.components)

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, label ) :

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def getRowBounds(self,unit=None):
        """
        Get the bounds of the row. If unit is specified, return the bounds in that unit.
        Otherwise, take unit from first sub-matrix.
        """
        if unit is None: unit = self[0].matrix.axes[-1].unit
        allbounds = [ c.getRowBounds(unit) for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))
    
    def getColumnBounds(self,unit=None):
        """
        Get the bounds of the column.  If unit is specified, return the bounds in that unit.
        Otherwise, take unit from the first sub-matrix.
        """
        if unit is None: unit = self[0].matrix.axes[-1].unit
        allbounds = [ c.getColumnBounds(unit) for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))

    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        for component in self.components:
            componentWarnings = component.check( info )
            if componentWarnings:
                warnings.append( warning.context('Component %s' % component.label, componentWarnings) )
        
        #if warnings:
        #    warnings = [warning.context('Section "%s": %s' % (self.id, form), warnings)]
        return warnings

    def convertUnits( self, unitMap ):

        for component in self: component.convertUnits( unitMap )
        
    def fix( self, **kw ): 
        warnings = []
        for comp in self.components: warnings += comp.fix( **kw )
        return warnings

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def addComponent(self, covariance):
        ''':param covariance: an instance of covariance (or inherited class)'''
        covariance.setAncestor(self)
        self.components.append(covariance)
        
    def getMatchingComponent(self,rowBounds=None,columnBounds=None):
        """

        :param rowBounds:
        :param columnBounds:
        :return:
        """
        gotRow = False
        for comp in self.components:
            if comp.getRowBounds() == rowBounds:
                gotRow = True
                if comp.getColumnBounds() == columnBounds: return comp
        if not gotRow: raise ValueError( "No component with row bounds matching %s found" % str(rowBounds) )
        raise ValueError( "No component with column bounds matching %s found" % str(columnBounds) )

    def makeSafeBounds(self):
        '''
        Go through all the components and make sure the bounds don't overlap.  If they do, it is likely a bug.
        '''
        import fudge.gnd.covariances.base as base
        itWorked = True
        for ic in range(len(self)-1):
            c0 = self.components[ic]
            if c0.getRowBounds() != c0.getColumnBounds():
                raise ValueError( "All components must have their row and column covarianceAxes matching.")
            c1 = self.components[ic+1]
            if c1.getRowBounds() != c1.getColumnBounds():
                raise ValueError( "All components must have their row and column covarianceAxes matching.")
            itWorked = c0.getRowBounds()[1] <= c1.getRowBounds()[0]
            if not itWorked: # uh oh
                # Try to fix c0
                if isinstance(c0,base.covarianceMatrix):
                    c0.removeExtraZeros()
                    itWorked = c0.getRowBounds()[1] <= c1.getRowBounds()[0]
            if not itWorked: # double uh oh
                # Try to fix c1
                if isinstance(c1,base.covarianceMatrix):
                    c1.removeExtraZeros()
                    itWorked = c0.getRowBounds()[1] <= c1.getRowBounds()[0]
        if not itWorked:
            raise ValueError("Bounds between components %i and %i overlap: %s" %(ic,ic+1, str( [c.getRowBounds() for c in self.components])))
                
    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combines all subsections into single uncertainty vector, converting to relative if requested.

        :returns: an XYs1d instance
        """
        return self.toCovarianceMatrix().getUncertaintyVector(theData=theData,relative=relative)

    def getCorrelationMatrix(self):
        return self.toCovarianceMatrix().getCorrelationMatrix()

    def toCovarianceMatrix( self, label="composed" ):
        """
        Sum all parts together to build a single matrix.
        """
        if len( self.components ) == 1: return self.components[0].toCovarianceMatrix()
        import fudge.gnd.covariances.base as base
        import numpy
        import xData.values as valuesModule
        import xData.array as arrayModule

        # Set up common data using first component
        firstCovMtx = self.components[0].toCovarianceMatrix()
        if not isinstance( firstCovMtx, base.covarianceMatrix):
            raise TypeError("Shoudd have gotten base.covarianceMatrix, instead got %s"%str(firstCovMtx.__class__))
        commonRowAxis = firstCovMtx.matrix.axes[2].copy([])#FIXME: unresolvedLinks are still unresolved!
        if firstCovMtx.matrix.axes[1].style=='link':
            commonColAxis = firstCovMtx.matrix.axes[2].copy([])#FIXME: unresolvedLinks are still unresolved!
        else:
            commonColAxis = firstCovMtx.matrix.axes[1].copy([])#FIXME: unresolvedLinks are still unresolved!
        commonMatrixAxis = firstCovMtx.matrix.axes[0] .copy([])#FIXME: unresolvedLinks are still unresolved!
        commonType = firstCovMtx.type

        # We need all the covariances to be either absolute or relative
        def make_common_type(cm):
            if commonType == 'relative': return cm.toRelative()
            else: return cm.toAbsolute()

        # We're going to have to merge grids, so we'll need this function to do the dirty work
        def add_values(v1,v2):
            v=set()
            v.update(v1.values)
            v.update(v2.values)
            return valuesModule.values(sorted(v))

        # First pass through components is to collect bins to set up the common grid + do assorted checking
        for c in self.components[1:]:
            cc = make_common_type(c.toCovarianceMatrix()) # a little recursion to take care of nested covariances
            if cc.type != commonType:
                raise ValueError( "Incompatible types in %s: %s vs. %s" % (self.__class__, commonType, cc.type) )
            if cc.matrix.axes[0].unit !=  commonMatrixAxis.unit: raise ValueError("covariance matrix components with different units?!? %s vs. %s"%(cc.matrix.axes[0].unit, commonMatrixAxis.unit))
            if cc.matrix.axes[1].style != 'link': cc.matrix.axes[1].convertToUnit(commonColAxis.unit)
            cc.matrix.axes[2].convertToUnit(commonRowAxis.unit)
            commonRowAxis.values.values = add_values(commonRowAxis.values, cc.matrix.axes[2].values)
            if cc.matrix.axes[1].style == 'link': commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[2].values)
            else:                                 commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[1].values)

        # Now sum up the components
        commonMatrix = numpy.mat( firstCovMtx.group( ( commonRowAxis.values.values, commonColAxis.values.values ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray() )
        for c in self.components[1:]:
            cc = make_common_type(c.toCovarianceMatrix()) # a little recursion to take care of nested covariances
            commonMatrix += numpy.mat( cc.group( ( commonRowAxis.values.values, commonColAxis.values.values ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray() )

        # Now create the instance of the resulting covarianceMatrix
        if all( [component.toCovarianceMatrix().matrix.axes[1].style == 'link' for component in self.components ] ):  commonColAxis = self.components[0].toCovarianceMatrix().matrix.axes[1].copy([])#FIXME: unresolvedLinks are still unresolved!
        newAxes=axesModule.axes(
                labelsUnits={0 : (commonMatrixAxis.label, commonMatrixAxis.unit),
                             1 : (commonColAxis.label, commonColAxis.unit),
                             2 : (commonRowAxis.label, commonRowAxis.unit)} )
        newAxes[2] = axesModule.grid(commonRowAxis.label, commonRowAxis.index, commonRowAxis.unit,
                                  style=axesModule.boundariesGridToken,
                                  values=commonRowAxis.values)
        newAxes[1] = axesModule.grid(commonColAxis.label, commonColAxis.index, commonColAxis.unit,
                              style=axesModule.linkGridToken,
                              values=linkModule.link(link=commonRowAxis.values, relative=True))
        newAxes[0] = axesModule.axis(commonMatrixAxis.label, commonMatrixAxis.index, commonMatrixAxis.unit)
        trigdata = commonMatrix[numpy.tri(commonMatrix.shape[0])==1.0].tolist()[0]
        gridded = griddedModule.gridded2d( axes=newAxes, array=arrayModule.full(shape=commonMatrix.shape, data=trigdata, symmetry=arrayModule.symmetryLowerToken))
        newCov = base.covarianceMatrix( label=label, type=commonType, matrix=gridded )
        newCov.setAncestor(self.ancestor)
        return newCov

    def toAbsolute( self, rowData=None, colData=None ):
        '''
        Rescales self (if it is a relative covariance) using XYs1d rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        result = copy.copy( self )
        result.components = []
        for c in self.components:
#            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toAbsolute() call with the correct row or column data
                result.addComponent( c.toCovarianceMatrix().toAbsolute( rowData, colData ) )
#            else: result.components.append( c.toAbsolute( rowData, colData ) )
        return result
        
    def toRelative( self, label="composed", rowData=None, colData=None ):
        '''
        Rescales self (if it is a absolute covariance) using XYs1d rowData and colData
        to convert self into a relative covariance matrix.
        
        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        result = copy.copy( self )
        result.components = []
        for c in self.components:
#            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toRelative() call with the correct row or column data
                result.addComponent( c.toCovarianceMatrix().toRelative( rowData, colData ) )
#            else: result.components.append( c.toRelative( rowData, colData ) )
        return result

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        for covariance in self.components:
            xmlString += covariance.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        """Translate <mixed> element from xml."""
        import fudge.gnd.covariances.base as base
        import fudge.gnd.covariances.summed as summed

        xPath.append( element.tag )
        mixed_ = cls( label = element.get( "label" ) )
        for child in element:
            formClass = {
                    base.covarianceMatrix.moniker: base.covarianceMatrix,
                    summed.summedCovariance.moniker: summed.summedCovariance,
                    }.get( child.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % child.tag)
            mixed_.addComponent( formClass.parseXMLNode( child, xPath, linkData ) )
        for i,form in enumerate( mixed_.components ): form.index = i
        xPath.pop()
        return mixed_
