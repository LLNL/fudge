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

from xData.ancestry import ancestry
from xData import gridded as griddedModule
from xData import axes as axesModule
import copy

__metaclass__ = type

class mixedForm( ancestry ):
    """
    Covariance for a single quantity, stored as several separate matrices that must be summed together.
    In general, the energy bounds for these matrices can overlap (unlike piecewise cross section data). 
    """

    moniker = 'mixed'

    def __init__(self, label = None, components=None):
        ancestry.__init__( self )
        self.__label = label
        self.components = components or [] #: a Python list containing instances of ``mixedForm``, ``summedCovariance``, and ``covarianceMatrix``

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
        '''Get the bounds of the row.  If unit is specified, return the bounds in that unit.'''
        allbounds = [ c.getRowBounds() for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))
    
    def getColumnBounds(self,unit=None):
        '''Get the bounds of the column.  If unit is specified, return the bounds in that unit.'''
        allbounds = [ c.getColumnBounds() for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))

    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        for component in self.components:
            componentWarnings = component.check( info )
            if componentWarnings:
                warnings.append( warning.context('Component %i' % component.index, componentWarnings) )
        
        #if warnings:
        #    warnings = [warning.context('Section "%s": %s' % (self.id, form), warnings)]
        return warnings
        
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
        gotRow = False
        for comp in self.components:
            if comp.getRowBounds() == rowBounds:
                gotRow = True
                if comp.getColumnBounds() == columnBounds: return comp
        if not gotRow: raise ValueError( "No component with row bounds matching %s found" % str(rowBounds) )
        raise ValueError( "No component with column bounds matching %s found" % str(columnBounds) )

    def shrinkToBounds(self,bounds):
        '''
        Return the subset of self's components whose row bounds fit within requested bounds.
        Only works on mixedForms that are made up of parts where the column's covarianceAxis is a mirror of the row's.
        '''
        import fudge.gnd.covariances.base as base
        def inRange( thisBounds, otherBounds ):
            return otherBounds[0] >= thisBounds[0] and otherBounds[1] <= thisBounds[1]
        newMixed = mixedForm()
        for c in self.components:
            if c.getRowBounds() != c.getColumnBounds(): raise ValueError( "All components must have their row and column covarianceAxes matching.")
            if isinstance(c,base.covarianceMatrix): 
                c = copy.copy(c)
                c.removeExtraZeros()
            if inRange(bounds,c.getRowBounds()): newMixed.addComponent(c)
        return newMixed

    def makeSafeBounds(self):
        '''
        Go through all the components and make sure the bounds don't overlap.  If they do, it is likely a bug.
        '''
        import fudge.gnd.covariances.base as base
        for ic in range(len(self)-1):
            itWorked = True
            c0 = self.components[ic]
            if c0.getRowBounds() != c0.getColumnBounds(): raise ValueError( "All components must have their row and column covarianceAxes matching.")
            c1 = self.components[ic+1]
            if c1.getRowBounds() != c1.getColumnBounds(): raise ValueError( "All components must have their row and column covarianceAxes matching.")
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
            print [c.getRowBounds() for c in self.components]
            raise ValueError("Bounds between components %i and %i overlap!" %(ic,ic+1))
                
    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combines all subsections into single uncertainty vector, converting to relative if requested.

        :returns: an XYs instance
        """
        return self.toCovarianceMatrix().getUncertaintyVector(theData=theData,relative=relative)

    def toCovarianceMatrix( self ):
        """
        Sum all parts together to build a single matrix.
        FIXME: currently breaks if the 'mixed' section contains a mixture of absolute/relative/correlation matrices.
        """
        if len( self.components ) == 1: return self.components[0].toCovarianceMatrix()
        import fudge.gnd.covariances.base as base
        import numpy
        
        # set up common data using first component
        firstCovMtx = self.components[0].toCovarianceMatrix()
        commonRowAxis = copy.copy( firstCovMtx.matrix.axes[2] )
        if firstCovMtx.matrix.axes[1].gridStyle=='link':
            commonColAxis = copy.copy( firstCovMtx.matrix.axes[2] )
        else:
            commonColAxis = copy.copy( firstCovMtx.matrix.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.matrix.axes[0] )
        commonType = firstCovMtx.type
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
            if cc.type != commonType:
                raise ValueError( "Incompatible types in %s: %s vs. %s" % (self.__class__, commonType, cc.type) )
            cc.convertAxesToUnits( ( commonRowAxis.unit, commonColAxis.unit, commonMatrixAxis.unit ) )
            commonRowAxis.data = commonRowAxis.data + cc.axes[0].data
            if cc.matrix.axes[1].mirrorOtherAxis: commonColAxis.data = commonColAxis.data + cc.matrix.axes[2].data
            else:                                 commonColAxis.data = commonColAxis.data + cc.matrix.axes[1].data
        commonRowAxis.data = sorted( set( commonRowAxis.data ) )
        commonColAxis.data = sorted( set( commonColAxis.data ) )
        
        # now sum up the components
        commonMatrix = numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        if all( [matrix.axes[1].mirrorOtherAxis for matrix in self.components] ):
            commonColAxis.mirrorOtherAxis = True
            commonColAxis.data = []
        gridded = griddedModule.gridded( axes=axesModule.axes( [commonRowAxis, commonColAxis, commonMatrixAxis] ),
                                         array=commonMatrix.tolist(), copyArray=False, label = 'unified' )
        return base.covarianceMatrix( type=commonType, matrix=gridded )

    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
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
                result.components.append( c.toCovarianceMatrix().toAbsolute( rowData, colData ) )
#            else: result.components.append( c.toAbsolute( rowData, colData ) )
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
        result = copy.copy( self )
        result.components = []
        for c in self.components:
#            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toRelative() call with the correct row or column data
                result.components.append( c.toCovarianceMatrix().toRelative( rowData, colData ) )
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
