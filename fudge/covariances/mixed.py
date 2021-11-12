# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from . import base, shortRangeSelfScalingVariance
from xData.ancestry import ancestry
from xData import gridded as griddedModule
from xData import axes as axesModule
from xData import link as linkModule

__metaclass__ = type

class mixedForm( ancestry, base.covariance ):
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
            c.setAncestor(self, attribute='label')

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

    @property
    def domainUnit( self ):
        return self[0].domainUnit

    def rowBounds(self, unit=None):
        """
        Get the bounds of the row. If unit is specified, return the bounds in that unit.
        Otherwise, take unit from first sub-matrix.
        """
        if unit is None: unit = self.domainUnit
        allbounds = [ c.rowBounds(unit) for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))
    
    def columnBounds(self, unit=None):
        """
        Get the bounds of the column.  If unit is specified, return the bounds in that unit.
        Otherwise, take unit from the first sub-matrix.
        """
        if unit is None: unit = self.domainUnit
        allbounds = [ c.columnBounds(unit) for c in self.components ]
        return (min([x[0] for x in allbounds]),max([x[1] for x in allbounds]))

    def check( self, info ): 
        from fudge import warning
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

    def group( self, groupBoundaries = ( None, None ), groupUnit = ( None, None ) ):
        return self.toCovarianceMatrix().group( groupBoundaries=groupBoundaries, groupUnit=groupUnit )

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def addComponent(self, covariance):
        """:param covariance: an instance of covariance (or inherited class)"""
        covariance.setAncestor(self, attribute='label')
        self.components.append(covariance)
        
    def getMatchingComponent(self,rowBounds=None,columnBounds=None):
        """

        :param rowBounds:
        :param columnBounds:
        :return:
        """
        gotRow = False
        for comp in self.components:
            if comp.rowBounds() == rowBounds:
                gotRow = True
                if comp.columnBounds() == columnBounds: return comp
        if not gotRow: raise ValueError( "No component with row bounds matching %s found" % str(rowBounds) )
        raise ValueError( "No component with column bounds matching %s found" % str(columnBounds) )

    def makeSafeBounds(self):
        """
        Go through all the components and make sure the bounds don't overlap.  If they do, it is likely a bug.
        """
        from fudge.covariances import covarianceMatrix
        itWorked = True
        for ic in range(len(self)-1):
            c0 = self.components[ic]
            if c0.rowBounds() != c0.columnBounds():
                raise ValueError( "All components must have their row and column covarianceAxes matching.")
            c1 = self.components[ic+1]
            if c1.rowBounds() != c1.columnBounds():
                raise ValueError( "All components must have their row and column covarianceAxes matching.")
            itWorked = c0.rowBounds()[1] <= c1.rowBounds()[0]
            if not itWorked: # uh oh
                # Try to fix c0
                if isinstance(c0,covarianceMatrix.covarianceMatrix):
                    c0.removeExtraZeros()
                    itWorked = c0.rowBounds()[1] <= c1.rowBounds()[0]
            if not itWorked: # double uh oh
                # Try to fix c1
                if isinstance(c1,covarianceMatrix.covarianceMatrix):
                    c1.removeExtraZeros()
                    itWorked = c0.rowBounds()[1] <= c1.rowBounds()[0]
        if not itWorked:
            raise ValueError("Bounds between components %i and %i overlap: %s" %
                             (ic,ic+1, str( [c.rowBounds() for c in self.components])))

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
        from .covarianceMatrix import covarianceMatrix
        from .shortRangeSelfScalingVariance import shortRangeSelfScalingVariance
        import numpy
        import xData.values as valuesModule
        import xData.xDataArray as arrayModule

        # Set up common data using first component
        summands = [term for term in self.components if not isinstance(term, shortRangeSelfScalingVariance)]
        firstCovMtx = summands[0].toCovarianceMatrix()
        commonRowAxis = firstCovMtx.matrix.axes[2].copy([])#FIXME: unresolvedLinks are still unresolved!
        if firstCovMtx.matrix.axes[1].style=='link':
            commonColAxis = firstCovMtx.matrix.axes[2].copy([])#FIXME: unresolvedLinks are still unresolved!
        else:
            commonColAxis = firstCovMtx.matrix.axes[1].copy([])#FIXME: unresolvedLinks are still unresolved!
        commonMatrixAxis = firstCovMtx.matrix.axes[0] .copy([])#FIXME: unresolvedLinks are still unresolved!
        commonType = firstCovMtx.type

        # We need all the covariances to be either absolute or relative
        def make_common_type(cm):
            cm2 = cm.toCovarianceMatrix()
            cm2.setAncestor(cm.ancestor)    # ensure cm2 has access to rowData, etc.
            if commonType == 'relative': return cm.toRelative()
            else: return cm.toAbsolute()

        # We're going to have to merge grids, so we'll need this function to do the dirty work
        def add_values(v1,v2):
            v=set()
            v.update(v1.values)
            v.update(v2.values)
            return valuesModule.values(sorted(v))

        # First pass through components is to collect bins to set up the common grid + do assorted checking
        for c in summands[1:]:
            cc = make_common_type(c)
            if cc.type != commonType:
                raise ValueError( "Incompatible types in %s: %s vs. %s" % (self.__class__, commonType, cc.type) )
            if cc.matrix.axes[0].unit !=  commonMatrixAxis.unit:
                raise ValueError("covariance matrix components with different units?!? %s vs. %s" %
                                 (cc.matrix.axes[0].unit, commonMatrixAxis.unit))
            if cc.matrix.axes[1].style != 'link': cc.matrix.axes[1].convertToUnit(commonColAxis.unit)
            cc.matrix.axes[2].convertToUnit(commonRowAxis.unit)
            commonRowAxis.values.values = add_values(commonRowAxis.values, cc.matrix.axes[2].values)
            if cc.matrix.axes[1].style == 'link':
                commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[2].values)
            else:
                commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[1].values)

        # Now sum up the components
        commonMatrix = numpy.mat( firstCovMtx.group( ( commonRowAxis.values.values, commonColAxis.values.values ),
            ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray() )
        for c in summands[1:]:
            cc = make_common_type(c)
            commonMatrix += numpy.mat( cc.group( ( commonRowAxis.values.values, commonColAxis.values.values ),
                ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray() )

        # Now create the instance of the resulting covarianceMatrix
        if all( [component.toCovarianceMatrix().matrix.axes[1].style == 'link' for component in summands ] ):
            commonColAxis = summands[0].toCovarianceMatrix().matrix.axes[1].copy([])#FIXME: unresolvedLinks are still unresolved!
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
        gridded = griddedModule.gridded2d( axes=newAxes, array=arrayModule.full(shape=commonMatrix.shape,
            data=trigdata, symmetry=arrayModule.symmetryLowerToken))
        newCov = covarianceMatrix( label=label, type=commonType, matrix=gridded )
        newCov.setAncestor(self.ancestor)
        return newCov

    '''
    def toAbsolute( self, rowData=None, colData=None ):
        """
        Rescales self (if it is a relative covariance) using XYs1d rowData and colData
        to convert self into an absolute covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"

        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already
                    absolute because this will throw an error.

        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        """
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
        """
        Rescales self (if it is a absolute covariance) using XYs1d rowData and colData
        to convert self into a relative covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"

        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already
                    relative because this will throw an error.

        :returns: a copy of self, but rescaled and with the type set to relativeToken
        """
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
    '''

    def toAbsolute(self, rowData=None, colData=None):
        """
        Rescales self (if it is a relative covariance) using XYs1d rowData and columnData
        to convert self into an absolute covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "column direction"

        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.

        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        """
        if rowData is None:
            rowData = self.findAttributeInAncestry('rowData').link.toPointwise_withLinearXYs(lowerEps=1e-8)
        if colData is None and self.findAttributeInAncestry('columnData') is not None:
            colData = self.findAttributeInAncestry('colData').link.toPointwise_withLinearXYs(lowerEps=1e-8)
        return self.toCovarianceMatrix().toAbsolute(rowData=rowData, colData=colData)

    def toRelative(self, rowData=None, colData=None):
        """
        Rescales self (if it is a absolute covariance) using XYs1d rowData and columnData
        to convert self into a relative covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs1d instance containing data to rescale covariance in the "column direction"

        .. note::   If the column axis is a link to row axis, only rowData is needed.

        :returns: a copy of self, but rescaled and with the type set to relativeToken
        """
        if rowData is None:
            rowData = self.findAttributeInAncestry('rowData').link.toPointwise_withLinearXYs(lowerEps=1e-8)
        if colData is None and self.findAttributeInAncestry('columnData') is not None:
            colData = self.findAttributeInAncestry('columnData').link.toPointwise_withLinearXYs(lowerEps=1e-8)
        return self.toCovarianceMatrix().toRelative(rowData=rowData, colData=colData)

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
        from fudge.covariances import covarianceMatrix, summed

        xPath.append( element.tag )
        mixed_ = cls( label = element.get( "label" ) )
        for child in element:
            formClass = {
                    covarianceMatrix.covarianceMatrix.moniker: covarianceMatrix.covarianceMatrix,
                    summed.summedCovariance.moniker: summed.summedCovariance,
                    shortRangeSelfScalingVariance.shortRangeSelfScalingVariance.moniker:
                        shortRangeSelfScalingVariance.shortRangeSelfScalingVariance,
                    }.get( child.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % child.tag)
            mixed_.addComponent( formClass.parseXMLNode( child, xPath, linkData ) )
        for i,form in enumerate( mixed_.components ): form.index = i
        xPath.pop()
        return mixed_
