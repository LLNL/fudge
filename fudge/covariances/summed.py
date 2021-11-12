# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from . import base
from xData.ancestry import ancestry
from xData import link as linkModule
from pqu import PQU

__metaclass__ = type

class summedCovariance( ancestry, base.covariance ):
    """ 
    Covariance matrix stored as sum/difference of matrices from other reactions,
    valid over a specified domain.
    Note: summed matrices may sum over other summed matrices, which opens up the possibility of
    a recursively-defined covariance. Recursion is avoided so long as the domains don't overlap.

    Each summand has a coefficient 'c' that scales the associated covariance when summing together:

        $$\sum_{0}^{N} c * Cov * c.T$$

    In current data libraries the coefficients are always +/-1.
    """

    moniker = 'sum'

    def __init__( self, label, domainMin, domainMax, domainUnit='eV', summands=None ):
        """
        :param label: label of associated style
        :param domainMin: lower domain limit for the sum
        :param domainMax: upper domain limit for the sum
        :param domainUnit: unit for domain min/max
        :param summands: list of 'summand' instances
        """
        ancestry.__init__( self )
        self.label = label
        self.__domainMin = domainMin
        self.__domainMax = domainMax
        self.__domainUnit = domainUnit
        self.summands = summands or []

    def __getitem__(self, idx):
        return self.summands[idx]

    def __len__(self):
        return len(self.summands)

    @property
    def label( self ) :

        return self.__label

    @label.setter
    def label( self, value ) :

        if value is not None:
            if not isinstance( value, str ): raise TypeError( 'label must be a string' )
        self.__label = value

    @property
    def domainMin( self ):
        return self.__domainMin

    @property
    def domainMax( self ):
        return self.__domainMax

    @property
    def domainUnit( self ):
        return self.__domainUnit

    def check( self, info ): return []

    def convertUnits( self, unitMap ):

        # don't descend into the list of summands: those should be converted independently
        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            self.__domainMin = PQU.PQU(self.domainMin, self.domainUnit).getValueAs(newUnit)
            self.__domainMax = PQU.PQU(self.domainMax, self.domainUnit).getValueAs(newUnit)
            self.__domainUnit = newUnit

    def fix( self, **kw ): return []

    def group( self, groupBoundaries = ( None, None ), groupUnit = ( None, None ) ):
        return self.toCovarianceMatrix().group( groupBoundaries=groupBoundaries, groupUnit=groupUnit )

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def rowBounds(self,unit=None):
        """Get the bounds of the row.  If unit is specified, return the bounds in that unit."""
        if unit is None:
            return (self.domainMin, self.domainMax)
        else:
            factor = PQU.PQU(1, self.domainUnit).getValueAs(unit)
            return (self.domainMin * factor, self.domainMax * factor)
    
    def columnBounds(self,unit=None):
        """Get the bounds of the column.  If unit is specified, return the bounds in that unit."""
        if unit is None:
            return (self.domainMin, self.domainMax)
        else:
            factor = PQU.PQU(1, self.domainUnit).getValueAs(unit)
            return (self.domainMin * factor, self.domainMax * factor)

    def toCorrelationMatrix( self ):
        return self.toCovarianceMatrix().toCorrelationMatrix()

    def toCovarianceMatrix( self, label="composed" ):
        """
        Sum the parts to construct the covariance matrix on the domain of self.
        Note: each part must be converted to absolute before summing, since the sum is over different reactions.

        :param label: label attached to the new covarianceMatrix
        :return: covarianceMatrix instance, on the union grid of all summands over the domain of self.
        """
        import numpy
        from .mixed import mixedForm
        from .covarianceMatrix import covarianceMatrix
        from .shortRangeSelfScalingVariance import shortRangeSelfScalingVariance
        from xData import values as valuesModule
        from xData import axes as axesModule
        from xData import xDataArray as arrayModule
        from xData import gridded as griddedModule

        domain = (self.domainMin, self.domainMax)

        if len( self.summands ) == 1:
            return self.summands[0].link.toCovarianceMatrix().domainSlice(rowDomainBounds=domain)

        def getAbsoluteMatrix(p):
            """
            helper function: for each summand 'p', first check if it overlaps
            with self's domain. If not return None, otherwise convert it to an absolute matrix
            on self's domain and return that
            """
            cm = p.link

            if isinstance(cm, (covarianceMatrix, summedCovariance)):
                cm = cm.toCovarianceMatrix()
            elif isinstance(cm,mixedForm):
                newMixed = mixedForm()
                for ic, c in enumerate(cm.components):
                    if c.rowBounds() != c.columnBounds():
                        raise ValueError("All components must have their row and column covarianceAxes matching.")
                    if (c.rowBounds()[-1] <= self.domainMin) or (c.rowBounds()[0] >= self.domainMax):
                        continue    # this part of 'mixed' matrix doesn't contribute to the sum
                    c_copy = c.copy()
                    # prune zero rows/columns covarianceMatrices, just in case
                    if isinstance(c_copy, covarianceMatrix):
                        if not c_copy.matrix.array.constructArray().any(): continue  # matrix is all zero
                        c_copy.removeExtraZeros()
                    elif isinstance(c_copy, mixedForm):
                        c_copy.makeSafeBounds()
                    newMixed.addComponent(c_copy)
                cm = newMixed.toCovarianceMatrix()

            # reduce the domain if necessary to match the summed matrix domain:
            cm = cm.domainSlice( self.rowBounds(), self.columnBounds() )

            if cm.type == 'absolute':
                return cm
            else:
                # FIXME: following only works for cross sections/multiplicities, not spectra
                meanValue = p.link.ancestor.rowData.link.toPointwise_withLinearXYs(lowerEps=1e-8)
                return cm.toAbsolute(rowData=meanValue)

        # Set up common data using first element in summands
        summands = [summand for summand in self.summands if not isinstance(summand.link, shortRangeSelfScalingVariance)]
        firstCovMtx = getAbsoluteMatrix(summands[0])
        commonRowAxis = firstCovMtx.matrix.axes[2].copy([])
        if firstCovMtx.matrix.axes[1].style=='link':
            commonColAxis = firstCovMtx.matrix.axes[2].copy([])
        else:
            commonColAxis = firstCovMtx.matrix.axes[1].copy([])
        commonMatrixAxis = firstCovMtx.matrix.axes[0].copy([])

        def add_values(v1,v2):
            # helper function for merging grids
            v=set()
            v.update(v1.values)
            v.update(v2.values)
            return valuesModule.values(sorted(v))

        # First pass through components is to collect bins to set up the common grid + do assorted checking
        for p in summands[1:]:
            cc = getAbsoluteMatrix( p )     # recursion to take care of nested covariances
            cc.matrix.axes[0].unit = commonMatrixAxis.unit
            cc.matrix.axes[1].convertToUnit(commonColAxis.unit)
            cc.matrix.axes[2].convertToUnit(commonRowAxis.unit)
            commonRowAxis.values.values = add_values(commonRowAxis.values, cc.matrix.axes[2].values)
            if cc.matrix.axes[1].style==axesModule.linkGridToken:
                commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[2].values)
            else:
                commonColAxis.values.values = add_values(commonColAxis.values, cc.matrix.axes[1].values)

        # Now sum up the components on common grid
        commonMatrix = summands[0].coefficient**2 * firstCovMtx.group(
                ( commonRowAxis.values.values, commonColAxis.values.values ),
                ( commonRowAxis.unit, commonColAxis.unit )
            ).matrix.array.constructArray()
        for p in summands[1:]:
            cc = getAbsoluteMatrix( p ) # recursion to take care of nested covariances
            commonMatrix += p.coefficient**2 * cc.group(
                    ( commonRowAxis.values.values, commonColAxis.values.values ),
                    ( commonRowAxis.unit, commonColAxis.unit )
                ).matrix.array.constructArray()
        
        # Now create the instance of the resulting covarianceMatrix
        newAxes = axesModule.axes(
            labelsUnits={0: (commonMatrixAxis.label, commonMatrixAxis.unit),
                         1: (commonColAxis.label, commonColAxis.unit),
                         2: (commonRowAxis.label, commonRowAxis.unit)})

        newAxes[2] = axesModule.grid(commonRowAxis.label, commonRowAxis.index, commonRowAxis.unit,
                                 style=axesModule.boundariesGridToken,
                                 values=commonRowAxis.values)
        newAxes[1] = axesModule.grid(commonColAxis.label, commonColAxis.index, commonColAxis.unit,
                                 style=axesModule.linkGridToken,
                                 values=linkModule.link(link=commonRowAxis.values, relative=True))
        newAxes[0] = axesModule.axis(commonMatrixAxis.label, commonMatrixAxis.index, commonMatrixAxis.unit)
        trigdata = commonMatrix[numpy.tri(commonMatrix.shape[0])==1.0].tolist()
        gridded = griddedModule.gridded2d( axes=newAxes,
            array=arrayModule.full(shape=commonMatrix.shape,data=trigdata,symmetry=arrayModule.symmetryLowerToken) )
        newCov = covarianceMatrix(label=label, type='absolute', matrix=gridded)
        newCov.setAncestor(self.ancestor)
        return newCov

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
        """

        :param indent:
        :param kwargs:
        :return:
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s" domainMin="%s" domainMax="%s" domainUnit="%s">' %
                ( indent, self.moniker, self.label, self.domainMin, self.domainMax, self.domainUnit ) ]
        xmlString.append( indent2 + '<!-- The matrix for this reaction equals the weighted sum of the following matrices: -->' )
        for pointer in self.summands:
            xmlString.append( pointer.toXML( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def getReferredCovariance( self, pointer ):
        """
        FIXME what's this for? It's not used anywhere in FUDGE... should it be removed?

        :param pointer:
        :return:
        """
        if pointer.link is not None: return pointer.link
        if 'covarianceSuite' in pointer.path:   return pointer.link.follow( pointer.path, self.getRootAncestor() )
        #elif'reactionSuite' in pointer.path:    return link.follow( pointer.path, None )
        else :                                  raise ValueError( "Need reference to root node of "+str( pointer.path ) )

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combine all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs1d instance 
        """
        covar = self.toCovarianceMatrix()
        if relative:
            if theData is None:
                theData = self.ancestor.rowData.link.toPointwise_withLinearXYs(lowerEps=1e-8)
            covar = covar.toRelative(rowData=theData)

        return covar.getUncertaintyVector()

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        label = element.get( "label" )
        lower = float( element.get("domainMin") )
        upper = float( element.get("domainMax") )
        summed_ = summedCovariance(label=label, domainMin=lower, domainMax=upper, domainUnit=element.get('domainUnit'))
        for summandElement in element:
            link_ = summand.parseXMLNode( summandElement, xPath, linkData )
            linkData['unresolvedLinks'].append( link_ )
            summed_.summands.append( link_ )
        xPath.pop()
        return summed_

class summand( linkModule.link ):
    """
    Stores one summand in the summed covariance. Consists of a link, an ENDF MF/MT and a coefficient
    """

    moniker = 'summand'

    def __init__( self, link=None, root=None, path=None, label=None, relative=False,
                  ENDF_MFMT=None, coefficient=None ):
        linkModule.link.__init__(self, link=link, root=root, path=path, label=label, relative=relative)
        self.ENDF_MFMT = ENDF_MFMT
        if coefficient is not None:
            coefficient = float(coefficient)
        self.coefficient = coefficient

    def __eq__(self, other):
        for attr in ('link', 'path', 'root', 'ENDF_MFMT', 'coefficient'):
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def toXMLList( self, indent = '', **kwargs ) :

        attributesStr = ""
        for attr in ('ENDF_MFMT', 'coefficient'):
            if getattr(self, attr):
                attributesStr += ' %s="%s"' % (attr, getattr(self, attr))

        return [ '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributesStr, self.build_href( **kwargs ) ) ]
