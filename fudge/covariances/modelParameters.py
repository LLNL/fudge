# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import link as linkModule
from xData import table as tableModule
from xData import xDataArray as arrayModule

from fudge import suites as suitesModule
from fudge import abstractClasses as abstractClassesModule
from fudge.resonances import scatteringRadius as scatteringRadiusModule

from . import enums as covarianceEnumsModule
from . import covarianceSection as sectionModule

class ParameterLink( linkModule.Link ):
    """
    Establishes a link between one or more rows of a parameterCovariance and corresponding parameter(s).
    Supports linking to specific parameters inside a table or list.

    For example, if we have a 2x4 table::

      |A B C D|
      |E F G H|

    and wish to give a 4x4 covariance matrix for elements in the 2nd and 4th column of the table, we can create
    a parameterLink pointing to the table, with
    'matrixStartIndex=0',  'nParameters=4',  'parameterStartIndex=1',  'parameterStride=2'.

    The corresponding covariance matrix rows would then correspond to 'B, D, F, H'.
    """

    moniker = 'parameterLink'

    def __init__( self, label, link, nParameters=1, matrixStartIndex=0, parameterStartIndex=0, parameterStride=1, **kwargs ):
        """
        :param label:               unique label for this parameter link
        :param link:                reference to parameter, table or list
        :param nParameters:         number of parameters linked to
        :param matrixStartIndex:    row of the covariance matrix corresponding to the start of this batch of parameters
        :param parameterStartIndex: for skipping one or more parameters at the start of a table or list
        :param parameterStride:     for stepping over parameters. Stride must be an integer
        :return:
        """
        linkModule.Link.__init__( self, link=link, **kwargs )
        self.label = label
        self.nParameters = int(nParameters)
        self.matrixStartIndex = int(matrixStartIndex)
        self.parameterStartIndex = int(parameterStartIndex)
        self.parameterStride = int(parameterStride)

    def toXML_strList( self, indent = '', **kwargs ) :

        attrString = ' label="%s" href="%s"' % (self.label, self.build_href(**kwargs))
        if self.nParameters != 1: attrString += ' nParameters="%d"' % self.nParameters
        if self.matrixStartIndex != 0: attrString += ' matrixStartIndex="%d"' % self.matrixStartIndex
        if self.parameterStartIndex != 0: attrString += ' parameterStartIndex="%d"' % self.parameterStartIndex
        if self.parameterStride != 1: attrString += ' parameterStride="%d"' % self.parameterStride
        return [ '%s<%s%s/>' % ( indent, self.moniker, attrString ) ]


class Parameters( suitesModule.Suite ):

    moniker = 'parameters'

    def __init__( self ):
        suitesModule.Suite.__init__( self, allowedClasses = (ParameterLink,) )

    @property
    def nParameters(self):
        return sum( [link.nParameters for link in self] )


class ParameterCovariance( suitesModule.Suite ):
    """
    For storing unresolved resonance parameter covariances. Very similar to covariances for cross section, nubar, etc.:
    they require an energy grid + matrix.

    Each average parameter (e.g. elastic width, capture width, etc.) has its own averageParameterCovariance section.
    """

    moniker = 'parameterCovariance'
    keyName = 'label'

    def __init__(self, label, rowData=None, columnData=None):

        suitesModule.Suite.__init__( self, [ParameterCovarianceMatrix] )
        self.label = label
        self.rowData = rowData
        self.columnData = columnData

    @property
    def crossTerm( self ):
        return self.columnData is not None and self.columnData.link != self.rowData.link

    @property
    def evaluated(self):
        """
        Helper method to grab evaluated style
        FIXME method should be inherited, but abstractClasses.Component defines methods that don't make sense for covariances
        """

        if not hasattr(self.rootAncestor, 'styles'):
            return self[0]  # For covarianceSection that is not part of a full covarianceSuite.

        evaluated = self.rootAncestor.styles.getEvaluatedStyle()
        try:
            return self[evaluated.label]
        except IndexError:
            return self[0]

    def check( self, info ):
        """ check each section """

        from fudge import warning
        warnings = []
        for form in self:
            formWarnings = form.check(info)
            if formWarnings:
                warnings.append(warning.Context("Form '%s':" % form.label, formWarnings))

        return warnings

    def toXML_strList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = [indent + '<%s label="%s"' % (self.moniker, self.label)]
        if self.crossTerm: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData', 'columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append(getattr(self, dataPointer).toXML(indent2, **kwargs))
        for form in self:
            xmlString += form.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Translate <section> element from xml."""

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))

        linkData['typeConversion'] = {'domainMin': float, 'domainMax': float}
        rowData_ = sectionModule.RowData.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        columnData_ = None
        if element[1].tag == "columnData":
            columnData_ = sectionModule.ColumnData.parseNodeUsingClass(element[1], xPath, linkData, **kwargs)
        del linkData['typeConversion']
        section_ = cls(element.get('label'), rowData_, columnData_)
        start = 2 if (columnData_ is not None) else 1
        for form in element[start:]:
            formClass = {
                ParameterCovarianceMatrix.moniker: ParameterCovarianceMatrix,
            }.get(form.tag)
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.add(formClass.parseNodeUsingClass(form, xPath, linkData, **kwargs))

        xPath.pop()

        return section_


class ParameterCovarianceMatrix( abstractClassesModule.Form ):
    """
    Store covariances (or correlations, depending on 'type') between model parameters
    """

    moniker = 'parameterCovarianceMatrix'
    keyName = 'label'

    def __init__(self, label, matrix, parameters_=None, type=covarianceEnumsModule.Type.relative):
        """

        :param label:
        :param type:        'relative' or 'absolute'
        :param parameters_: list of parameterLinks that relate rows/columns to parameters
        :param matrix:      xData.array instance containing the covariance or correlation
        :return:
        """

        abstractClassesModule.Form.__init__( self )
        self.label = label
        self.type = covarianceEnumsModule.Type.checkEnumOrString(type)
        if parameters_ is not None:
            self.parameters = parameters_
        else:
            self.parameters = Parameters()
        self.parameters.setAncestor( self )
        matrix.setAncestor( self )
        self.matrix = matrix

    def check( self, info ):
        import numpy
        from fudge import warning
        warnings = []

        params = self.parameters.nParameters
        if self.matrix.shape != (params,params):
            warnings.append( warning.ParameterCovarianceMismatch(params, self.matrix.shape) )

        else:
            matrix = self.matrix.constructArray()
            matrixWarnings = []

            params = []
            for parameterLink in self.parameters:
                if isinstance(parameterLink.link, tableModule.Table):
                    for row in parameterLink.link.data:
                        params += row
                elif isinstance(parameterLink.link, scatteringRadiusModule.ScatteringRadius):
                    params.append(parameterLink.link.form.value)
                else:
                    raise NotImplementedError("Model parameter covariance linking to %s" % type(parameterLink.link))
            params = numpy.array(params)

            if self.type == covarianceEnumsModule.Type.absolute:
                # params may contain zeros, e.g. for L/J columns in BreitWigner table. Replace with 1
                params2 = params[:]
                params2[ params2==0 ] = 1
                matrix = matrix / params2 / params2[:,numpy.newaxis]

            uncertainty = numpy.sqrt(matrix.diagonal())
            for idx, unc in enumerate(uncertainty):
                if unc < info['minRelUnc']:
                    # FIXME should probably ignore fake resonances (e.g. with negative energies)
                    warnings.append( warning.VarianceTooSmall( idx, unc**2, self ) )
                elif unc > info['maxRelUnc']:
                    warnings.append( warning.VarianceTooLarge( idx, unc**2, self ) )

            vals, vecs = numpy.linalg.eigh( matrix )
            if min(vals) < info['negativeEigenTolerance']:
                matrixWarnings.append( warning.NegativeEigenvalues( len(vals[vals<0]), min(vals), self ))
            minpos, maxpos = min(vals[vals>=0]),max(vals[vals>=0])
            # Check that the condition number of the matrix is reasonable
            if minpos/maxpos < info['eigenvalueRatioTolerance'] and matrix.size != 1:
                matrixWarnings.append( warning.BadEigenvalueRatio( minpos/maxpos, self ) )

            if matrixWarnings:
                warnings.append( warning.Context("Model parameter covariances", matrixWarnings ) )

        return warnings

    def convertUnits( self, unitMap ):

        #raise NotImplementedError()
        pass

    def fix( self, **kw ):
        '''assemble some useful info, to be handed down to children's fix() functions'''
        info = {}
        info['rowENDF_MFMT'] = None
        info['columnENDF_MFMT'] = None
        info.update( kw )
        return self.matrix.fix( **info )

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

        #print(self.parameters)
        #for p in self.parameters: print(p)

        _array = self.matrix.constructArray()

        def covariance_to_correlation(__matrix):
            """Convert a covariance matrix to a correlation matrix"""
            diag = np.sqrt(__matrix.diagonal())
            corr = __matrix / diag / diag[:, np.newaxis]
            # now fix diagonal + remove any NaN (from div/0):
            for i in range(len(corr)):
                corr[i, i] = 1.0  # must be exactly 1
            corr[np.isnan(corr)] = 0
            return corr

        _array = covariance_to_correlation(_array)

        x, y = range(_array.shape[0]+1), range(_array.shape[1]+1)
        X, Y = np.meshgrid( x, y )
        XY = np.hstack((X.ravel()[:,np.newaxis], Y.ravel()[:,np.newaxis]))

        ax = plt.subplot(1,1,1)
        if title is None:
            title = str( self.toXLink() )
        plt.suptitle(title)

        qc = QuadMesh(
            meshWidth=len(x)-1,
            meshHeight=len(y)-1,
            coordinates=XY,
            antialiased=True,
            shading='flat',
            transOffset=ax.transData)

        qc.set_array(_array)
        ax.add_collection(qc, autolim=True)

        if xlim is None:
            ax.set_xlim( x[0], x[-1] )
        else:
            ax.set_xlim( xlim[0], xlim[1] )
        if ylim is None:
            ax.set_ylim( y[0], y[-1] )
        else:
            ax.set_ylim( ylim[0], ylim[1] )
        if xlog:
            ax.set_xscale( 'log' )
        if ylog:
            ax.set_yscale( 'log' )

        xlabel = 'row index'
        ylabel = 'col index'

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        cbar = plt.colorbar(qc)
        if scalelabel is not None:
            cbar.set_label(scalelabel)
        else:
            cbar.set_label(str(self.type)+' correlation' )
        plt.show()

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmllist = [ '%s<%s label="%s" type="%s">' % (indent, self.moniker, self.label, self.type) ]
        xmllist += self.parameters.toXML_strList( indent2, **kwargs )
        xmllist += self.matrix.toXML_strList( indent2, **kwargs )
        xmllist[-1] += ('</%s>' % self.moniker)
        return xmllist

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        matrix_ = arrayModule.ArrayBase.parseNodeUsingClass(element.find( arrayModule.Full.moniker ), xPath, linkData, **kwargs)
        PC = cls(element.get('label'), matrix_, type=element.get('type'))
        PC.parameters.parseNode(element.find( Parameters.moniker ), xPath, linkData, **kwargs)

        xPath.pop()

        return PC


class AverageParameterCovariance(sectionModule.CovarianceSection):
    """
    For storing unresolved resonance parameter covariances. Very similar to covariances for cross section, nubar, etc.:
    they require an energy grid + matrix.

    Each average parameter (e.g. elastic width, capture width, etc.) has its own averageParameterCovariance section.
    """

    moniker = 'averageParameterCovariance'
