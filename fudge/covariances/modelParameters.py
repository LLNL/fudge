# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import link as linkModule
from xData import xDataArray as arrayModule
from fudge import suites as suitesModule
from fudge import abstractClasses as abstractClassesModule
from . import covarianceSection as sectionModule, tokens as tokensModule

__metaclass__ = type

class parameterLink( linkModule.link ):
    """
    Establishes a link between one or more rows of a parameterCovariance and corresponding parameter(s).
    Supports linking to specific parameters inside a table or list.

    For example, if we have a 2x4 table
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
        linkModule.link.__init__( self, link=link, **kwargs )
        self.label = label
        self.nParameters = int(nParameters)
        self.matrixStartIndex = int(matrixStartIndex)
        self.parameterStartIndex = int(parameterStartIndex)
        self.parameterStride = int(parameterStride)

    def toXMLList( self, indent = '', **kwargs ) :

        attrString = ' label="%s" href="%s"' % (self.label, self.build_href(**kwargs))
        if self.nParameters != 1: attrString += ' nParameters="%d"' % self.nParameters
        if self.matrixStartIndex != 0: attrString += ' matrixStartIndex="%d"' % self.matrixStartIndex
        if self.parameterStartIndex != 0: attrString += ' parameterStartIndex="%d"' % self.parameterStartIndex
        if self.parameterStride != 1: attrString += ' parameterStride="%d"' % self.parameterStride
        return [ '%s<%s%s/>' % ( indent, self.moniker, attrString ) ]


class parameters( suitesModule.suite ):

    moniker = 'parameters'

    def __init__( self ):
        suitesModule.suite.__init__( self, allowedClasses = (parameterLink,) )

    @property
    def nParameters(self):
        return sum( [link.nParameters for link in self] )


class parameterCovariance( suitesModule.suite ):
    """
    For storing unresolved resonance parameter covariances. Very similar to covariances for cross section, nubar, etc.:
    they require an energy grid + matrix.

    Each average parameter (e.g. elastic width, capture width, etc.) has its own averageParameterCovariance section.
    """

    moniker = 'parameterCovariance'

    def __init__(self, label, rowData=None, columnData=None):

        suitesModule.suite.__init__( self, [parameterCovarianceMatrix] )
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
        FIXME method should be inherited, but abstractClasses.component defines methods that don't make sense for covariances
        """

        if not hasattr(self.getRootAncestor(), 'styles'):
            return self[0]  # For covarianceSection that is not part of a full covarianceSuite.

        evaluated = self.getRootAncestor().styles.getEvaluatedStyle()
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
                warnings.append(warning.context("Form '%s':" % form.label, formWarnings))

        return warnings

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = [indent + '<%s label="%s"' % (self.moniker, self.label)]
        if self.crossTerm: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData', 'columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append(getattr(self, dataPointer).toXML(indent2, **kwargs))
        for form in self:
            xmlString += form.toXMLList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        """Translate <section> element from xml."""

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))
        linkData['typeConversion'] = {'domainMin': float, 'domainMax': float}
        rowData_ = sectionModule.rowData.parseXMLNode(element[0], xPath, linkData)
        columnData_ = None
        if element[1].tag == "columnData":
            columnData_ = sectionModule.columnData.parseXMLNode(element[1], xPath, linkData)
        del linkData['typeConversion']
        section_ = cls(element.get('label'), rowData_, columnData_)
        start = 2 if (columnData_ is not None) else 1
        for form in element[start:]:
            formClass = {
                parameterCovarianceMatrix.moniker: parameterCovarianceMatrix,
            }.get(form.tag)
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.add(formClass.parseXMLNode(form, xPath, linkData))
        xPath.pop()
        return section_


class parameterCovarianceMatrix( abstractClassesModule.form ):
    """ 
    Store covariances (or correlations, depending on 'type') between model parameters
    """

    moniker = 'parameterCovarianceMatrix'

    def __init__( self, label, matrix, parameters_=None, type='relativeCovariance' ):
        """

        :param label:
        :param type:        'relativeCovariance', 'absoluteCovariance' or 'correlation'
        :param parameters_: list of parameterLinks that relate rows/columns to parameters
        :param matrix:      xData.array instance containing the covariance or correlation
        :return:
        """

        ancestryModule.ancestry.__init__( self )
        self.label = label
        self.type = type
        if parameters_ is not None:
            self.parameters = parameters_
        else:
            self.parameters = parameters()
        self.parameters.setAncestor( self )
        matrix.setAncestor( self )
        self.matrix = matrix

    def check( self, info ):
        import numpy
        from fudge import warning
        warnings = []

        params = self.parameters.nParameters
        if self.matrix.shape != (params,params):
            warnings.append( warning.parameterCovarianceMismatch(params, self.matrix.shape) )

        else:
            matrix = self.matrix.constructArray()
            matrixWarnings = []

            params = []
            for parameterLink in self.parameters:
                for row in parameterLink.link.data:
                    params += row
            params = numpy.array(params)

            if self.type == tokensModule.absoluteToken:
                matrix = matrix / params / params[:,numpy.newaxis]

            uncertainty = numpy.sqrt(matrix.diagonal())
            for idx, unc in enumerate(uncertainty):
                if unc < info['minRelUnc']:
                    # FIXME should probably ignore fake resonances (e.g. with negative energies)
                    warnings.append( warning.varianceTooSmall( idx, unc**2, self ) )
                elif unc > info['maxRelUnc']:
                    warnings.append( warning.varianceTooLarge( idx, unc**2, self ) )

            vals, vecs = numpy.linalg.eigh( matrix )
            if min(vals) < info['negativeEigenTolerance']:
                matrixWarnings.append( warning.negativeEigenvalues( len(vals[vals<0]), min(vals), self ))
            minpos, maxpos = min(vals[vals>=0]),max(vals[vals>=0])
            # Check that the condition number of the matrix is reasonable
            if minpos/maxpos < info['eigenvalueRatioTolerance'] and matrix.size != 1:
                matrixWarnings.append( warning.badEigenvalueRatio( minpos/maxpos, self ) )

            if matrixWarnings:
                warnings.append( warning.context("Model parameter covariances", matrixWarnings ) )

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

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmllist = [ '%s<%s label="%s" type="%s">' % (indent, self.moniker, self.label, self.type) ]
        xmllist += self.parameters.toXMLList( indent2, **kwargs )
        xmllist += self.matrix.toXMLList( indent2, **kwargs )
        xmllist[-1] += ('</%s>' % self.moniker)
        return xmllist

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        matrix_ = arrayModule.full.parseXMLNode( element.find( arrayModule.full.moniker ), xPath, linkData )
        PC = cls(element.get('label'), matrix_, type=element.get('type'))
        PC.parameters.parseXMLNode( element.find( parameters.moniker ), xPath, linkData )
        xPath.pop()
        return PC


class averageParameterCovariance(sectionModule.covarianceSection):
    """
    For storing unresolved resonance parameter covariances. Very similar to covariances for cross section, nubar, etc.:
    they require an energy grid + matrix.

    Each average parameter (e.g. elastic width, capture width, etc.) has its own averageParameterCovariance section.
    """

    moniker = 'averageParameterCovariance'
