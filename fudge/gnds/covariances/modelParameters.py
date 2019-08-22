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

from xData import ancestry as ancestryModule
from xData import link as linkModule
from xData import array as arrayModule
from fudge.gnds import suites as suitesModule
from fudge.gnds import abstractClasses as abstractClassesModule
from . import section as sectionModule

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

        self.updateXPath()
        if self.root is not None:    # external link
            href = '%s#%s' % ( self.root, self.path )
        else:                               # link within this file
            href = self.path

        attrString = ' label="%s" href="%s"' % (self.label, href)
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

    def check( self, info ):
        """ check each section """

        from fudge.gnds import warning
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
        from fudge.gnds import warning
        warnings = []

        if self.parameters.nParameters != self.matrix.nrows:
            warnings.append( "need real warning here" ) # FIXME
        matrixWarnings = self.matrix.check( info )
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


class averageParameterCovariance( sectionModule.section ):
    """
    For storing unresolved resonance parameter covariances. Very similar to covariances for cross section, nubar, etc.:
    they require an energy grid + matrix.

    Each average parameter (e.g. elastic width, capture width, etc.) has its own averageParameterCovariance section.
    """

    moniker = 'averageParameterCovariance'
