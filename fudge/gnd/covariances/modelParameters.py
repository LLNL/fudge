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
from fudge.gnd import suites as suitesModule
from fudge.core.math import matrix as gndMatrix

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
        attrString = ' label="%s" xlink:href="%s"' % (self.label, self.path)
        if self.nParameters != 1: attrString += ' nParameters="%d"' % self.nParameters
        if self.matrixStartIndex != 0: attrString += ' matrixStartIndex="%d"' % self.matrixStartIndex
        if self.parameterStartIndex != 0: attrString += ' parameterStartIndex="%d"' % self.parameterStartIndex
        if self.parameterStride != 1: attrString += ' parameterStride="%d"' % self.parameterStride
        return [ '%s<%s%s/>' % ( indent, self.moniker, attrString ) ]


class parameters( suitesModule.suite ):

    moniker = 'parameters'

    def __init__( self ):
        suitesModule.suite.__init__( self, allowedClasses = (parameterLink, loopOverResonanceParameters) )

    @property
    def nParameters(self):
        return sum( [link.nParameters for link in self] )


class parameterCovariance( ancestryModule.ancestry ):   # FIXME should inherit from 'form' in abstractBaseClasses
    """ 
    Store covariances (or correlations, depending on 'type') between model parameters
    """

    moniker = 'parameterCovariance'

    def __init__( self, label, matrix, parameters_=None, type='relativeCovariance', ENDFconversionFlags='' ):
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
        self.ENDFconversionFlags = ENDFconversionFlags

    def check( self, info ):
        from fudge.gnd import warning
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

        attrs = ''
        if self.ENDFconversionFlags: attrs = ' ENDFconversionFlags="%s"' % self.ENDFconversionFlags
        xmllist = [ '%s<%s label="%s" type="%s"%s>' %
                    (indent, self.moniker, self.label, self.type, attrs) ]
        xmllist += self.parameters.toXMLList( indent2, **kwargs )
        xmllist += self.matrix.toXMLList( indent2, **kwargs )
        xmllist[-1] += ('</%s>' % self.moniker)
        return xmllist

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        matrix_ = arrayModule.full.parseXMLNode( element.find( arrayModule.full.moniker ), xPath, linkData )
        PC = parameterCovariance(element.get('label'), matrix_, type=element.get('type'),
            ENDFconversionFlags=element.get('ENDFconversionFlags',''))
        PC.parameters.parseXMLNode( element.find( parameters.moniker ), xPath, linkData )
        xPath.pop()
        return PC


class loopOverResonanceParameters( linkModule.link ):       # FIXME get rid of this, replace with parameterLink
    """ 
    For resonance region covariances, we need a compact way to express many model inputs.
    Simplest is to specify a loop over the resonances 
    """
    moniker = 'loopOverResonanceParameters'

    def __init__(self, link=None, root=None, path=None, relative=False, nResonances=0, parametersPerResonance=''):

        linkModule.link.__init__(self, link=link, root=root, path=path, relative=relative)
        self.attributes = {
            'nResonances': nResonances, #: an int, the number of resonances
            'parametersPerResonance': parametersPerResonance, #: dunno, got to ask Caleb or Bret
        }

    @property
    def nResonances(self): return self.attributes['nResonances']

    @property
    def parametersPerResonance(self): return self.attributes['parametersPerResonance']

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        LORPs = super(loopOverResonanceParameters,cls).parseXMLNode( element, xPath, linkData )
        LORPs.attributes['nResonances'] = int( LORPs.attributes['nResonances'] )

        return LORPs


class resonanceParameterCovariance( parameterCovariance ):      # FIXME get rid of this, replace with parameterCovariance
    """
    In the resonance region, covariances are given between resonance parameters (energy and widths).
    Generally, the dimension of the matrix is 3*(number of resonances) for light targets, and 4*(nres)
    for heavy targets (where the fission width must be given).
    
    We also allow including the scattering radius in the covariance, although ENDF files currently only
    have room to list the uncertainty (variance) on the scattering radius. 
    """

    def __init__(self, label=None, matrix=None, inputParameters=None, type=None, ENDFconversionFlags=''):
        parameterCovariance.__init__(self, label, matrix, inputParameters, type)
        self.tag = 'resonanceParameterCovariance' #: usually set to 'resonanceParameterCovariance'
        self.ENDFconversionFlags = ENDFconversionFlags

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """Translate <resonanceParameterCovariance> element from xml."""

        xPath.append( element.tag )
        params = []
        for param in element[0]:
            if param.tag=='loopOverResonanceParameters':
                param=loopOverResonanceParameters.parseXMLNode(param, xPath, linkData)
            elif param.tag=='parameter': pass
            params.append(param)
        Matrix = gndMatrix.parseXMLNode( element[1], xPath, linkData )
        RPCs = resonanceParameterCovariance( matrix=Matrix, inputParameters=params, **dict(element.items()) )
        xPath.pop()
        return RPCs
