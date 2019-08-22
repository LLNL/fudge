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

from xData import ancestry as ancestryModule
from xData import link as linkModule
from fudge.core.math import matrix as gndMatrix

__metaclass__ = type

class inputParameter( ancestryModule.ancestry ):
    """
    For use within a modelParameterCovariance: the rows/columns of that matrix each correspond to one input
    model parameter. Parameters must each be inputParameter or multipleInputParameter instances.
    """

    moniker = 'inputParameter'

    def __init__(self, name, path, unit=None):
        ancestryModule.ancestry.__init__( self )
        self.name = name #: just a Python str
        self.path = path #: a fudge.gnd.link 
        self.unit = unit #: an acceptable PQU unit

    def toXMLList( self, indent = '', **kwargs ) :

        if self.unit: unit = ' unit="%s"' % self.unit
        else: unit = ''
        return [ '%s<parameter name="%s"%s xlink:href="%s"/>' % ( indent, self.name, unit, self.path.toXLink() ) ]

class modelParameterCovariance( ancestryModule.ancestry ):
    """ 
    Express covariance between input parameters for a model 
    """

    moniker = 'modelParameterCovariance'

    def __init__(self, label=None, inputParameters=None, matrix=None, type=None, **kwargs):

        ancestryModule.ancestry.__init__( self )
        self.label = label #: a str
        self.inputParameters = inputParameters #: list of inputParameter instances to relate rows/columns to parameters 
        for param in self.inputParameters: param.setAncestor( self )
        self.matrix = matrix #: the actual fudge.core.math.matrix instance with the covariance
        #self.matrix.setAncestor( self )
        self.type = type  #: dunno, got to ask Caleb or Bret
        self.tag = 'modelParameterCovariance' #: usually 'modelParameterCovariance'
        self.attributes = kwargs #: a Python dict
        
    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        matrixWarnings = self.matrix.check( info )
        if matrixWarnings:
            warnings.append( warning.context("Model parameter covariances", matrixWarnings ) )
        return warnings
    
    def fix( self, **kw ): 
        '''assemble some useful info, to be handed down to children's fix() functions'''
        info = {}
        info['rowENDF_MFMT'] = None
        info['columnENDF_MFMT'] = None
        info.update( kw )
        return self.matrix.fix( **info )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrStr = ''.join( [' %s="%s"' % (key, self.attributes[key]) for key in self.attributes
            if bool(self.attributes[key]) ] )
        xmllist = [indent+( '<%s label="%s" type="%s" %s>' % (self.tag, self.label, self.type, attrStr) )]
        xmllist.extend( [indent+'  <inputParameters>',
            indent2+'<!-- Each row of this matrix corresponds to a model parameter. Parameters may be listed singly,',
            indent2+'  as for scattering radii, or in the case of resonance parameters they may be given all together.',
            indent2+'  In that case, rows of the matrix correspond to a loop over parameters for each resonance,',
            indent2+'  with resonances sorted by energy. -->'] )
        for inputParam in self.inputParameters: xmllist += inputParam.toXMLList( indent2, **kwargs )
        xmllist[-1] += '</inputParameters>'
        xmllist += self.matrix.toXMLList( indent2, **kwargs )
        xmllist[-1] += ('</%s>' % self.tag)
        return xmllist

class loopOverResonanceParameters( linkModule.link ):
    """ 
    For resonance region covariances, we need a compact way to express many model inputs.
    Simplest is to specify a loop over the resonances 
    """
    moniker = 'loopOverResonanceParameters'

    def __init__(self, link=None, root=None, path=None, nResonances=0, parametersPerResonance=''):

        linkModule.link.__init__(self, link=link, root=root, path=path)
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

class resonanceParameterCovariance( modelParameterCovariance ):
    """
    In the resonance region, covariances are given between resonance parameters (energy and widths).
    Generally, the dimension of the matrix is 3*(number of resonances) for light targets, and 4*(nres)
    for heavy targets (where the fission width must be given).
    
    We also allow including the scattering radius in the covariance, although ENDF files currently only
    have room to list the uncertainty (variance) on the scattering radius. 
    """

    def __init__(self, label=None, inputParameters=None, matrix=None, type=None, **kwargs):
        modelParameterCovariance.__init__(self, label, inputParameters, matrix, type)
        self.tag = 'resonanceParameterCovariance' #: usually set to 'resonanceParameterCovariance'
        self.attributes = kwargs #: a Python dict

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
        RPCs = resonanceParameterCovariance( inputParameters=params, matrix=Matrix, **dict(element.items()) )
        xPath.pop()
        return RPCs
