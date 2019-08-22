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

"""
Some evaluations use short-range self-scaling variance components  (ENDF LB=8 or 9)
to represent part or all of the uncertainty.
When processing these sections to generate a multi-group covariance library,
the size of the variance depends either directly or inversely on the size of the processed group.
These sections only produce a diagonal, no cross-correlations.
"""

from . import tokens

from xData import ancestry as ancestryModule
from xData import array as arrayModule
from xData import gridded as griddedModule

from pqu import PQU

__metaclass__ = type

# define types of processed width dependences:
inverseToken = 'inverse'
directToken = 'direct'

class shortRangeSelfScalingVariance(ancestryModule.ancestry):

    moniker = 'shortRangeSelfScalingVariance'

    def __init__( self, label, type=tokens.absoluteToken, dependenceOnProcessedGroupWidth=inverseToken, matrix=None ):

        ancestryModule.ancestry.__init__(self)
        self.__label = label
        self.dependenceOnProcessedGroupWidth = dependenceOnProcessedGroupWidth
        self.type = type
        self.matrix = matrix

    @property
    def label(self): return self.__label

    @label.setter
    def label(self, value):
        if value is not None:
            if not isinstance(value, str): raise TypeError( 'label must be a string' )
        self.__label = value

    @property
    def dependenceOnProcessedGroupWidth(self):
        return self.__dependenceOnProcessedGroupWidth

    @dependenceOnProcessedGroupWidth.setter
    def dependenceOnProcessedGroupWidth(self, value):
        if value not in (inverseToken, directToken):
            raise TypeError("Unrecognized dependenceOnProcessedGroupWidth flag '%s'" % value)
        self.__dependenceOnProcessedGroupWidth = value

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, type):
        if type not in (tokens.absoluteToken, tokens.relativeToken):
            raise TypeError("Unrecognized type '%s'" % type)
        self.__type = type

    @property
    def matrix(self):
        return self.__matrix

    @matrix.setter
    def matrix(self, matrix):
        if not isinstance(matrix, griddedModule.gridded2d):
            raise TypeError("Expected gridded2d instance, got '%s'" % type(matrix))
        self.__matrix = matrix

    @property
    def gridded2d(self): return self.matrix     # convenience method

    def check( self, info ):
        """
        Check that matrix is diagonal and all elements >= 0
        """
        from fudge.gnds import warning
        warnings = []

        if not isinstance( self.matrix.array, arrayModule.diagonal ):
            warnings.append( warning.invalidShortRangeVarianceData( type(self.matrix.array) ) )
        eigenvals = self.matrix.array.constructArray().diagonal()
        if any( eigenvals < 0 ):
            warnings.append( warning.negativeEigenvalues(negativeCount=(eigenvals < 0).sum(),
                    worstCase=eigenvals.min(), obj=self) )
        return warnings

    def getRowBounds(self,unit=None):
        """Get the bounds of the row.  If unit is specified, return the bounds in that unit."""
        factor = 1
        if unit:
            factor = PQU.PQU(1, self.matrix.axes[-1].unit).getValueAs(unit)
        return( self.matrix.axes[-1].values[0] * factor, self.matrix.axes[-1].values[-1] * factor )

    def getColumnBounds(self,unit=None):
        """Get the bounds of the column.  If unit is specified, return the bounds in that unit."""
        if self.matrix.axes[-2].style=='link': return self.getRowBounds(unit)
        factor = 1
        if unit:
            factor = PQU.PQU(1, self.matrix.axes[-1].unit).getValueAs(unit)
        return( self.matrix.axes[-1].values[0] * factor, self.matrix.axes[-1].values[-1] * factor )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ''
        if self.label is not None: attrs = ' label="%s"' % self.label
        xmlString = [ '%s<%s%s type="%s" dependenceOnProcessedGroupWidth="%s">' %
                      ( indent, self.moniker, attrs, self.type, self.dependenceOnProcessedGroupWidth ) ]
        xmlString += self.matrix.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        label = element.get('label')
        if label is not None:
            xPath.append( '%s[@label="%s"]' % (element.tag, label))
        else:
            xPath.append( element.tag )
        matrix_ = griddedModule.gridded2d.parseXMLNode( element.find(griddedModule.gridded2d.moniker), xPath, linkData )
        srssv = cls( label = label, type=element.get('type'),
            dependenceOnProcessedGroupWidth=element.get('dependenceOnProcessedGroupWidth'), matrix=matrix_ )
        xPath.pop()
        return srssv
