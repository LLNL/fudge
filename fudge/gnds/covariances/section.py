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

"""A covarianceSuite is organized into sections, where each section contains either
    - a covariance matrix for a single reaction quantity (cross section, multiplicity, etc), or
    - a covariance matrix between two different quantities (off-diagonal block)
    """
from xData.ancestry import ancestry
from xData import link as linkModule
from fudge.gnds import suites as suitesModule
from . import covarianceMatrix
from .mixed import mixedForm
from .summed import summedCovariance
from pqu import PQU as PQUModule

__metaclass__ = type

class section( suitesModule.suite ):
    """
    A covarianceSuite contains sections, where each section represents either a self-covariance for one quantity,
    or a cross-covariance between two quantities

    More generally, the covarianceSuite can be thought of as a single covariance matrix with all covariance data
    for a target/projectile. It is broken into sections, where each section holds a chunk of the full matrix.

    Within each section, covariance data can take multiple forms: :py:class:`covarianceMatrix` is the most common,
    but 'summed', 'mixed' are also possible. 

    Valid values in the :py:attr:`forms` dictionary are:
        * mixedForm
        * summedCovariance
        * covarianceMatrix
    """

    moniker = 'section'

    def __init__(self, label, rowData=None, columnData=None):
        """ each section needs a unique id, pointers to the central values (row and column),
        and one or more forms """

        suitesModule.suite.__init__( self, [covarianceMatrix, mixedForm, summedCovariance] )
        self.label = label #: a str label that gets used on plots, etc.
        self.rowData = rowData #: xData.link.link pointing to the corresponding data for the covariance row
        self.columnData = columnData #: xData.link.link pointing to the corresponding data for the covariance column

    @property
    def crossTerm(self):
        return self.columnData is not None and self.columnData.link != self.rowData.link

    def check( self, info ):
        """ check each section """

        from fudge.gnds import warning
        warnings = []
        for form in self:
            formWarnings = form.check( info )
            if formWarnings:
                warnings.append( warning.context( "Form '%s':" % form.label, formWarnings ) )

        return warnings
    
    def fix( self, **kw ): 
        """assemble some useful info, to be handed down to children's check() functions"""
        info = {}
        warnings = []
        if self.rowData is None:    info['rowENDF_MFMT'] = None
        else:                       info['rowENDF_MFMT'] = self.rowData['ENDF_MFMT']
        if self.columnData is None: info['columnENDF_MFMT'] = None
        else:                       info['columnENDF_MFMT'] = self.columnData['ENDF_MFMT']
        info.update( kw )
        for form in self: warnings += form.fix( **info )
        return warnings

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [indent+'<%s label="%s"' % (self.moniker, self.label)]
        if self.crossTerm: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData','columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append( getattr(self, dataPointer).toXML( indent2, **kwargs ) )
        for form in self:
            xmlString += form.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        """Translate <section> element from xml."""

        xPath.append( '%s[@label="%s"]' % (element.tag, element.get('label') ) )
        linkData['typeConversion'] = {'domainMin':float, 'domainMax':float}
        rowData_ = rowData.parseXMLNode( element[0], xPath, linkData )
        columnData_ = None
        if element[1].tag=="columnData":
            columnData_ = columnData.parseXMLNode( element[1], xPath, linkData )
        del linkData['typeConversion']
        section_ = cls( element.get('label'), rowData_, columnData_ )
        start = 2 if (columnData_ is not None) else 1
        for form in element[start:]:
            formClass = {
                    covarianceMatrix.moniker: covarianceMatrix,
                    mixedForm.moniker: mixedForm,
                    summedCovariance.moniker: summedCovariance,
                    }.get( form.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.add( formClass.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return section_

class rowData( linkModule.link ):

    moniker = 'rowData'

class columnData( linkModule.link ):

    moniker = 'columnData'
