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
from xData import link as linkModule
from fudge.core.math import matrix as gndMatrix
from pqu import PQU

__metaclass__ = type

class summedCovariance( ancestry ):
    """ 
    Covariance matrix stored as sum/difference of other matrices. 
    """

    moniker = 'sum'

    def __init__(self, label = None, index=None, lowerBound=None, upperBound=None,
            coefficients=None, pointerList=None):
        ancestry.__init__( self )
        self.__label = label
        self.index = index #: an int, only needed within 'mixed' section
        self.lowerBound = lowerBound #: lower bound of row/column direction.  (type: PQU)
        self.upperBound = upperBound #: upper bound of row/column direction.  (type: PQU)
        self.pointerList = pointerList or [] #a list of pointers (type: xData.link.link).  Each pointer has a dictionary with entry "coefficient" to use to weight the covariance referred to

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = value

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def getRowBounds(self,unit=None):
        '''Get the bounds of the row.  If unit is specified, return the bounds in that unit.'''
        if unit is None: return (self.lowerBound,self.upperBound)
        else: return (self.lowerBound.getValueAs(unit),self.upperBound.getValueAs(unit))
    
    def getColumnBounds(self,unit=None):
        '''Get the bounds of the column.  If unit is specified, return the bounds in that unit.'''
        if unit is None: return (self.lowerBound,self.upperBound)
        else: return (self.lowerBound.getValueAs(unit),self.upperBound.getValueAs(unit))

    def toCorrelationMatrix( self ):
        return self.toCovarianceMatrix().toCorrelationMatrix()

    def toCovarianceMatrix( self ): 
        '''
        Sum the parts to construct the covariance matrix.  
        Note, each part must be converted to an absolute covariance before summing.
        '''
        if len( self.pointerList ) == 1: return self.pointerList[0].link.forms['eval'].toCovarianceMatrix()
        from fudge.core.utilities.brb import uniquify
        import numpy, copy
        from .base import covarianceMatrix
        from .mixed import mixedForm
        
        # utility function to get a component as an absolute matrix over the correct row/column bounds
        # need special coding if mixed since only need the part of a mixedForm that overlaps with the current 
        # covariance
        def __get_abs_cov_mtx( ptr ):
            c = ptr.link.forms['eval']
            if isinstance( c, mixedForm ): c=c.shrinkToBounds(self.getRowBounds())
            return c.toCovarianceMatrix()
        
        # set up common data using first element in pointerList
        firstCovMtx = __get_abs_cov_mtx(self.pointerList[0]) #.link.forms['eval'].toCovarianceMatrix().toAbsolute()
        commonRowAxis = copy.copy( firstCovMtx.matrix.axes[0] )
        if firstCovMtx.matrix.axes[1].mirrorOtherAxis: commonColAxis = copy.copy( firstCovMtx.matrix.axes[0] )
        else:                                   commonColAxis = copy.copy( firstCovMtx.matrix.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.matrix.axes[2] )
        commonType = firstCovMtx.type
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for p in self.pointerList[1:]:
            cc = __get_abs_cov_mtx(p) #.link.forms['eval'].toCovarianceMatrix().toAbsolute() # a little recursion to take care of nested covariances
            if cc.type != commonType: raise ValueError( "Incompatable types in "+str(self.__class__)+": "+str(commonType)+' vs. '+str(cc.type) )
            cc.convertAxesToUnits( ( commonRowAxis.unit, commonColAxis.unit, commonMatrixAxis.unit ) )
            commonRowAxis.data = commonRowAxis.data + cc.matrix.axes[0].data
            if cc.matrix.axes[1].mirrorOtherAxis: commonColAxis.data = commonColAxis.data + cc.matrix.axes[0].data
            else:                          commonColAxis.data = commonColAxis.data + cc.matrix.axes[1].data
        commonRowAxis.data.sort()
        commonColAxis.data.sort()
        commonRowAxis.data = uniquify( commonRowAxis.data )
        commonColAxis.data = uniquify( commonColAxis.data )
        
        # now sum up the components
        commonMatrix = self.pointerList[0]['coefficient'] * numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for p in self.pointerList[1:]:
            cc = p.link.forms['eval'].toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += p['coefficient'] * numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        return covarianceMatrix( type=commonType, axes=[ commonRowAxis, commonColAxis, commonMatrixAxis ], matrix=gndMatrix.matrix( commonMatrix.tolist(), form=gndMatrix.symmetricFormToken ) )
        

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s' % ( indent, self.moniker ) ]
        if self.label is not None: xmlString[0] += ' label="%s"' % self.label
        if self.index is not None: xmlString[0] += ' index="%s"' % self.index
        lowerBound = self.lowerBound.toString( keepPeriod = False )
        upperBound = self.upperBound.toString( keepPeriod = False )
        xmlString[0] += ' lowerBound="%s" upperBound="%s">' % ( lowerBound, upperBound )
        xmlString.append( indent2 + '<!-- The matrix for this reaction equals the weighted sum of the following matrices: -->' )
        for pointer in self.pointerList:
            xmlString.append( pointer.toXML( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString
        
    def getReferredCovariance( self, pointer ):
        if pointer.link is not None: return pointer.link
        if 'covarianceSuite' in pointer.path:   return pointer.link.follow( pointer.path, self.getRootAncestor() )
        #elif'reactionSuite' in pointer.path:    return link.follow( pointer.path, None )
        else :                                  raise ValueError( "Need reference to root node of "+str( pointer.path ) )

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combine all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs instance 
        """
        #xys = [ self.getReferredCovariance( p ).getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        xys = [ p.link.forms['eval'].getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        uncert = coefficients[0] * xys[0]
        for i, xy in enumerate( xys[1:] ):
            uncert, xy = uncert.mutualify( 1e-8, 0, 0, xy, 1e-8, 0, 0 )
            uncert += coefficients[ i ] * xy
        return uncert

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        label = element.get( "label" )
        index = element.get("index")
        if index is not None: index=int(index)
        lower = PQU.PQU( element.get("lowerBound") )
        upper = PQU.PQU( element.get("upperBound") )
        summed_ = summedCovariance( label = label, index=index, lowerBound=lower, upperBound=upper )
        for summandElement in element:
            link_ = summand.parseXMLNode( summandElement, xPath, linkData )
            link_['coefficient'] = float( link_['coefficient'] )
            linkData['unresolvedLinks'].append( link_ )
            summed_.pointerList.append( link_ )
        xPath.pop()
        return summed_

class summand( linkModule.link ):
    moniker = 'summand'
