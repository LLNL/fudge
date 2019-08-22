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

    def __init__(self, label = None, lowerBound=None, upperBound=None,
            coefficients=None, pointerList=None):
        ancestry.__init__( self )
        self.__label = label
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
        if len( self.pointerList ) == 1: return self.pointerList[0].link['eval'].toCovarianceMatrix()
        import numpy, copy
        from .mixed import mixedForm
        from xData import values as valuesModule
        from xData import axes as axesModule
        from xData import array as arrayModule
        from xData import gridded as griddedModule

        # utility function to get a component as an absolute matrix over the correct row/column bounds
        # need special coding if mixed since only need the part of a mixedForm that overlaps with the current 
        # covariance
        def __get_abs_cov_mtx( ptr ):
            c = ptr.link['eval']
            if isinstance( c, mixedForm ): c=c.shrinkToBounds(self.getRowBounds())
            return c.toCovarianceMatrix()
        
        # set up common data using first element in pointerList
        firstCovMtx = __get_abs_cov_mtx(self.pointerList[0]) #.link['eval'].toCovarianceMatrix().toAbsolute()
        commonRowAxis = copy.copy( firstCovMtx.matrix.axes[2] )
        if firstCovMtx.matrix.axes[1].style=='link': commonColAxis = copy.copy( firstCovMtx.matrix.axes[2] )
        else:                                        commonColAxis = copy.copy( firstCovMtx.matrix.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.matrix.axes[0] )
        commonType = firstCovMtx.type
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        
        # We're going to have to merge grids, so we'll need this function to do the dirty work
        def add_values(v1,v2):
            v=set()
            v.update(v1.values)
            v.update(v2.values)
            return valuesModule.values(sorted(v))

        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for p in self.pointerList[1:]:
            cc = __get_abs_cov_mtx(p) #.link['eval'].toCovarianceMatrix().toAbsolute() # a little recursion to take care of nested covariances
            if cc.type != commonType: raise ValueError( "Incompatable types in "+str(self.__class__)+": "+str(commonType)+' vs. '+str(cc.type) )
            cc.matrix.axes[0].unit = commonMatrixAxis.unit
            cc.matrix.axes[1].convertToUnit(commonColAxis.unit)
            cc.matrix.axes[2].convertToUnit(commonRowAxis.unit)
            commonRowAxis._values = add_values(commonRowAxis.values, cc.matrix.axes[2].values)
            if cc.matrix.axes[1].style=='link': commonColAxis.__values = add_values(commonColAxis.values, cc.matrix.axes[2].values)
            else:                               commonColAxis.__values = add_values(commonColAxis.values, cc.matrix.axes[1].values)

        # now sum up the components
        commonMatrix = self.pointerList[0]['coefficient'] * firstCovMtx.group( ( commonRowAxis.values, commonColAxis.values ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray()
        for p in self.pointerList[1:]:
            cc = p.link['eval'].toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += p['coefficient'] * cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.array.constructArray()
        
        # now create the instance of the resulting covarianceMatrix
        newAxes=axesModule.axes(
                labelsUnits={0 : (commonMatrixAxis.label, commonMatrixAxis.unit),
                             1 : (commonColAxis.label, commonColAxis.unit),
                             2 : (commonRowAxis.label, commonRowAxis.unit)} )
        newAxes[0]=commonMatrixAxis
        newAxes[1]=commonColAxis
        newAxes[2]=commonRowAxis
        trigdata = commonMatrix[numpy.tri(commonMatrix.shape[0])==1.0].tolist()[0]
        return griddedModule.gridded( axes=newAxes, array=arrayModule.full(shape=commonMatrix.shape,data=trigdata,symmetry=arrayModule.symmetryLowerToken) )

    def toXMLList( self, indent = '', **kwargs ) :
        """

        :param indent:
        :param kwargs:
        :return:
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s' % ( indent, self.moniker ) ]
        if self.label is not None: xmlString[0] += ' label="%s"' % self.label
        lowerBound = self.lowerBound.toString( keepPeriod = False )
        upperBound = self.upperBound.toString( keepPeriod = False )
        xmlString[0] += ' lowerBound="%s" upperBound="%s">' % ( lowerBound, upperBound )
        xmlString.append( indent2 + '<!-- The matrix for this reaction equals the weighted sum of the following matrices: -->' )
        for pointer in self.pointerList:
            xmlString.append( pointer.toXML( indent2, **kwargs ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString
        
    def getReferredCovariance( self, pointer ):
        """

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
        #xys = [ self.getReferredCovariance( p ).getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        xys = [ p.link['eval'].getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
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
        lower = PQU.PQU( element.get("lowerBound") )
        upper = PQU.PQU( element.get("upperBound") )
        summed_ = summedCovariance( label = label, lowerBound=lower, upperBound=upper )
        for summandElement in element:
            link_ = summand.parseXMLNode( summandElement, xPath, linkData )
            link_['coefficient'] = float( link_['coefficient'] )
            linkData['unresolvedLinks'].append( link_ )
            summed_.pointerList.append( link_ )
        xPath.pop()
        return summed_

class summand( linkModule.link ):
    moniker = 'summand'
