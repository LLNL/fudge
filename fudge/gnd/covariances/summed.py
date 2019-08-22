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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

from fudge.core.ancestry import ancestry
from . import tokens
from fudge.core.math import matrix as gndMatrix
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

class summedCovariance( ancestry ):
    """ 
    Covariance matrix stored as sum/difference of other matrices. 
    """

    moniker = tokens.summedFormToken

    def __init__(self, index=None, lowerBound=None, upperBound=None, 
            coefficients=None, pointerList=None):
        ancestry.__init__( self, tokens.summedFormToken, None, attribute = 'index' )
        self.index = index #: an int, only needed within 'mixed' section
        self.lowerBound = lowerBound #: lower bound of row/column direction.  (type: PhysicalQuantityWithUncertainty)
        self.upperBound = upperBound #: upper bound of row/column direction.  (type: PhysicalQuantityWithUncertainty)
        self.pointerList = pointerList or [] #a list of pointers (type: fudge.gnd.link).  Each pointer has a dictionary with entry "coefficient" to use to weight the covariance referred to

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def toCorrelationMatrix( self ):
        return self.toCovarianceMatrix().toCorrelationMatrix()

    def toCovarianceMatrix( self ): 
        '''
        Sum the parts to construct the covariance matrix.  
        Note, each part must be converted to an absolute covariance before summing.
        '''
        if len( self.pointerList ) == 1: return self.pointerList[0].link.getNativeData().toCovarianceMatrix()
        from fudge.core.utilities.brb import uniquify
        import numpy, copy
        
        # set up common data using first element in pointerList
        firstCovMtx = self.pointerList[0].link.getNativeData().toCovarianceMatrix().toAbsolute()
        commonRowAxis = copy.copy( firstCovMtx.axes[0] )
        if firstCovMtx.axes[1].mirrorOtherAxis: commonColAxis = copy.copy( firstCovMtx.axes[0] )
        else:                                   commonColAxis = copy.copy( firstCovMtx.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.axes[2] )
        commonType = firstCovMtx.type
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for p in self.pointerList[1:]:
            cc = p.link.getNativeData().toCovarianceMatrix().toAbsolute() # a little recursion to take care of nested covariances
            if cc.type != commonType: raise ValueError( "Incompatable types in "+str(self.__class__)+": "+str(commonType)+' vs. '+str(cc.type) )
            cc.convertAxesToUnits( ( commonRowAxis.unit, commonColAxis.unit, commonMatrixAxis.unit ) )
            commonRowAxis.data = commonRowAxis.data + cc.axes[0].data
            if cc.axes[1].mirrorOtherAxis: commonColAxis.data = commonColAxis.data + cc.axes[0].data
            else:                          commonColAxis.data = commonColAxis.data + cc.axes[1].data
        commonRowAxis.data.sort()
        commonColAxis.data.sort()
        commonRowAxis.data = uniquify( commonRowAxis.data )
        commonColAxis.data = uniquify( commonColAxis.data )
        
        # now sum up the components
        commonMatrix = self.pointerList[0]['coefficient'] * numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for p in self.pointerList[1:]:
            cc = p.link.getNativeData().toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += p['coefficient'] * numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        return covarianceMatrix( type=commonType, axes=[ commonRowAxis, commonColAxis, commonMatrixAxis ], matrix=gndMatrix.matrix( commonMatrix.tolist(), form=gndMatrix.symmetricFormToken ) )
        

    def toXMLList(self, flags=None, indent=''):
        indent2 = indent+'  '
        xmlString = [ indent+'<%s' % self.moniker ]
        if self.index!=None: xmlString[0] += ' index="%s"' % self.index
        xmlString[0] += ' lowerBound="%s" upperBound="%s">' % (
                self.lowerBound, self.upperBound)
        xmlString.append(indent2+'<!-- The matrix for this reaction equals the weighted sum of the following matrices: -->')
        for pointer in self.pointerList:
            xmlString.append( pointer.toXML( indent2 ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString
        
    def getReferredCovariance( self, pointer ):
        if pointer.link is not None: return pointer.link
        if 'covarianceSuite' in pointer.path:   return link.follow( pointer.path, self.getRootParent() )
        #elif'reactionSuite' in pointer.path:    return link.follow( pointer.path, None )
        else :                                  raise ValueError( "Need reference to root node of "+str( pointer.path ) )

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combine all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs instance 
        """
        #xys = [ self.getReferredCovariance( p ).getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        xys = [ p.link.getNativeData().getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        uncert = coefficients[0] * xys[0]
        for i, xy in enumerate( xys[1:] ):
            uncert, xy = uncert.mutualify( 1e-8, 0, 0, xy, 1e-8, 0, 0 )
            uncert += coefficients[ i ] * xy
        return uncert

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        from fudge.gnd import link
        index = element.get("index")
        if index is not None: index=int(index)
        lower = PhysicalQuantityWithUncertainty( element.get("lowerBound") )
        upper = PhysicalQuantityWithUncertainty( element.get("upperBound") )
        summed_ = summedCovariance( index, lower, upper )
        for summand in element:
            link_ = link.parseXMLNode( summand, linkData )
            link_['coefficient'] = float( link_['coefficient'] )
            linkData['unresolvedLinks'].append( link_ )
            summed_.pointerList.append( link_ )
        xPath.pop()
        return summed_

    def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
        from fudge.legacy.converting import endfFormats
        endf = []
        if not inCovarianceGroup:
            # print header for this subsection (contains one NL sub-subsection)
            XMF1,XLFS1,MAT1,NC,NI = 0,0,0,1,0
            rowdat, coldat = targetInfo['dataPointer']
            MT1 = map(int, rowdat['ENDF_MFMT'].split(',')) [1]
            endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        # header:
        LTY=0
        endf.append( endfFormats.endfHeadLine(0,0,0,LTY,0,0) )
        NCI = len(self.pointerList)
        endf.append( endfFormats.endfHeadLine(self.lowerBound.getValueAs('eV'),
            self.upperBound.getValueAs('eV'),0,0, 2*NCI,NCI) )
        mtList = [ map(int, a['ENDF_MFMT'].split(','))[1] for a in self.pointerList]
        coefficients = [a['coefficient'] for a in self.pointerList]
        endf += endfFormats.endfDataList( [i for j in zip(coefficients,mtList)
            for i in j] )
        return endf
