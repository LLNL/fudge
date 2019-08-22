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
from fudge.core.math import matrix as gndMatrix
from . import tokens, covarianceMatrix
from .summed import summedCovariance


class mixedForm( ancestry ):
    """
    Covariance for a single quantity, stored as several separate matrices that must be summed together.
    In general, the energy bounds for these matrices can overlap (unlike piecewise cross section data). 
    """

    moniker = tokens.mixedFormToken

    def __init__(self, components=None):
        ancestry.__init__( self, 'mixedForm', None )
        self.components = components or [] #: a Python list containing instances of ``mixedForm``, ``summedCovariance``, and ``covarianceMatrix``

    def __getitem__(self, idx):
        return self.components[idx]

    def __len__(self):
        return len(self.components)

    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        for component in self.components:
            componentWarnings = component.check( info )
            if componentWarnings:
                warnings.append( warning.context('Component %i' % component.index, componentWarnings) )
        
        #if warnings:
        #    warnings = [warning.context('Section "%s": %s' % (self.id, form), warnings)]
        return warnings
        
    def fix( self, **kw ): 
        warnings = []
        for comp in self.components: warnings += comp.fix( **kw )
        return warnings

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def addComponent(self, covariance):
        ''':param covariance: an instance of covariance (or inherited class)'''
        covariance.setParent(self)
        self.components.append(covariance)

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combines all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs instance
        """

        xys = [ sec.getUncertaintyVector( theData=theData, relative=relative ) for sec in self.components ]
        uncert = xys[0]
        for xy in xys[1:]:
            uncert, xy = uncert.mutualify( 1e-8, 1e-8, 0, xy, 1e-8, 1e-8, 0 )
            uncert += xy
        return uncert
        
    def toCovarianceMatrix( self ): 
        if len( self.components ) == 1: return self.components[0].toCovarianceMatrix()
        from fudge.core.utilities.brb import uniquify
        import numpy, copy
        
        # set up common data using first component
        firstCovMtx = self.components[0].toCovarianceMatrix()
        commonRowAxis = copy.copy( firstCovMtx.axes[0] )
        if firstCovMtx.axes[1].mirrorOtherAxis: commonColAxis = copy.copy( firstCovMtx.axes[0] )
        else:                                   commonColAxis = copy.copy( firstCovMtx.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.axes[2] )
        commonType = firstCovMtx.type
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
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
        commonMatrix = numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        return covarianceMatrix( type=commonType, axes=[ commonRowAxis, commonColAxis, commonMatrixAxis ], matrix=gndMatrix.matrix( commonMatrix.tolist(), form=gndMatrix.symmetricFormToken ) )

    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        import copy
        result = copy.copy( self )
        result.components = []
        for c in self.components:
            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toAbsolute() call with the correct row or column data
                result.components.append( c.toCovarianceMatrix().toAbsolute( rowData, colData ) )
            else: result.components.append( c.toAbsolute( rowData, colData ) )
        return result
        
    def toRelative( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a absolute covariance) using XYs rowData and colData
        to convert self into a relative covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        import copy
        result = copy.copy( self )
        result.components = []
        for c in self.components:
            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toRelative() call with the correct row or column data
                result.components.append( c.toCovarianceMatrix().toRelative( rowData, colData ) )
            else: result.components.append( c.toRelative( rowData, colData ) )
        return result
        
        
    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<%s>' % self.moniker]
        for covariance in self.components:
            xmlString += covariance.toXMLList(flags, indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath=[], linkData={} ):
        """Translate <mixed> element from xml."""

        xPath.append( element.tag )
        mixed_ = cls()
        for child in element:
            formClass = {
                    tokens.covarianceFormToken: covarianceMatrix,
                    tokens.summedFormToken: summedCovariance,
                    }.get( child.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % child.tag)
            mixed_.addComponent( formClass.parseXMLNode( child, xPath, linkData ) )
        for i,form in enumerate( mixed_.components ): form.index = i
        xPath.pop()
        return mixed_
    
    def toENDF6(self, flags, targetInfo):
        from fudge.legacy.converting import endfFormats
        endf = []
        XMF1,XLFS1,MAT1 = 0,0,0
        NI = len([cov for cov in self.components if hasattr(cov,'matrix')])
        NC = len(self.components) - NI
        rowdat, coldat = targetInfo['dataPointer']
        MF, MT1 = map(int, rowdat['ENDF_MFMT'].split(','))
        if MF in (31,33):
            endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        for cov in self.components:
            endf += cov.toENDF6(flags, targetInfo, inCovarianceGroup=True)
        return endf
