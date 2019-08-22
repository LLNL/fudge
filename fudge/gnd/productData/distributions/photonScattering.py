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
Photon coherent and incoherent scattering factors.
"""

import xData.ancestry as ancestryModule
import xData.standards as standardsModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule

import base as baseModule

class scatteringFactor :        # This class is designed to look like a module and not a real class.

    class XYs1d( XYsModule.XYs1d ) :

        def __init__( self, **kwargs ) :

            XYsModule.XYs1d.__init__( self, **kwargs )

        def toXMLList( self, indent = "", **kwargs ) :

            return( XYsModule.XYs1d.toXMLList( self, indent, **kwargs ) )

    class regions1d( regionsModule.regions1d ) :

        def __init__( self, **kwargs ) :

            regionsModule.regions1d.__init__( self, **kwargs )

            def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

                accuracy = 1e-3
                if( len( self ) > 0 ) : accuracy = self[0].getAccuracy( )
                xys = regionsModule.regions1d.toPointwise_withLinearXYs( self, accuracy, lowerEps, upperEps )
                return( scatteringFactor.XYs1d( xys.axes, xys, accuracy = xys.getAccuracy( ) ) )

    @staticmethod
    def defaultAxes( factorLabel, energyUnit = 'eV' ) :

        axes = axesModule.axes( rank = 2 )
        axes[0] = axesModule.axis( factorLabel, 0, "" )
        axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
        return( axes )

class coherentFunctionBase( ancestryModule.ancestry ) :

    def __init__( self, data ) :

        if( not( isinstance( data, ( scatteringFactor.XYs1d, scatteringFactor.regions1d ) ) ) ) :
            raise TypeError( "Needed %s scatteringFactor subform." % self.moniker )

        self.data = data
        data.setAncestor( self )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

class scatteringFunction( coherentFunctionBase ) :

    moniker = 'scatteringFactor'
    ENDFMT = 502

class imaginaryAnomalousFactor( coherentFunctionBase ) :

    moniker = 'imaginaryAnomalousFactor'
    ENDFMT = 505

class realAnomalousFactor( coherentFunctionBase ) :

    moniker = 'realAnomalousFactor'
    ENDFMT = 506

class coherent :

    class form( baseModule.form ) :

        moniker = 'coherentScattering'
        subformAttributes = ( 'formFactor', 'anomalousScatteringFactor_realPart', 'anomalousScatteringFactor_imaginaryPart' )

        def __init__( self, label, productFrame, formFactor, realPart, imaginaryPart, makeCopy = True ) :

            if( not( isinstance( formFactor, scatteringFunction ) ) ) :
                raise Exception( 'Instance is class "%s" and not scatteringFunction' % formFactor.__class__ )

            if( realPart is not None ) :
                if( not( isinstance( realPart, realAnomalousFactor ) ) ) :
                    raise Exception( 'Instance is class "%s" and not realAnomalousFactor' % realPart.__class__ )

            if( imaginaryPart is not None ) :
                if( not( isinstance( imaginaryPart, imaginaryAnomalousFactor ) ) ) :
                    raise Exception( 'Instance is class "%s" and not imaginaryAnomalousFactor' % imaginaryPart.__class__ )

            baseModule.form.__init__( self, label, productFrame, ( formFactor, realPart, imaginaryPart ) )

        def calculateAverageProductData( self, style, indent = '', **kwargs ) :

            raise Exception( 'Not supported for %s' % self.moniker )

        def toXMLList( self, indent = "", **kwargs ) :

            indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

            attributeStr = ''
            if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
            if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
            xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
            if( self.formFactor is not None ) : xmlString += self.formFactor.toXMLList( indent2, **kwargs )
            if( self.anomalousScatteringFactor_realPart is not None ) :
                xmlString += self.anomalousScatteringFactor_realPart.toXMLList( indent2, **kwargs )
            if( self.anomalousScatteringFactor_imaginaryPart is not None ) :
                xmlString += self.anomalousScatteringFactor_imaginaryPart.toXMLList( indent2, **kwargs )
            xmlString[-1] += '</%s>' % self.moniker
            return( xmlString )

class incoherent :

    class form( baseModule.form ) :

        moniker = 'incoherentScattering'
        subformAttributes = ( 'scatteringFunction', )

        def __init__( self, label, productFrame, scatteringFunction, makeCopy = True ) :

            if( not( isinstance( scatteringFunction, ( incoherent.XYs1d, incoherent.regions1d ) ) ) ) :
                raise Exception( 'Instance is class "%s" and not scatteringFunction' % scatteringFunction.__class__ )
            if( makeCopy ) : scatteringFunction = scatteringFunction.copy( )

            baseModule.form.__init__( self, label, productFrame, ( scatteringFunction, ) )

        def calculateAverageProductData( self, style, indent = '', **kwargs ) :

            raise Exception( 'Not supported for %s' % self.moniker )

    class XYs1d( scatteringFactor.XYs1d ) :

        ENDFMT = 504

        def __init__( self, **kwarg ) :

            scatteringFactor.XYs1d.__init__( self, **kwarg )

    class regions1d( scatteringFactor.regions1d ) :

        ENDFMT = 504

        def __init__( self, **kwarg ) :

            scatteringFactor.regions1d.__init__( self, **kwarg )
