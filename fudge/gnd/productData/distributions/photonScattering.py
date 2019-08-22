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

"""
Photon coherent and incoherent scattering factors.
"""

import base as baseModule

import xData.standards as standardsModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule

class scatteringFactor :

    class pointwise( XYsModule.XYs ) :

        def __init__( self, **kwargs ) :

            XYsModule.XYs.__init__( self, **kwargs )

        def toXMLList( self, indent = "", **kwargs ) :

            return( XYsModule.XYs.toXMLList( self, indent, **kwargs ) )

        @staticmethod
        def defaultAxes( factorLabel, energyUnit = 'eV' ) :

            axes = axesModule.axes( rank = 2 )
            axes[0] = axesModule.axis( factorLabel, 0, "" )
            axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
            return( axes )

    class piecewise( regionsModule.regions ) :

        def __init__( self, **kwargs ) :

            regionsModule.regions.__init__( self, **kwargs )

            def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

                accuracy = 1e-3
                if( len( self ) > 0 ) : accuracy = self[0].getAccuracy( )
                xys = regionsModule.regions.toPointwise_withLinearXYs( self, accuracy, lowerEps, upperEps )
                return( scatteringFactor.pointwise( xys.axes, xys, accuracy = xys.getAccuracy( ) ) )

        @staticmethod
        def defaultAxes( factorLabel, energyUnit = 'eV' ) :

            axes = axesModule.axes( rank = 2 )
            axes[0] = axesModule.axis( 'energy_in', 0, energyUnit )
            axes[1] = axesModule.axis( factorLabel, 1, "" )
            return( axes )

class scatteringFunction :

    class pointwise( scatteringFactor.pointwise ) :

        ENDFMT = 502

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'scatteringFactor'
            scatteringFactor.pointwise.__init__( self, **kwargs )

    class piecewise( scatteringFactor.piecewise ) :

        ENDFMT = 502

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'scatteringFactor'
            scatteringFactor.piecewise.__init__( self, **kwargs )

class imaginaryAnomalousFactor :

    class pointwise( scatteringFactor.pointwise ) :

        ENDFMT = 505

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'imaginaryAnomalousFactor'
            scatteringFactor.pointwise.__init__( self, **kwargs )

    class piecewise( scatteringFactor.piecewise ) :

        ENDFMT = 505

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'imaginaryAnomalousFactor'
            scatteringFactor.piecewise.__init__( self, **kwargs )

class realAnomalousFactor :

    class pointwise( scatteringFactor.pointwise ) :

        ENDFMT = 506

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'realAnomalousFactor'
            scatteringFactor.pointwise.__init__( self, **kwargs )

    class piecewise( scatteringFactor.piecewise ) :

        ENDFMT = 506

        def __init__( self, **kwargs ) :

            kwargs['label'] = 'realAnomalousFactor'
            scatteringFactor.piecewise.__init__( self, **kwargs )

class coherent :

    class form( baseModule.form ) :

        moniker = 'coherentScattering'
        subformAttributes = ( 'formFactor', 'anomalousScatteringFactor_realPart', 'anomalousScatteringFactor_imaginaryPart' )

        def __init__( self, label, productFrame, formFactor, realPart, imaginaryPart, makeCopy = True ) :

            if( not( isinstance( formFactor, ( scatteringFunction.pointwise, scatteringFunction.piecewise ) ) ) ) :
                raise Exception( 'Instance is class "%s" and not scatteringFunction' % formFactor.__class__ )

            if( realPart is not None ) :
                if( not( isinstance( realPart, ( realAnomalousFactor.pointwise, realAnomalousFactor.piecewise ) ) ) ) :
                    raise Exception( 'Instance is class "%s" and not realAnomalousFactor' % realPart.__class__ )

            if( imaginaryPart is not None ) :
                if( not( isinstance( imaginaryPart, ( imaginaryAnomalousFactor.pointwise, imaginaryAnomalousFactor.piecewise ) ) ) ) :
                    raise Exception( 'Instance is class "%s" and not imaginaryAnomalousFactor' % imaginaryPart.__class__ )

            if( makeCopy ) :
                formFactor = formFactor.copy( )
                realPart = realPart.copy( )
                imaginaryPart = imaginaryPart.copy( )

            baseModule.form.__init__( self, label, productFrame, ( formFactor, realPart, imaginaryPart ) )

        def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

            raise Exception( 'Not supported for %s' % self.moniker )

        def process( self, processInfo, tempInfo, verbosityIndent ) :

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

            if( not( isinstance( scatteringFunction, ( incoherent.pointwise, incoherent.piecewise ) ) ) ) :
                raise Exception( 'Instance is class "%s" and not scatteringFunction' % scatteringFunction.__class__ )
            if( makeCopy ) : scatteringFunction = scatteringFunction.copy( )

            baseModule.form.__init__( self, label, productFrame, ( scatteringFunction, ) )

        def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

            raise Exception( 'Not supported for %s' % self.moniker )

        def process( self, processInfo, tempInfo, verbosityIndent ) :

            raise Exception( 'Not supported for %s' % self.moniker )

    class pointwise( scatteringFactor.pointwise ) :

        ENDFMT = 504

        def __init__( self, **kwarg ) :

            scatteringFactor.pointwise.__init__( self, **kwarg )

    class piecewise( scatteringFactor.piecewise ) :

        ENDFMT = 504

        def __init__( self, **kwarg ) :

            scatteringFactor.piecewise.__init__( self, **kwarg )
