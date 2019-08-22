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

""" angular distribution classes """

import math
import fudge
import miscellaneous
from fudge.core.math import matrix
from fudge.core.utilities import brb
from fudge.gnd import tokens as tokensModule
from fudge.core.math import fudgemath

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule

noExtrapolationToken = 'noExtrapolation'

__metaclass__ = type

class subform( baseModule.subform ) :
    """Abstract base class for angular forms."""

    pass

class pointwise( subform, multiD_XYsModule.multiD_XYs ) :

    def __init__( self, axes, productFrame ) :

        raise 'this is probably no longer used. If it is used see Bret.'
        multiD_XYsModule.multiD_XYs.__init__( self, axes = axes )
        subform.__init__( self, productFrame )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

# Do we handle corresponding points (see ENDF n-009_F_019.endf MT=16)?????????????

        def calculateAverageEnergy( self, Ein ) :
            """For internal use only."""

            W_XYs_LegendreSeries = self.energyAngularForm.getW_XYs_LegendreSeriesAtV( Ein )
            EpC0s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries[0] ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            EpC1s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries.getCoefficientSafely( 1 ) ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            averageEp = XYsModule.XYs( EpC0s, axes = XYsModule.XYs.defaultAxes( ), accuracy = 1e-6 ).integrateWithWeight_x( )

            if( self.productFrame == standardsModule.frames.centerOfMassToken ) :
                A_Ein = self.massRatio * Ein
                averageSqrt_Ep = XYsModule.XYs( EpC0s, axes = XYsModule.XYs.defaultAxes( ), accuracy = 1e-6 ).integrateWithWeight_sqrt_x( )
                averageEp += A_Ein + 2. * math.sqrt( A_Ein ) * averageSqrt_Ep

            return( self.multiplicity.getValue( Ein ) * averageEp )

        def calculateAverageMomentum( self, Ein ) :
            """For internal use only."""

            W_XYs_LegendreSeries = self.energyAngularForm.getW_XYs_LegendreSeriesAtV( Ein )
            EpC1s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries.getCoefficientSafely( 1 ) ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            averagePp = XYsModule.XYs( EpC1s, axes = XYsModule.XYs.defaultAxes( ), accuracy = 1e-6 ).integrateWithWeight_sqrt_x( )
            if( self.productFrame == standardsModule.frames.centerOfMassToken ) : averagePp += math.sqrt( self.massRatio * Ein )
            return( multiplicity.getValue( E ) * math.sqrt( 2. * self.productMass ) * averagePp )

        class calculateDepositionInfo :
            """For internal use only."""

            def __init__( self, productFrame, productMass, massRatio, energyAngularForm, multiplicity ) :

                self.productFrame = productFrame
                self.productMass = productMass
                self.massRatio = massRatio
                self.energyAngularForm = energyAngularForm
                self.multiplicity = multiplicity

            def evaluateAtX( self, E ) :

                return( self._evaluateAtX( self, E ) )

            def setEvaluateAtX( self, evaluateAtX ) :

                self._evaluateAtX = evaluateAtX

            def setTolerances( self, relativeTolerance, absoluteTolerance ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

        energyUnit = self.axes[0].getUnit( )
        momentumDepositionUnit = energyPrimeUnit + '/c'
        massUnit = energyPrimeUnit + '/c**2'
        energyPrimeUnit = self.axes[1].getUnit( )
        energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy

        depData = []
        multiplicity = tempInfo['multiplicity']

        productFrame = self.getProductFrame( )
        projectileMass = tempInfo['reactionSuite'].projectile.getMass( massUnit )
        targetMass = tempInfo['reactionSuite'].target.getMass( massUnit )
        productMass = tempInfo['product'].getMass( massUnit )
        massRatio = projectileMass * productMass / ( projectileMass + targetMass )**2

        calculationData = calculateDepositionInfo( productFrame, productMass, massRatio, self, multiplicity )

        calculationData.setEvaluateAtX( calculateAverageEnergy )
        Es = [ W_XY_Ls.value for W_XY_Ls in self ]
        depEnergy = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in depEnergy ] )
        calculationData.setTolerances( energyAccuracy, absoluteTolerance )
        depEnergy = fudgemath.thickenXYList( depEnergy, calculationData )
        axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyPrimeUnit )
        depData.append( energyDepositionModule.pointwise( data = depEnergy, axes = axes,
                label = processInfo.style.label, accuracy = energyAccuracy ) )

        calculationData.setEvaluateAtX( calculateAverageMomentum )
        depMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ Ep for E, Ep in depMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        depMomentum = fudgemath.thickenXYList( depMomentum, calculationData )
        axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
                label = processInfo.style.label, accuracy = momentumAccuracy ) )

        return( depData )

    def check( self, info ) :
        from fudge.gnd import warning
        from pqu import PQU
        warnings = []
        axes = axesModule.axes()
        axes[0] = self.axes[0]
        axes[1] = self.axes[1]
#        unitbase = self.axes[0].interpolation.isQualifierUnitBase()

        for energy_in in self:
            integral = XYsModule.XYs( [ ( eout.value, eout.coefficients[0] ) for eout in energy_in ], 
                    axes = axes, accuracy = 0.001 ).integrate( )
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU(energy_in.value,axes[-1].unit),
                    energy_in.index, integral, self.toXLink() ) )
            minProb = min( [xy.toPointwise_withLinearXYs( 0.001 ).rangeMin() for xy in energy_in] )
            if minProb < 0:
                warnings.append( warning.negativeProbability( PQU.PQU(energy_in.value,axes[-1].unit),
                    value = minProb, obj = energy_in ) )
            if unitbase:    # check for more than one L0 coefficient = 0 at the high end of distribution
                LeqZeroCoefs = [eout[0] for eout in energy_in]
                for idx, coef in enumerate(LeqZeroCoefs[::-1]):
                    if coef != 0: break
                if idx > 1:
                    warnings.append( warning.extraOutgoingEnergy( PQU.PQU(energy_in.value,axes[-1].unit),
                        obj = energy_in ) )
        return warnings

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        newComponents = []
        crossSection = tempInfo['crossSection']
        multiplicity = tempInfo['multiplicity']
        projectile, target, product = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target, tempInfo['product'].particle
        projectileName, productGroupName = projectile.name, product.name

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            TM_1, TM_E = transferMatrices.Legendre_TransferMatrix( processInfo, projectileName, productGroupName, crossSection, self, multiplicity,
                tempInfo['masses'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        import energyAngular

        return( LegendreSeries.V_W_XYs_LegendreSeries.toPointwise_withLinearXYs( self, accuracy, cls = energyAngular.pointwise ) ) 

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent2, **kwargs )
        for index, data in enumerate( self ) : xmlString += data.toXMLList( tag = self.axes[0].getLabel( ), indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energy_outUnit = 'eV', C_lUnit = '1/eV' ) :

        axes = axesModule.axes( rank = 3 )
        axes[0] = axesModule.axis( 'C_l',        0, C_lUnit )
        axes[1] = axesModule.axis( 'energy_out', 1, energy_outUnit )
        axes[2] = axesModule.axis( 'energy_in',  2, energyUnit )
        return( axes )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        axes = axesModule.parseXMLNode( element[0], xPath )
        pw = pointwise( axes, element.get( 'productFrame' ) )
        for energy_in in element[1:]:
            w_xys = LegendreSeries.W_XYs_LegendreSeries( value = float( energy_in.get( "value" ) ) )
            for energy_out in energy_in:
                xys = LegendreSeries.XYs_LegendreSeries( map( float, energy_out.text.split() ),
                        value = float( energy_out.get( "value" ) ), parent = w_xys )
                w_xys.append( xys, copy=False )
            pw.append( w_xys, copy=False )
        xPath.pop()
        return pw

class LLNLPointwise( subform ) :
    """This is a temporary class, to be removed once testing is done and all coefficients in endl99, H1(n,2n)p data are filled in."""

    moniker = 'LLNLLegendrePointwise'

    def __init__( self, axes ) :

        subform.__init__( self )
        self.data = []
        self.axes = axes
        self.label = None

    def __getitem__( self, l ) :

        return( self.data[l] )

    def __len__( self ) :

        return( len( self.data ) )

    def append( self, EEpP ) :

        if( not( isinstance( EEpP, multiD_XYsModule.multiD_XYs ) ) ) : raise Exception( 'EEpP is an instance of %s' % brb.getType( EEpP ) )
        self.data.append( EEpP )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        energyUnit = self.axes[1].getUnit( )
        momentumDepositionUnit = energyPrimeUnit + '/c'
        energyPrimeUnit = self.axes[2].getUnit( )
        energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy

        depData = []
        multiplicity = tempInfo['multiplicity']

        Legendre_l0, Legendre_l1 = self[0], None
        if( len( self ) > 1 ) : Legendre_l1 = self[1]

        if( isinstance( multiplicity, multiplicityModule.pointwise ) ) :
            Es = [ EpP.value for EpP in Legendre_l0 ]
            doIt = False
            for E, m in multiplicity.getFormByToken( tokensModule.pointwiseFormToken ) :
                if( E < Es[0] ) : continue
                if( E not in Es ) :
                    doIt = True
                    Es.append( E )
            if( doIt ) :
                Es.sort( )
                Legendre_l0p, Legendre_l1p = LLNLPointwiseEEpP( 0 ), LLNLPointwiseEEpP( 1 )
                indexE, EpP1 = 0, None
                for EpP2 in Legendre_l0 :
                    while( ( indexE < len( Es ) ) and ( EpP2.value > Es[indexE] ) ) :
                        if( EpP1.value != Es[indexE] ) :
                            EpP = XYs.pointwiseXY_C.unitbaseInterpolate( Es[indexE], EpP1.value, EpP1, EpP2.value, EpP2 )
                            axes = EpP1.axes.copy( )
                            Legendre_l0p.append( XYs.XYs( axes, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
                        indexE += 1
                    Legendre_l0p.append( EpP2 )
                    EpP1 = EpP2
                Legendre_l0 = Legendre_l0p
                if( Legendre_l1 is not None ) :
                    for EpP2 in Legendre_l1 :
                        while( ( indexE < len( Es ) ) and ( EpP2.value > Es[indexE] ) ) :
                            if( EpP1 is None ) : raise Exception( 'multiplicity energy = %s before Legendre l = 0 distribution energy = %s' 
                                % ( Es[indexE], EpP2.value ) )
                            if( EpP1.value != Es[indexE] ) :
                                EpP = XYs.pointwiseXY_C.unitbaseInterpolate( Es[indexE], EpP1.value, EpP1, EpP2.value, EpP2 )
                                axes = EpP1.axes.copy( )
                                Legendre_l1p.append( XYs.XYs( axes, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
                            indexE += 1
                        Legendre_l1p.append( EpP2 )
                        EpP1 = EpP2
                    Legendre_l1 = Legendre_l1p

        depEnergy = [ [ EpP.value, multiplicity.getValue( EpP.value ) * miscellaneous.calculateDepositionEnergyFromEpP( EpP.value, EpP ) ] for EpP in Legendre_l0 ]
        if( depEnergy[0][0] > tempInfo['EMin'] ) : depEnergy.insert( 0, [ tempInfo['EMin'], 0. ] ) # Special case for bad data

        axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyPrimeUnit )
        depData.append( energyDepositionModule.pointwise( data = depEnergy, axes = axes,
                label = processInfo.style.label, accuracy = energyAccuracy ) )

        if( Legendre_l1 is not None ) :
            const = math.sqrt( 2. * tempInfo['product'].getMass( energyPrimeUnit + '/c**2' ) )
            if( const == 0. ) : const = 1               # For gammas.
            depMomentum = []
            EpP = Legendre_l1[0]
            for EpP in Legendre_l1 :
                if( const == 0. ) :                     # For gammas.
                    depMomentum.append( [ EpP.value, const * multiplicity.getValue( EpP.value ) * EpP.integrateWithWeight_x( ) ] )
                else :
                    depMomentum.append( [ EpP.value, const * multiplicity.getValue( EpP.value ) * EpP.integrateWithWeight_sqrt_x( ) ] )
        else :
            depMomentum = [ [ Legendre_l0[0].value, 0. ], [ Legendre_l0[-1].value, 0. ] ]
        axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        if( depMomentum[0][0] > tempInfo['EMin'] ) : depMomentum.insert( 0, [ tempInfo['EMin'], 0. ] ) # Special case for bad data
        depData.append( momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
                label = processInfo.style.label, accuracy = momentumAccuracy ) )
        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        import energy
        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        Legendre = self.data

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.name
            TM_1, TM_E = transferMatrices.ELEpP_TransferMatrix( processInfo, projectile, product, tempInfo['crossSection'],
                Legendre, tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent2, **kwargs )
        for l in self : xmlString += l.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        for lOrderData in element :
            if( lOrderData.tag == axesModule.axes.moniker ) :
                axes = axesModule.axes.parseXMLNode( element[0], xPath )
                break
        lpw = LLNLPointwise( axes )
        for lOrderData in element :
            if( lOrderData.tag == axesModule.axes.moniker ) : continue
            lpw.append( multiD_XYsModule.multiD_XYs.parseXMLNode( lOrderData ) )
        xPath.pop( )
        return( lpw )

    @staticmethod
    def defaultAxes( ) :

        axes = axesModule.axes( rank = 4 )
        axes[0] = axesModule.axis( 'c_l',        0, '1/MeV' )
        axes[1] = axesModule.axis( 'energy_out', 1, 'MeV' )
        axes[2] = axesModule.axis( 'energy_in',  2, 'MeV' )
        axes[3] = axesModule.axis( 'l',          3, '' )
        return( axes )

class LLNLPointwiseEEpP :

    def __init__( self, l ) :

        raise "this should no longer be used, but should not be removed until LLNLPointwise.calculateDepositionData is changed as it uses it."
        self.l = l
        self.EpP = []

    def __getitem__( self, index ) :

        return( self.EpP[index] )

    def __len__( self ) :

        return( len( self.EpP ) )

    def append( self, EpP ) :

        if( not( isinstance( EpP, XYsModule.XYs ) ) ) : raise Exception( 'EpP is an instance of %s' % brb.getType( EpP ) )
        if( len( self ) > 0 ) :
            if( EpP.value is None ) : raise Exception( "EpP's value is None" )
            if( EpP.value <= self.EpP[-1].value ) : raise Exception( "EpP.value = %s <= prior's values = %s " % ( EpP.value, self.EpP[-1].value ) )
        self.EpP.append( EpP )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        kwargs['oneLine'] = True

        xmlString = [ '%s<l index="%s">' % ( indent, self.l ) ]
        for energy in self : xmlString += energy.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</l>'
        return( xmlString )

class groupedBase( subform ) :

    xData = matrix.matrixToken

    def __init__( self, axes, transferMatrices, productFrame ) :

        subform.__init__( self, base.groupedFormToken, productFrame )
        self.axes = axes
        self.transferMatrices = transferMatrices

    def __len__( self ) :

        return( len( self.transferMatrices ) )

    def __getitem__( self, index ) :

        return( self.transferMatrices[index] )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        energyEpCl = self.transferMatrices[0]
        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = ' size="%d,%d"' % ( energyEpCl.nrows, energyEpCl.ncols ) ) ]
        xmlString += self.axes.toXMLList( indent2, **kwargs )
        for l, energyEpCl in enumerate( self.transferMatrices ) :
            xmlString.append( '%s<l value="%d">' % ( indent2, l ) )
            xmlString.extend( energyEpCl.toXMLList( indent3, **kwargs ) )
            xmlString[-1] += '</l>'
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        axes = axesModule.parseXMLNode( element[0], xPath )
        frame = element.get( 'productFrame' )
        transferMatrices = []
        for lvalue in element[1:]:
            transferMatrices.append( matrix.parseXMLNode( lvalue[0], xPath, linkData ) )
        GB = cls( axes, transferMatrices, frame )
        xPath.pop()
        return GB

class grouped( groupedBase ) :

    def __init__( self, axes, transferMatrices, productFrame ) :

        groupedBase.__init__( self, axes, transferMatrices, productFrame )

class energyConservationComponent( baseModule.form ) :

    moniker = 'LegendreEnergyConservation'

    def __init__( self ) :

        baseModule.form.__init__( self )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        LC = energyConservationComponent( nativeData = element.get( 'nativeData' ) )
        for subformElement in element:
            subformClass = {
                    baseModule.groupedFormToken: energyConservationGrouped,
                }.get( subformElement.tag )
            if subformClass is not None:
                LC.addForm( subformClass.parseXMLNode( subformElement, xPath, linkData ) )
            else: raise Exception( "encountered unknown LegendreEnergyConservation subform: %s" % subformElement.tag )
        xPath.pop()
        return LC

class energyConservationGrouped( groupedBase ) :

    def __init__( self, axes, transferMatrices, productFrame ) :

        groupedBase.__init__( self, axes, transferMatrices, productFrame )

class form( baseModule.form ) :

    moniker = 'Legendre'
    subformAttributes = ( 'LegendreSubform', )

    def __init__( self, label, productFrame, LegendreSubform, makeCopy = True ) :

        if( not( isinstance( LegendreSubform, subform ) ) ) : raise TypeError( 'instance is not a Legendre subform' )
#        if( makeCopy ) : LegendreSubform = LegendreSubform.copy()
        baseModule.form.__init__( self, label, productFrame, ( LegendreSubform, ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        raise 'hell - I do not think this is used anymore - BRB'
        return( self.forms[0].calculateDepositionData( processInfo, tempInfo, verbosityIndent ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.forms[0].process( processInfo, tempInfo, verbosityIndent ) )

    @staticmethod
    def parseXMLNode( LegendreElement, xPath, linkData ) :

        xPath.append( LegendreElement.tag )
        subformElement = LegendreElement[0]
        subformClass = {
                LLNLPointwise.moniker : LLNLPointwise,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Legendre subform: %s" % subformElement.tag )
        LegendreSubForm = subformClass.parseXMLNode( subformElement, xPath, linkData )
        LegendreComponent = component( LegendreElement.get( 'label' ), LegendreElement.get( 'productFrame' ), LegendreSubform )
        xPath.pop( )
        return( LegendreComponent )
