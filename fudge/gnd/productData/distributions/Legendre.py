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

""" angular distribution classes """

import math
import fudge
import base, miscellaneous
from fudge.core.math.xData import axes, LegendreSeries, XYs
from fudge.core.math import matrix
from fudge.core.utilities import brb
from fudge.gnd import tokens
from fudge.core.math import fudgemath

noExtrapolationToken = 'noExtrapolation'

__metaclass__ = type

class component( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.LegendreComponentToken, nativeData )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        return( self.forms[self.nativeData].calculateDepositionData( processInfo, tempInfo ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.forms[self.nativeData].process( processInfo, tempInfo, verbosityIndent ) )

    def toENDF6( self, MT, endfMFList, flags, tempInfo ) :

        self.forms[self.nativeData].toENDF6( MT, endfMFList, flags, tempInfo )

    @staticmethod
    def parseXMLNode( LegendreElement, xPath=[], linkData={} ):

        xPath.append( LegendreElement.tag )
        LC = component( nativeData = LegendreElement.get('nativeData') )
        for formElement in LegendreElement:
            formClass = {
                    base.LegendrePointwiseFormToken: pointwise,
                    base.LLNLLegendrePointwiseFormToken: LLNLPointwise,
                    base.groupedFormToken: grouped,
                    }.get( formElement.tag )
            if formClass is not None:
                LC.addForm( formClass.parseXMLNode( formElement, xPath, linkData ) )
            else: raise Exception(" encountered unknown Legendre form: %s" % formElement.tag)
        xPath.pop()
        return LC

class form( base.form ) :
    """Abstract base class for Legendre forms."""

    pass

class pointwise( form, LegendreSeries.V_W_XYs_LegendreSeries ) :

    moniker = base.LegendrePointwiseFormToken

    def __init__( self, axes, productFrame ) :

        form.__init__( self, base.LegendrePointwiseFormToken, productFrame )
        LegendreSeries.V_W_XYs_LegendreSeries.__init__( self, axes )

    def calculateDepositionData( self, processInfo, tempInfo ) :

# Do we handle corresponding points (see ENDF n-009_F_019.endf MT=16)?????????????

        def calculateAverageEnergy( self, Ein ) :
            """For internal use only."""

            W_XYs_LegendreSeries = self.energyAngularForm.getW_XYs_LegendreSeriesAtV( Ein )
            EpC0s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries[0] ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            EpC1s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries.getCoefficientSafely( 1 ) ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            averageEp = XYs.XYs( XYs.XYs.defaultAxes( ), EpC0s, 1e-6 ).integrateWithWeight_x( )

            if( self.productFrame == axes.centerOfMassToken ) :
                A_Ein = self.massRatio * Ein
                averageSqrt_Ep = XYs.XYs( XYs.XYs.defaultAxes( ), EpC0s, 1e-6 ).integrateWithWeight_sqrt_x( )
                averageEp += A_Ein + 2. * math.sqrt( A_Ein ) * averageSqrt_Ep

            return( self.multiplicity.getValue( Ein ) * averageEp )

        def calculateAverageMomentum( self, Ein ) :
            """For internal use only."""

            W_XYs_LegendreSeries = self.energyAngularForm.getW_XYs_LegendreSeriesAtV( Ein )
            EpC1s = [ [ XYs_LegendreSeries.value, XYs_LegendreSeries.getCoefficientSafely( 1 ) ] for XYs_LegendreSeries in W_XYs_LegendreSeries ]
            averagePp = XYs.XYs( XYs.XYs.defaultAxes( ), EpC1s, 1e-6 ).integrateWithWeight_sqrt_x( )
            if( self.productFrame == axes.centerOfMassToken ) : averagePp += math.sqrt( self.massRatio * Ein )
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

        depData = []
        multiplicity = tempInfo['multiplicity']
        energyAccuracy, momentumAccuracy = 1e-4, 1e-2
        energyUnit = self.axes[0].getUnit( )
        energyPrimeUnit = self.axes[1].getUnit( )
        massUnit = energyPrimeUnit + '/c**2'
        momentumDepositionUnit = energyPrimeUnit + '/c'

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
        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyPrimeUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )

        calculationData.setEvaluateAtX( calculateAverageMomentum )
        depMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ Ep for E, Ep in depMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        depMomentum = fudgemath.thickenXYList( depMomentum, calculationData )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def check( self, info ) :
        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []
        ax = axes.axes()
        ax[0] = self.axes[1]; ax[1] = self.axes[2]
        unitbase = self.axes[0].interpolation.isQualifierUnitBase()

        for energy_in in self:
            integral = XYs.XYs( ax, [(eout.value, eout.coefficients[0]) for eout in energy_in], 0.001 ).integrate()
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU(energy_in.value,ax[0].unit),
                    energy_in.index, integral, self.toXLink() ) )
            minProb = min( [xy.toPointwise_withLinearXYs( 0.001 ).yMin() for xy in energy_in] )
            if minProb < 0:
                warnings.append( warning.negativeProbability( PQU(energy_in.value,ax[0].unit),
                    value = minProb, obj = energy_in ) )
            if unitbase:    # check for more than one L0 coefficient = 0 at the high end of distribution
                LeqZeroCoefs = [eout[0] for eout in energy_in]
                for idx, coef in enumerate(LeqZeroCoefs[::-1]):
                    if coef != 0: break
                if idx > 1:
                    warnings.append( warning.extraOutgoingEnergy( PQU(energy_in.value,ax[0].unit),
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
        projectileName, productGroupName = projectile.getName( ), product.getName( )

        if( 'LLNL_MC' in processInfo['styles'] ) :
            import energyAngular

            energyAngularComponentForm = self.toPointwise_withLinearXYs( 1e-3 )
            energyAngularComponent = energyAngular.component( energyAngularComponentForm.moniker )
            energyAngularComponent.addForm( energyAngularComponentForm )
            newComponents.append( energyAngularComponent )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            TM_1, TM_E = transferMatrices.Legendre_TransferMatrix( processInfo, projectileName, productGroupName, crossSection, self, multiplicity,
                tempInfo['masses'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        import energyAngular

        return( LegendreSeries.V_W_XYs_LegendreSeries.toPointwise_withLinearXYs( self, accuracy, cls = energyAngular.linear ) ) 

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for index, data in enumerate( self ) : xmlString += data.toXMLList( tag = self.axes[0].getLabel( ), indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, tempInfo ) :

        from fudge.legacy.converting import gndToENDF6, endfFormats
        EInInterpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[0] )
        independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
        if( dependent == axes.flatToken ) :
            LEP = 1  # interpolation for Eout
        else :
            LEP = 2
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 1, LEP, 1, len( self ) ) ]
        ENDFDataList += endfFormats.endfInterpolationList( [ len( self ), EInInterpolation ] )
        for energy_in in self :
            NA, data = 0, []
            for w_xys_LegendreSeries in energy_in :
                NA = max( len( w_xys_LegendreSeries) , NA )
                data += [ w_xys_LegendreSeries.value ] + w_xys_LegendreSeries.coefficients
            ENDFDataList.append( endfFormats.endfContLine( 0, energy_in.value, 0, NA - 1, len( data ), len( data ) / ( NA + 1 ) ) )
            ENDFDataList += endfFormats.endfDataList( data )
        LAW = 1
        gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, tempInfo, LAW, self.productFrame, ENDFDataList )

    @staticmethod
    def defaultAxes( energyInterpolation, energyFunctionInterpolation, energyInterpolationQualifier = None, energyUnit = 'eV',
        energy_outInterpolation = axes.linearToken, energy_outUnit = 'eV', C_lInterpolation = axes.linearToken, C_lUnit = '1/eV' ) :

        axes_ = axes.axes( dimension = 3 )
        axes_[0] = axes.axis( 'energy_in',  0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, 
            energyFunctionInterpolation, energyInterpolationQualifier ) )
        axes_[1] = axes.axis( 'energy_out', 1, energy_outUnit, interpolation = axes.interpolationXY( energy_outInterpolation, C_lInterpolation ) )
        axes_[2] = axes.axis( 'C_l',        2, C_lUnit )
        return( axes_ )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        axes_ = axes.parseXMLNode( element[0], xPath )
        pw = pointwise( axes_, element.get( 'productFrame' ) )
        for energy_in in element[1:]:
            w_xys = LegendreSeries.W_XYs_LegendreSeries(
                    axes.referenceAxes(pw), value=float( energy_in.get("value") ), parent=pw)
            for energy_out in energy_in:
                xys = LegendreSeries.XYs_LegendreSeries( map( float, energy_out.text.split() ),
                        value = float( energy_out.get( "value" ) ), parent = w_xys )
                w_xys.append( xys, copy=False )
            pw.append( w_xys, copy=False )
        xPath.pop()
        return pw

class LLNLPointwise( form ) :
    """This is a temporary class, to be removed once testing is done and all coefficients in endl99, H1(n,2n)p data are filled in."""

    xData = base.LLNLLegendrePointwiseFormToken

    def __init__( self, axes, productFrame ) :

        form.__init__( self, base.LLNLLegendrePointwiseFormToken, productFrame )
        self.data = []
        self.axes = axes

    def __getitem__( self, l ) :

        return( self.data[l] )

    def __len__( self ) :

        return( len( self.data ) )

    def append( self, EEpP ) :

        l = 0
        if( len( self.data ) > 0 ) : l = self.data[-1].l + 1
        if( not( isinstance( EEpP, LLNLPointwiseEEpP ) ) ) : raise Exception( 'EEpP is an instance of %s' % brb.getType( EEpP ) )
        if( EEpP.l != l ) : raise Exception( 'EEpP.l = %s != l = %s' % ( EEpP.l, l ) )
        self.data.append( EEpP )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        depData = []
        multiplicity = tempInfo['multiplicity']
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3
        energyUnit = self.axes[1].getUnit( )
        energyPrimeUnit = self.axes[2].getUnit( )
        momentumDepositionUnit = energyPrimeUnit + '/c'

        Legendre_l0, Legendre_l1 = self[0], None
        if( len( self ) > 1 ) : Legendre_l1 = self[1]

        if( multiplicity.getNativeDataToken( ) == tokens.pointwiseFormToken ) :
            Es = [ EpP.value for EpP in Legendre_l0 ]
            doIt = False
            for E, m in multiplicity.getFormByToken( tokens.pointwiseFormToken ) :
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
                            axes_ = EpP1.axes.copy( )
                            Legendre_l0p.append( XYs.XYs( axes_, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
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
                                axes_ = EpP1.axes.copy( )
                                Legendre_l1p.append( XYs.XYs( axes_, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
                            indexE += 1
                        Legendre_l1p.append( EpP2 )
                        EpP1 = EpP2
                    Legendre_l1 = Legendre_l1p

        depEnergy = [ [ EpP.value, multiplicity.getValue( EpP.value ) * miscellaneous.calculateDepositionEnergyFromEpP( EpP.value, EpP ) ] for EpP in Legendre_l0 ]
        if( depEnergy[0][0] > tempInfo['EMin'] ) : depEnergy.insert( 0, [ tempInfo['EMin'], 0. ] ) # Special case for bad data

        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyPrimeUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )

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
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        if( depMomentum[0][0] > tempInfo['EMin'] ) : depMomentum.insert( 0, [ tempInfo['EMin'], 0. ] ) # Special case for bad data
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )
        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        import energy
        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        Legendre = self.data
        if( 'LLNL_MC' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sEqual probably binning pointwise %s' % ( verbosityIndent, self.moniker )
            raise Exception( 'Not implemented: %s' % self.moniker )
            nBins = 32
            component = energy.component( nativeData = base.equalProbableBinsFormToken )
            E_Mus = [ [ E, fudge2dEqualProbableBinning.equalProbableBins( nBins, muP ) ] for E, muP in Legendre.data[0][1] ]
            variables = fudge.gnd.miscellaneous.copyVariables( Legendre.variables, [ 0, 1 ] )
            variables[1]['interpolation'] = 'epb'
            component.addForm( energy.equalProbableBins( variables, E_Mus ) )
            newComponents.append( component )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.getToken( )
            TM_1, TM_E = transferMatrices.ELEpP_TransferMatrix( processInfo, projectile, product, tempInfo['crossSection'],
                Legendre, tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    def toENDF6( self, MT, endfMFList, flags, tempInfo ) :

        from fudge.legacy.converting import gndToENDF6, endfFormats
        LAW, LEP = 1, 2
        EInInterpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[1] )
        E_ins = [ [ EpCl.value, {} ] for EpCl in self[0].EpP ]
        for l_EEpCl in self :
            for indexE, EpCl in enumerate( l_EEpCl ) :
                if( EpCl.value != E_ins[indexE][0] ) :
                    raise Exception( "E_in = %s not in list E_ins" % EpCl.value )
                for Ep, Cl in EpCl :
                    if( Ep not in E_ins[indexE][1] ) :
                        E_ins[indexE][1][Ep] = []
                        for l in xrange( l_EEpCl.l ) : E_ins[indexE][1][Ep].append( 0. )
                    E_ins[indexE][1][Ep].append( Cl )
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 1, LEP, 1, len( E_ins ) ) ]
        ENDFDataList += endfFormats.endfInterpolationList( [ len( E_ins ), EInInterpolation ] )
        for Es in E_ins :
            NA, data = 0, []
            for key in sorted( Es[1] ) :
                LegendreSeries = Es[1][key]
                NA = max( len( LegendreSeries ) , NA )
                data += [ key ] + LegendreSeries
            ENDFDataList.append( endfFormats.endfContLine( 0, Es[0], 0, NA - 1, len( data ), len( data ) / ( NA + 1 ) ) )
            ENDFDataList += endfFormats.endfDataList( data )
        gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, tempInfo, LAW, self.productFrame, ENDFDataList )

    def toXMLList( self, indent = '' ) :

        indent2 = indent + '  '
        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for l in self : xmlString += l.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        axes_ = axes.parseXMLNode( element[0], xPath )
        frame = element.get( 'productFrame' )
        lpw = LLNLPointwise( axes_, frame )
        for lvalue in element[1:]:
            LegendreL = LLNLPointwiseEEpP( int(lvalue.get("index") ) )
            for energy_in in lvalue:
                # doesn't work yet:
                #LegendreL.EpP.append( XYs.XYs.parseXMLNode( energy_in ) )
                axes_ = axes.referenceAxes( parent=lpw )
                data = map(float, energy_in.text.split())
                data = zip( data[::2], data[1::2] )
                LegendreL.EpP.append( XYs.XYs( axes_, data, float(energy_in.get("accuracy")),
                    index=int(energy_in.get("index")), value=float(energy_in.get("value")), parent=LegendreL ) )
            lpw.append( LegendreL )
        xPath.pop()
        return lpw

    @staticmethod
    def defaultAxes( ) :

        axes_ = axes.axes( dimension = 4 )
        axes_[0] = axes.axis( 'l',          0, '',      interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[1] = axes.axis( 'energy_in',  1, 'MeV',   interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, axes.unitBaseToken ) )
        axes_[2] = axes.axis( 'energy_out', 2, 'MeV',   interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[3] = axes.axis( 'c_l',        3, '1/MeV' )
        return( axes_ )

class LLNLPointwiseEEpP :

    def __init__( self, l ) :

        self.l = l
        self.EpP = []

    def __getitem__( self, index ) :

        return( self.EpP[index] )

    def __len__( self ) :

        return( len( self.EpP ) )

    def append( self, EpP ) :

        if( not( isinstance( EpP, XYs.XYs ) ) ) : raise Exception( 'EpP is an instance of %s' % brb.getType( EpP ) )
        if( len( self ) > 0 ) :
            if( EpP.value is None ) : raise Exception( "EpP's value is None" )
            if( EpP.value <= self.EpP[-1].value ) : raise Exception( "EpP.value = %s <= prior's values = %s " % ( EpP.value, self.EpP[-1].value ) )
        self.EpP.append( EpP )

    def toXMLList( self, indent = '' ) :

        indent2 = indent + '  '
        xmlString = [ '%s<l index="%s">' % ( indent, self.l ) ]
        for energy in self : xmlString += energy.toXMLList( tag = 'energy_in', indent = indent2, pairsPerLine = 100, oneLine = True )
        xmlString[-1] += '</l>'
        return( xmlString )

class groupedBase( form ) :

    xData = matrix.matrixToken

    def __init__( self, axes_, transferMatrices, productFrame ) :

        form.__init__( self, base.groupedFormToken, productFrame )
        self.axes = axes_
        self.transferMatrices = transferMatrices

    def __len__( self ) :

        return( len( self.transferMatrices ) )

    def __getitem__( self, index ) :

        return( self.transferMatrices[index] )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        indent2 = indent + '  '
        energyEpCl = self.transferMatrices[0]
        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = 'size="%d,%d"' % ( energyEpCl.nrows, energyEpCl.ncols ) ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for l, energyEpCl in enumerate( self.transferMatrices ) :
            xmlString.append( '%s<l value="%d">' % ( indent2, l ) )
            xmlString.extend( energyEpCl.toXMLList( indent = indent2 + '  ' ) )
            xmlString[-1] += '</l>'
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseXMLNode( cls, element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        axes_ = axes.parseXMLNode( element[0], xPath )
        frame = element.get( 'productFrame' )
        transferMatrices = []
        for lvalue in element[1:]:
            transferMatrices.append( matrix.parseXMLNode( lvalue[0] ) )
        GB = cls( axes_, transferMatrices, frame )
        xPath.pop()
        return GB

class grouped( groupedBase ) :

    def __init__( self, axes_, transferMatrices, productFrame ) :

        groupedBase.__init__( self, axes_, transferMatrices, productFrame )

class energyConservationComponent( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.LegendreEnergyConservationComponentToken, nativeData )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        LC = energyConservationComponent( nativeData = element.get('nativeData') )
        for formElement in element:
            formClass = {
                    base.groupedFormToken: energyConservationGrouped,
                    }.get( formElement.tag )
            if formClass is not None:
                LC.addForm( formClass.parseXMLNode( formElement, xPath, linkData ) )
            else: raise Exception(" encountered unknown LegendreEnergyConservation form: %s" % formElement.tag)
        xPath.pop()
        return LC

class energyConservationGrouped( groupedBase ) :

    def __init__( self, axes_, transferMatrices, productFrame ) :

        groupedBase.__init__( self, axes_, transferMatrices, productFrame )
