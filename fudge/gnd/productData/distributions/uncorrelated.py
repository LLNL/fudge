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

"""Uncorrelated double differential distribution classes"""

import math
import miscellaneous
import fudge
from fudge.core.math import fudgemath
from fudge.core.utilities import brb

import xData.standards as standardsModule
import xData.ancestry as ancestryModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

import angular as angularModule
import energy as energyModule

from . import base as baseModule

__metaclass__ = type

class subform( ancestryModule.ancestry ) :

    def __init__( self, data, dataSubform ) :

        if( not isinstance( data, dataSubform ) ) : raise TypeError( "Needed %s distribution subform" % self.moniker )
        ancestryModule.ancestry.__init__( self )
        self.data = data
        self.data.ancestor = self

    @property
    def productFrame( self ) :

        return( self.getAncestor( ).productFrame )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

class angularSubform( subform ) :

    moniker = 'angular'

    def __init__( self, data ) :

        subform.__init__( self, data, angularModule.subform )

class energySubform( subform ) :

    moniker = 'energy'

    def __init__( self, data ) :

        subform.__init__( self, data, energyModule.subform )

class form( baseModule.form ) :
    """
    Contains uncorrelated distributions for outgoing angle and outgoing energy. Use when the correlations
    between angle and energy are unknown.
    """

    moniker = 'uncorrelated'
    subformAttributes = ( 'angularSubform', 'energySubform' )

    def __init__( self, label, productFrame, _angularSubform, _energySubform, makeCopy = True ) :

        if( not isinstance( _angularSubform, angularSubform ) ) : raise TypeError( "Needed angular distribution subform, got %s" %str(_angularSubform.__class__) )
        if( not isinstance( _energySubform, energySubform ) ) : raise TypeError( "Needed energy distribution subform, got %s" % str(_energySubform.__class__) )

        baseModule.form.__init__( self, label, productFrame, ( _angularSubform, _energySubform ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.gnd.productData import multiplicity as multiplicityModule

        multiplicity = tempInfo['multiplicity']['eval']     # BRB - FIXME; Handle gamma data with 1 point for multiplicity weight
                                                    # until it is fixed.

        if( isinstance( self.angularSubform.data, angularModule.piecewise ) ) :
            raise Exception( 'piecewise angular is currently not supported' )
        if( isinstance( self.energySubform.data, energyModule.piecewise ) ) :
            energyUnit = tempInfo['incidentEnergyUnit']
            momentumDepositionUnit = energyUnit + '/c'

            axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit,
                    energyDepositionUnit = energyUnit )
            energyDeposition = energyDepositionModule.piecewise( axes = axes, label = processInfo.style.label )

            axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit,
                    momentumDepositionUnit = momentumDepositionUnit )
            momentumDeposition = momentumDepositionModule.piecewise( axes = axes, label = processInfo.style.label )

            EMin, EMax = tempInfo['EMin'], tempInfo['EMax']
            for i1, energySubform in enumerate( self.energySubform.data ) :
                tempInfo['EMax'] = energySubform.domainMax( )
                if( i1 == ( len( self.energySubform.data ) - 1 ) ) : tempInfo['EMax'] = EMax
                energyDep, momentumDep = calculateDepositionData( self.productFrame, self.angularSubform.data, energySubform,
                        processInfo, tempInfo, verbosityIndent )
                energyDeposition.append( energyDep )
                momentumDeposition.append( momentumDep )
                tempInfo['EMin'] = tempInfo['EMax']
            tempInfo['EMin'], tempInfo['EMax'] = EMin, EMax
            return( ( energyDeposition, momentumDeposition ) )

        return( calculateDepositionData( self.productFrame, self.angularSubform.data, self.energySubform.data, processInfo, tempInfo, verbosityIndent ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

        return( self.energySubform.data.getSpectrumAtEnergy( energy ) )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. uncorrelated contains both energy and angular form,
        need to check both to find desired entity
        """

        if self.energySubform.moniker == entityName: return self.energySubform
        elif self.angularSubform.moniker == entityName: return self.angularSubform
        return ancestryModule.ancestry.findEntity( self, entityName, attribute, value )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        angularSubform = self.angularSubform
        energySubform = self.energySubform
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            crossSection = tempInfo['crossSection']
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, target, product = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target, tempInfo['product']
            projectileName, productName = processInfo.getProjectileName( ), product.particle.name
            if( isinstance( energySubform, energyModule.constant ) ) :
                if( product.name == 'gamma' ) :
                    Ep = float( energySubform.data )
                    TM_1, TM_E = transferMatrices.discreteGammaAngularData( processInfo, projectileName, productName, Ep, crossSection,
                        angularSubform, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
                else :
                    raise Exception( 'See Bret' )
            else :
                if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
                    totalMass = energySubform.numberOfProductsMasses.getValueAs( massUnit )
                    Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
                    TM_1, TM_E = transferMatrices.NBodyPhaseSpace( processInfo, projectileName, productName, crossSection, \
                        energySubform.numberOfProducts, tempInfo['masses'], totalMass, Q, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
                else :
                    TM_1, TM_E = transferMatrices.uncorrelated_EMuP_EEpP_TransferMatrix( processInfo, projectileName, productName, tempInfo['masses'], \
                        crossSection, angularSubform, energySubform, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        import angularEnergy

        angularSubform = self.angularSubform.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        energySubform = self.energySubform.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        grid = angularSubform.domainGrid( )
        for e in energySubform.domainGrid( ) :
            if( e not in grid ) : grid.append( e )
        grid.sort( )
        axes = angularEnergy.pointwise.defaultAxes( energyUnit = energySubform.axes[2].unit, energy_outUnit = energySubform.axes[1].unit,
            probabilityUnit = energySubform.axes[0].unit )
        f_E_mu_Ep = angularEnergy.pointwise( axes = axes )
        for e in grid :
            f_mu_Ep = angularEnergy.pdfOfMuAndEp.pointwise( value = e )
            af = angularSubform.getAtEnergy( e )   # FIXME why doesn't getSpectrumAtEnergy work here?
            ef = energySubform.getSpectrumAtEnergy( e )
            for mu, P in af :
                efp = P * ef
                efp.__class__ = angularEnergy.pdfOfEp.pointwise # FIXME better way to assign class?
                efp.value = mu
                f_mu_Ep.append( efp )
            f_E_mu_Ep.append( f_mu_Ep )
        return( f_E_mu_Ep )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]

        xmlString += self.angularSubform.toXMLList( indent2, **kwargs )
        xmlString += self.energySubform.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker

        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <uncorrelated> element from xml. Must contain one <angular> and one <energy> component."""
        
        xPath.append( element.tag )
        #angularSubform = angularModule.component.parseXMLNode( element.find('angular'), xPath, linkData )
        #energySubform = energyModule.component.parseXMLNode( element.find('energy'), xPath, linkData )
        for child in element :
            if( child.tag == angularModule.form.moniker ) :
                _angularSubform = angularSubform( angularModule.form.parseXMLNode( child[0], xPath, linkData ) )
            elif( child.tag == energyModule.form.moniker ) :
                _energySubform = energySubform( energyModule.form.parseXMLNode( child[0], xPath, linkData ) )
            else :
                raise ValueError( 'unsupported tag = "%s"' % child.tag )
        uncorrelated = form( element.get( 'label' ), element.get( 'productFrame' ),
                _angularSubform, _energySubform, makeCopy=False )
        xPath.pop( )
        return( uncorrelated )

def calculateDepositionData( productFrame, angularSubform, energySubform, processInfo, tempInfo, verbosityIndent ) :

    def fixLimits( multiplicityLimits, Es ) :

        if( Es[0] < multiplicityLimits[0] ) :
            while( len( Es ) and ( Es[0] < multiplicityLimits[0] ) ) : del Es[0]
            if( len( Es ) and ( Es[0] > multiplicityLimits[0] ) ) : Es.insert( 0, multiplicityLimits[0] )
        if( Es[-1] > multiplicityLimits[1] ) :
            while( len( Es ) and ( Es[-1] > multiplicityLimits[1] ) ) : del Es[-1]
            if( len( Es ) and ( Es[-1] < multiplicityLimits[1] ) ) : Es.append( multiplicityLimits[1] )

    def calculateAverageEnergy( self, Ein ) :

        averageEp = self.energySubform.averageEp( Ein )
        if( self.productFrame == standardsModule.frames.centerOfMassToken ) :
            A_Ein = self.massRatio * Ein
            if( not( self.angularSubform.isIsotropic( ) ) ) : averageEp += 2 * math.sqrt( A_Ein * averageEp ) * self.angularSubform.averageMu( Ein )
            averageEp += A_Ein
        averageEp *= self.multiplicity.getValue( Ein )
        return( averageEp )

    def calculateAverageMomentum( self, Ein ) :

        pp, averageMu = 0., self.angularSubform.averageMu( E )
        if( averageMu != 0. ) : pp = energySubform.sqrtEp_AverageAtE( E ) * averageMu
        if( self.productFrame == standardsModule.frames.centerOfMassToken ) : pp += math.sqrt( self.massRatio * Ein )
        return( multiplicity.getValue( E ) * math.sqrt( 2. * self.productMass ) * pp )

    def calculateAverageMomentumForPhoton( self, Ein ) :

        muAverage = angularSubform.averageMu( E )
        if( muAverage == 0. ) : return( 0. )
        return( multiplicity.getValue( E ) * energySubform.averageEp( E ) * muAverage )

    class calculateDepositionInfo :

        def __init__( self, productFrame, productMass, massRatio, angularSubform, energySubform, multiplicity ) :

            self.productFrame = productFrame
            self.productMass = productMass
            self.massRatio = massRatio
            self.angularSubform = angularSubform
            self.energySubform = energySubform
            self.multiplicity = multiplicity

        def evaluateAtX( self, E ) :

            return( self._evaluateAtX( self, E ) )

        def setEvaluateAtX( self, evaluateAtX ) :

            self._evaluateAtX = evaluateAtX

        def setTolerances( self, relativeTolerance, absoluteTolerance ) :

            self.relativeTolerance = relativeTolerance
            self.absoluteTolerance = absoluteTolerance

    class NBodyPhaseSpace :

        def __init__( self, energySubform, massUnit, projectileMass, targetMass, productMass, Q ) :

            self.energySubform = energySubform
            self.massUnit = massUnit
            self.projectileMass = projectileMass
            self.targetMass = targetMass
            self.productMass = productMass
            self.Q = Q

        def averageEp( self, Ein ) :

            return( self.energySubform.averageEp( Ein, self.massUnit, self.projectileMass, self.targetMass, self.productMass, self.Q ) )

        def getEnergyArray( self, EMin, EMax ) :

            return( [ EMin, EMax ] )

    energyUnit = tempInfo['incidentEnergyUnit']
    momentumDepositionUnit = energyUnit + '/c'
    massUnit = energyUnit + '/c**2'
    energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy

    depData = []

    projectileMass = tempInfo['reactionSuite'].projectile.getMass( massUnit )
    targetMass = tempInfo['reactionSuite'].target.getMass( massUnit )
    productMass = tempInfo['product'].getMass( massUnit )
    massRatio = projectileMass * productMass / ( projectileMass + targetMass )**2

    if( tempInfo['product'].name == 'gamma' ) : productFrame = standardsModule.frames.labToken  # All gamma data treated as in lab frame.
    if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
        Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
        energySubform = NBodyPhaseSpace( energySubform, massUnit, projectileMass, targetMass, productMass, Q )
    multiplicity = tempInfo['multiplicity']

    calculationData = calculateDepositionInfo( productFrame, productMass, massRatio, angularSubform, energySubform, multiplicity )

    Es = energySubform.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] )
    multiplicityLimits = multiplicity.getEnergyLimits( tempInfo['EMin'], tempInfo['EMax'] )
    fixLimits( multiplicityLimits, Es )
    calculationData.setEvaluateAtX( calculateAverageEnergy )
    depEnergy = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
    absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in depEnergy ] )
    calculationData.setTolerances( energyAccuracy, absoluteTolerance )
    depEnergy = fudgemath.thickenXYList( depEnergy, calculationData )
    axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
    depData.append( energyDepositionModule.pointwise( data = depEnergy, axes = axes,
            label = processInfo.style.label, accuracy = energyAccuracy ) )

    if( isinstance( angularSubform, angularModule.isotropic ) and ( productFrame != standardsModule.frames.centerOfMassToken ) ) :
        depMomentum = [ [ depEnergy[0][0], 0. ], [ depEnergy[-1][0], 0. ] ]
    else :
        if( tempInfo['product'].name == 'gamma' ) :
            calculationData.setEvaluateAtX( calculateAverageMomentumForPhoton )
        else :
            calculationData.setEvaluateAtX( calculateAverageMomentum )
        for E in angularSubform.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] ) :
            if( E not in Es ) : Es.append( E )
        Es.sort( )
        fixLimits( multiplicityLimits, Es )
        depMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in depMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        depMomentum = fudgemath.thickenXYList( depMomentum, calculationData )
    axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
    depData.append( momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
            label = processInfo.style.label, accuracy = momentumAccuracy ) )

    return( depData )
