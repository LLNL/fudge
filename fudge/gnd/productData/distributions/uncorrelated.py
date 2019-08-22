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

"""Uncorrelated double differential distribution classes."""

import math
import miscellaneous
import fudge
from fudge.core.math import fudgemath
from fudge.core.utilities import brb

from PoPs import IDs as IDsPoPsModule

import xData.standards as standardsModule
import xData.ancestry as ancestryModule

import angular as angularModule
import energy as energyModule

from . import base as baseModule

__metaclass__ = type

class subform( ancestryModule.ancestry ) :

    ancestryMembers = ( 'data', )

    def __init__( self, data, dataSubform ) :

        if( not isinstance( data, dataSubform ) ) : raise TypeError( "Needed %s distribution subform" % self.moniker )
        ancestryModule.ancestry.__init__( self )
        self.data = data
        self.data.setAncestor( self )

    @property
    def productFrame( self ) :

        return( self.getAncestor( ).productFrame )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def copy( self ) :

        return( self.__class__( self.data.copy( ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        data = self.data.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( data is None ) : return( None )
        return( self.__class__( data ) )

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
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angularSubform, _energySubform ) :

        if( not isinstance( _angularSubform, angularSubform ) ) : raise TypeError( "Needed angular distribution subform, got %s" %str(_angularSubform.__class__) )
        if( not isinstance( _energySubform, energySubform ) ) : raise TypeError( "Needed energy distribution subform, got %s" % str(_energySubform.__class__) )

        baseModule.form.__init__( self, label, productFrame, ( _angularSubform, _energySubform ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        if( isinstance( self.angularSubform.data, angularModule.regions2d ) ) :
            raise Exception( 'regions2d angular is currently not supported' )

        if( isinstance( self.energySubform.data, energyModule.regions2d ) ) :
            aveEnergies = []
            aveMomenta = []
            EMax = kwargs['EMax']
            for i1, energySubform in enumerate( self.energySubform.data ) :
                kwargs['EMax'] = energySubform.domainMax
                if( i1 == ( len( self.energySubform.data ) - 1 ) ) : kwargs['EMax'] = EMax
                aveEnergy, aveMomentum = calculateAverageProductData( self.productFrame, self.angularSubform.data, energySubform,
                        style, indent, **kwargs )
                aveEnergies.append( aveEnergy )
                aveMomenta.append( aveMomentum )
                kwargs['EMin'] = kwargs['EMax']
            return( aveEnergies, aveMomenta )

        aveEnergy, aveMomentum = calculateAverageProductData( self.productFrame, self.angularSubform.data, self.energySubform.data, style, indent, **kwargs )
        return( [ aveEnergy ], [ aveMomentum ] )

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

    def processMC( self, style, tempInfo, indent ) :

        _angularSubform = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        _energySubform = self.energySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( ( _angularSubform is not None ) or ( _energySubform is not None ) ) :
            if( _angularSubform is None ) : _angularSubform = self.angularSubform.copy( )
            if( _energySubform is None ) : _energySubform = self.energySubform.copy( )
            _form = form( style.label, self.productFrame, _angularSubform, _energySubform )
        else :
            _form = None
        return( _form )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']

        angularSubform = self.angularSubform.data
        energySubform = self.energySubform.data
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )

        crossSection = tempInfo['crossSection']
        product = tempInfo['product']
        if( isinstance( energySubform, energyModule.discreteGamma ) ) :
            if( product.id == IDsPoPsModule.photon ) :
                Ep = float( energySubform.value )
                TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( style, tempInfo, Ep, crossSection,
                        angularSubform, multiplicity = product.multiplicity.evaluated,
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            else :
                raise Exception( 'See Bret' )
        else :
            if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
                totalMass = energySubform.numberOfProductsMasses.getValueAs( massUnit )
                Q = tempInfo['reaction'].getQ( energyUnit, final = False, groundStateQ = True )
                TM_1, TM_E = transferMatricesModule.NBodyPhaseSpace( style, tempInfo, crossSection, 
                        energySubform.numberOfProducts, totalMass, Q, tempInfo['multiplicity'], 
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            else :
                TM_1, TM_E = transferMatricesModule.uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, 
                        self.productFrame, angularSubform, energySubform, tempInfo['multiplicity'], 
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        import angularEnergy

        angularSubform = self.angularSubform.toPointwise_withLinearXYs( **kwargs )
        energySubform = self.energySubform.toPointwise_withLinearXYs( **kwargs )
        grid = angularSubform.domainGrid
        for e in energySubform.domainGrid :
            if( e not in grid ) : grid.append( e )
        grid.sort( )
        axes = angularEnergy.defaultAxes( energyUnit = energySubform.axes[-1].unit )
        f_E_mu_Ep = angularEnergy.XYs3d( axes = axes )
        for e in grid :
            f_mu_Ep = angularEnergy.XYs2d( value = e )
            af = angularSubform.getAtEnergy( e )   # FIXME why doesn't getSpectrumAtEnergy work here?
            ef = energySubform.getSpectrumAtEnergy( e )
            for mu, P in af :
                efp = P * ef
                efp.__class__ = angularEnergy.XYs1d # FIXME better way to assign class?
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
        uncorrelated = form( element.get( 'label' ), element.get( 'productFrame' ), _angularSubform, _energySubform )
        xPath.pop( )
        return( uncorrelated )

def calculateAverageProductData( productFrame, angularSubform, energySubform, style, indent, **kwargs ) :

    def fixLimits( EMin, EMax, Es ) :

        if( EMin not in Es ) : Es.append( EMin )
        if( EMax not in Es ) : Es.append( EMax )
        Es.sort( )
        while( Es[0] < EMin ) : del Es[0]
        while( Es[-1] > EMax ) : del Es[-1]

    def calculateAverageEnergy( self, Ein ) :

        averageEp = self.energySubform.averageEp( Ein )
        if( self.productFrame == standardsModule.frames.centerOfMassToken ) :
            A_Ein = self.massRatio * Ein
            if( not( self.angularSubform.isIsotropic( ) ) ) : averageEp += 2 * math.sqrt( A_Ein * averageEp ) * self.angularSubform.averageMu( Ein )
            averageEp += A_Ein
        averageEp *= self.multiplicity.evaluate( Ein )
        return( averageEp )

    def calculateAverageMomentum( self, Ein ) :

        pp, averageMu = 0., self.angularSubform.averageMu( E )
        if( averageMu != 0. ) : pp = energySubform.sqrtEp_AverageAtE( E ) * averageMu
        if( self.productFrame == standardsModule.frames.centerOfMassToken ) : pp += math.sqrt( self.massRatio * Ein )
        return( multiplicity.evaluate( E ) * math.sqrt( 2. * self.productMass ) * pp )

    def calculateAverageMomentumForPhoton( self, Ein ) :

        muAverage = angularSubform.averageMu( E )
        if( muAverage == 0. ) : return( 0. )
        return( multiplicity.evaluate( E ) * energySubform.averageEp( E ) * muAverage )

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

        def getEnergyArray( self, EMin = None, EMax = None ) :

            return( [ EMin, EMax ] )

    energyUnit = kwargs['incidentEnergyUnit']
    massUnit = kwargs['massUnit']
    multiplicity = kwargs['multiplicity']               # BRB - FIXME; Handle gamma data with 1 point for multiplicity weight until it is fixed.
    energyAccuracy = kwargs['energyAccuracy']
    momentumAccuracy = kwargs['momentumAccuracy']
    projectileMass = kwargs['projectileMass']
    targetMass = kwargs['targetMass']
    product = kwargs['product']
    productMass = kwargs['productMass']

    EMin = kwargs['EMin']
    EMax = kwargs['EMax']

    massRatio = projectileMass * productMass / ( projectileMass + targetMass )**2

    if( product.id == IDsPoPsModule.photon ) : productFrame = standardsModule.frames.labToken  # All gamma data treated as in lab frame.
    if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
        Q = kwargs['reaction'].getQ( energyUnit, final = False, groundStateQ = True )
        energySubform = NBodyPhaseSpace( energySubform, massUnit, projectileMass, targetMass, productMass, Q )

    calculationData = calculateDepositionInfo( productFrame, productMass, massRatio, angularSubform, energySubform, multiplicity )

    Es = energySubform.getEnergyArray( )
    if( Es[0] is None ) : Es[0] = EMin
    if( Es[-1] is None ) : Es[-1] = EMax
    if( EMin < Es[0] ) : EMin = Es[0]
    if( EMax > Es[-1] ) : EMax = Es[-1]
    Es = sorted( set( energySubform.getEnergyArray( EMin, EMax ) + multiplicity.domainGrid ) )
    fixLimits( EMin, EMax, Es )
    calculationData.setEvaluateAtX( calculateAverageEnergy )
    aveEnergy = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
    absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in aveEnergy ] )
    calculationData.setTolerances( energyAccuracy, absoluteTolerance )
    aveEnergy = fudgemath.thickenXYList( aveEnergy, calculationData )

    if( isinstance( angularSubform, angularModule.isotropic ) and ( productFrame == standardsModule.frames.labToken ) ) :
        aveMomentum = [ [ aveEnergy[0][0], 0. ], [ aveEnergy[-1][0], 0. ] ]
    else :
        if( product.id == IDsPoPsModule.photon ) :
            calculationData.setEvaluateAtX( calculateAverageMomentumForPhoton )
        else :
            calculationData.setEvaluateAtX( calculateAverageMomentum )
        Es = sorted( set( Es + angularSubform.getEnergyArray( EMin, EMax ) ) )
        fixLimits( EMin, EMax, Es )
        aveMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in aveMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        aveMomentum = fudgemath.thickenXYList( aveMomentum, calculationData )

    return( aveEnergy, aveMomentum )
