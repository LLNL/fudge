# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Angular/energy double differential distribution classes.

This module contains the following classes:
        
    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | Form          | Class representing GNDS eneary/angular data (i.e., P(E',mu|E)).       |
    +---------------+-----------------------------------------------------------------------+
    | Subform       | Base class for :py:class:`LLNLPointwise`.                             |
    +---------------+-----------------------------------------------------------------------+
    | XYs1d         | Class used to store f(E') at a specified mu.                          |
    +---------------+-----------------------------------------------------------------------+
    | XYs2d         | Class used to store P(E',mu) at a specified E.                        |
    +---------------+-----------------------------------------------------------------------+
    | XYs3d         | Class used to store P(mu,E'|E).                                       |
    +---------------+-----------------------------------------------------------------------+
"""

import math

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.processing import group as groupModule

from . import miscellaneous as miscellaneousModule

from . import base as baseModule
from . import angular as angularModule
from . import angularEnergyMC as angularEnergyMCModule


def defaultAxes( energyUnit, probabilityLabel = 'P(mu,energy_out|energy_in)' ) :
    """
    Returns an :py:class:`Axes` instance for GNDS angular/energy data (i.e., P(mu,E'|E)).

    :param energyUnit:          Unit for the energy data.
    :param probabilityLabel:    Label for the probability axis.

    :return:            An :py:class:`Axes` instance.
    """

    axes = axesModule.Axes(4)
    axes[3] = axesModule.Axis( 'energy_in', 3, energyUnit )
    axes[2] = axesModule.Axis( 'mu', 2, '' )
    axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.Axis( probabilityLabel, 0, '1/' + energyUnit )
    return( axes )

class XYs1d( XYs1dModule.XYs1d ) :
    """
    Class to store P(E') at a specified mu. mu is specified via the outerDamainValue attribute.
    """

    def averageEnergy( self ) :
        """
        The method calculates the average outgoing energy of *self*.
        """

        return( self.integrateWithWeight_x( ) )
    
    def averageMomentum( self ) :
        """
        The method calculates the average outgoing momentum.
        """

        return( self.integrateWithWeight_sqrt_x( ) )

    def toLinearXYsClass( self ) :
        """
        This method return the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )
    
class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    Class to store P(mu,E') at a specified E. E is specified via the outerDamainValue attribute.
    """

    def  __add__( self, other ) :
        """
        This methods adds *self* to *other*. *Other* can be a float of another :py:class:`XYs2d` instance.
        """

        functionNd = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, outerDomainValue = self.outerDomainValue,
                label = self.label, interpolationQualifier = self.interpolationQualifier )

        if( isinstance( other, self.__class__ ) ) :
            x2s = list( set( [ function.outerDomainValue for function in self ] + [ function.outerDomainValue for function in other ] ) )
            for x2 in sorted( x2s ) :
                function = self.evaluate( x2 ) + other.evaluate( x2 )
                function.outerDomainValue = x2
                functionNd.append( function )
            raise NotImplementedError("Adding two angularEnergy.XYs2d instances")
        else :
            try :
                value = float( other )
            except :
                raise ValueError( 'Other must be a number or a subclass of "%s".' % other.__class__ )
            for function in self.functionals : functionNd.append( function + value )

        return( functionNd )

    __radd__ = __add__

    def averageEnergy( self ) :
        """
        The method calculates the average outgoing energy of *self*.
        """

        EpOfMu = [ [ pdfOfEpAtMu.outerDomainValue, pdfOfEpAtMu.averageEnergy( ) ] for pdfOfEpAtMu in self ]
        return XYs1dModule.XYs1d(EpOfMu).integrate()

    def averageMomentum( self ) :
        """
        The method calculates the average outgoing momentum.
        """

        MpOfMu = [ [ pdfOfMpAtMu.outerDomainValue, pdfOfMpAtMu.averageMomentum( ) ] for pdfOfMpAtMu in self ]
        return XYs1dModule.XYs1d(MpOfMu).integrate()

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs1d, ) )

class Subform( baseModule.Subform ) :
    """Abstract base class for angularEnergy forms."""

    pass

class XYs3d( Subform, multiD_XYsModule.XYs3d ) :
    """
    Class used to store the probability P(mu,E'|E).
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        Subform.__init__( self )

    def check( self, info ) :
        """
        Check for incomplete angular/energy distributions and any negative probabilities.
        Ignore normalization for this double-differential distribution.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning
        from pqu import PQU
        warnings = []

        for index, energy_in in enumerate(self):
            integral = energy_in.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.UnnormalizedDistribution( PQU.PQU( energy_in.outerDomainValue, self.axes[0].unit ), index, integral, self.toXLink() ) )
            if( energy_in.domainMin != -1 ) or ( energy_in.domainMax != 1 ) :
                warnings.append( warning.IncompleteDistribution( PQU.PQU( energy_in.outerDomainValue, self.axes[0].unit ), energy_in.domainMin, energy_in.domainMax, energy_in ) )
            for mu in energy_in:
                if mu.domainMin < 0:
                    warnings.append(warning.ValueOutOfRange("Negative outgoing energy for energy_in=%s!"
                        % PQU.PQU(energy_in.outerDomainValue, self.axes[0].unit), mu.domainMin, 0, 'inf', self.toXLink()))
                if mu.rangeMin < 0:
                    warnings.append(warning.NegativeProbability(mu.rangeMin, PQU.PQU(energy_in.outerDomainValue, self.axes[-1].unit), mu=mu.outerDomainValue, obj=mu))

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        multiplicity = kwargs['multiplicity']
        productMass = kwargs['productMass']

        aveEnergy = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.outerDomainValue
            aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageEnergy( ) ] )

        const = math.sqrt( 2. * productMass )
        aveMomentum = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.outerDomainValue
            momemtum = const * multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageMomentum( )
            if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
            aveMomentum.append( [ energy, momemtum ] )

        return( aveEnergy, aveMomentum )

    def evaluate( self, domainValue, epsilon = 0 ) :
        """
        This methods evaluates P(mu,E'|E) at the requested E = *domainValue*.

        :param domainValue:     The energy of the projectile for which P(mu,E') is requested.
        :epsilon:               If the outerDomainValue of a child 2-d function is, relatively, within this value, return that function.

        :return:                A 2-d function representing *self* at projectile energy *domainValue*.
        """

        position, function1, function2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( domainValue )
        print( type( function1 ), type( function2 ), position, frac, interpolationQualifier )
        if( position is None ) : raise Exception( "No data to interpolate" )

        if( frac <= epsilon ) :                 # If close to first point pick it.
            function = function1.copy( )
        elif( ( 1 - frac ) <= epsilon ) :       # If close to second point pick it.
            function = function2.copy( )
        else :
            if( position in ( '=', '<', '>' ) ) :
                if( position == '=' ) :
                    function = function1.copy( )
                elif( position == '>' ) :
                    function = function1.copy( )
                else :
                    index = { '<' : 0, '>' : -1 }[position]
                    raise Exception( "evaluation point = %s %s than %s" % ( domainValue, { '<' : 'less', '>' : 'greater' }[position], self[index].outerDomainValue ) )
            else :
                if( not( isinstance( function1, XYs1dModule.XYs1d ) ) ) :      # FIXME, accuracy, lowerEps and upperEps should not be hardwired.
                    function1 = function1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( not( isinstance( function2, XYs1dModule.XYs1d ) ) ) :
                    function2 = function2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )

                if( not( isinstance( function1, XYs2d ) ) ) : raise Exception( 'function1 must be an XYs2d instance' )
                if( not( isinstance( function2, XYs2d ) ) ) : raise Exception( 'function2 must be an XYs2d instance' )

                mus = sorted( set( [ function1d.outerDomainValue for function1d in function1 ] + [ function1d.outerDomainValue for function1d in function2 ] ) )
                function = XYs2d( axes = function1.axes )
                for mu in mus :
                    function1d1 = function1.evaluate( mu )
                    function1d2 = function2.evaluate( mu )
                    if( not( isinstance( function1d1, XYs1d ) ) or not( isinstance( function1d2, XYs1d ) ) ) :
                        raise Exception( 'function1d1 and function1d2 must be an XYs1d instance' )
                    xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( domainValue, function1.outerDomainValue, function1d1, function2.outerDomainValue, function1d2, 1 )
                    xy = XYs1d( xy, outerDomainValue = mu )
                    function.append( xy )

        function.outerDomainValue = domainValue
        for ii in function : print( '   ii  ', ii.outerDomainValue, ii.integrate( ) )
        function.normalize( insitu = True )
        for ii in function : print( '   ii2 ', ii.outerDomainValue, ii.integrate( ) )
        return( function )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param reaction_suite:      The :py:class:`ReactionSuite` instance for this distribution.
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.
        :param muOut:               The outgoing mu range to integrate over.
        :param phiOut:              The outgoing phi range to integrate over.
        :param frame:               The frame the outgoing energy, mu and phi range specify.
        :param LegendreOrder:       The parameter is not used.

        :return:                    A float representing the value of the integration.
        """

        if( energyIn < self.domainMin ) : return( 0.0 )
        if( energyIn > self.domainMax ) : energyIn = self.domainMax
        function2d = self.evaluate( energyIn )

        muMin, muMax = miscellaneousModule.domainLimits( muOut, function2d.domainMin, function2d.domainMax )
        if( muMax is None ) :
            function1d = function2d.evaluate( muMin )
            energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, function1d.domainMin, function1d.domainMax )
            if( energyOutMax is None ) :
                muEnergyOutEvaluate = function1d.evaluate( energyOutMin )
            else :
                muEnergyOutEvaluate = function1d.integrate( energyOutMin, energyOutMax )
        else :
            xys = []
            for function1d in function2d :
                energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, function1d.domainMin, function1d.domainMax )
                if( energyOutMax is None ) :
                    value = function1d.evaluate( energyOutMin )
                else :
                    value = function1d.integrate(energyOutMin, energyOutMax)
                xys.append( [ function1d.outerDomainValue, value ] )
            xys = XYs1d( xys )
            muEnergyOutEvaluate = xys.integrate( muMin, muMax )

        phiEvaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( phiEvaluate * muEnergyOutEvaluate )

    def normalize( self, insitu = True ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate( )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        productFrame = tempInfo['productFrame']

        if( verbosity > 2 ) : print( '%sGrouping %s' % ( indent, self.moniker ) )
        TM_1, TM_E = transferMatricesModule.ENDFEMuEpP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as an :py:class:`angularEnergyMCModule.Angular` representation
        and an :py:class:`angularEnergyMCModule.AngularEnergy` representation. This is, P(mu,E'|E) becomes
        P(mu|E) * P(E'|E,mu) with the mu data in P(mu|E) and the E' data in P(E'|E,mu) being represented with
        :py:class:`angularEnergyMCModule.Angular` and :py:class:`angularEnergyMCModule.AngularEnergy`, respectively.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`angularEnergyMCModule.Angular` and :py:class:`angularEnergyMCModule.AngularEnergy` instances.
        """

        angular = angularModule.XYs2d( axes = self.axes )
        for PofEpGivenMu in self :
            data = angularModule.XYs1d( [ [ PofEp.outerDomainValue, PofEp.integrate() ] for PofEp in PofEpGivenMu ] )
            angular.append(angularModule.Xs_pdf_cdf1d.fromXYs(angularModule.XYs1d(data), outerDomainValue=PofEpGivenMu.outerDomainValue, thinEpsilon=1e-14))

        xys3d = angularEnergyMCModule.XYs3d( axes = self.axes )
        for PofEpGivenMu in self :
            xys2d = angularEnergyMCModule.XYs2d( outerDomainValue = PofEpGivenMu.outerDomainValue )
            for PofEp in PofEpGivenMu :
                _PofEp = PofEp.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
                xys2d.append(angularEnergyMCModule.Xs_pdf_cdf1d.fromXYs(_PofEp, PofEp.outerDomainValue, thinEpsilon=1e-14))
            xys3d.append( xys2d )

        return( angularEnergyMCModule.Angular( angular ), angularEnergyMCModule.AngularEnergy( xys3d ) )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs2d, ) )

class Form(  baseModule.Form ) :
    """
    This class stores a P(mu,E'|E) distribution.
    The following table list the primary members of this class:

    +-----------------------+-----------------------------------------------------------+
    | Member                | Description                                               |
    +=======================+===========================================================+
    | angularEnergySubform  | This member stores the actual data as a 3d function.      |
    +-----------------------+-----------------------------------------------------------+
        
    :param label:                   The label for this form.
    :param productFrame:            The frame the product data are specified in.
    :param angularEnergySubform:    A :py:class:`Subform` instance.
    """

    moniker = 'angularEnergy'
    subformAttributes = ( 'angularEnergySubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularEnergySubform ) :

        if( not( isinstance( angularEnergySubform, Subform ) ) ) : raise TypeError( 'instance is not an angularEnergy subform' )
        baseModule.Form.__init__( self, label, productFrame, ( angularEnergySubform, ) )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        return( self.angularEnergySubform.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        if self.productFrame == xDataEnumsModule.Frame.centerOfMass:
            raise Exception( 'center of mass calculation not supported for %s' % self.moniker )

        kwargs['productFrame'] = self.productFrame
        aveEnergy, aveMomentum = self.angularEnergySubform.calculateAverageProductData( style, indent = indent, **kwargs )
        return( [ aveEnergy ], [ aveMomentum ] )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.angularEnergySubform.convertUnits( unitMap )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in.
        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    XYs1d instance for the energy spectrum.
        """

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)

        raise NotImplementedError('energySpectrumAtEnergy is currently not supported for an angularEnergy distribution.')

        if self.productFrame == frame:
            function2d = self.angularEnergySubform.evaluate(energyIn)
        else:
            raise TypeError('%s to %s translation not supported.' % (self.productFrame, frame))

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        This method call *fixDomains* on the *angularEnergySubform* member.
        """

        return self.angularEnergySubform.fixDomains(energyMin, energyMax, domainToFix)
        
    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param reaction_suite:      The :py:class:`ReactionSuite` instance for this distribution.
        :param energyIn:            The energy of the projectile.
        :param energyOut:           The outgoing energy range to integrate over.
        :param muOut:               The outgoing mu range to integrate over.
        :param phiOut:              The outgoing phi range to integrate over.
        :param frame:               The frame the outgoing energy, mu and phi range specify.
        :param LegendreOrder:       The parameter is not used.

        :return:                    A float representing the value of the integration.
        """

        return( self.angularEnergySubform.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an :py:class:`angularEnergyMCModule.Form` instance representing *self*.

    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        angular, angularEnergy = self.angularEnergySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( angularEnergyMCModule.Form( style.label, self.productFrame, angular, angularEnergy ) )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        tempInfo['productFrame'] = self.productFrame
        return( self.angularEnergySubform.processMultiGroup( style, tempInfo, indent ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance *cls*.

        :param cls:                 Form class to return.
        :param LegendreElement:     Node to parse.
        :param xPath:               List containing xPath to current node, useful mostly for debugging.
        :param linkData:            dict that collects unresolved links.
        :param kwargs:              A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
                     
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        formElement = element[0]
        subformClass = {
            XYs3d.moniker :  XYs3d,
        }[ formElement.tag ]
        if subformClass is None: raise Exception("encountered unknown angularEnergy subform: %s" % formElement.tag)
        subform = XYs3d.parseNodeUsingClass(formElement, xPath, linkData, **kwargs)
        AEC = cls( element.get( 'label' ), element.get('productFrame'), subform )
        xPath.pop()
        return AEC
