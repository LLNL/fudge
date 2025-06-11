# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""     
This module contains classes and functions supporting LLNL I = 1 and 3 distribution data. The distribution is
the product P(mu|E) * P(E'|E,mu) where P(mu|E) is the LLNL I = 1 data dn P(E'|E,mu) is the LLNL I = 3 data .
        
    This module contains the following classes:

    +-------------------------------------------+-----------------------------------------------------------------------+
    | Class                                     | Description                                                           |
    +===========================================+=======================================================================+
    | XYs1d                                     | Class to store the energy f(E') for a mu value.                       |
    +-------------------------------------------+-----------------------------------------------------------------------+
    | XYs2d                                     | Class to store the energy P(E'|mu) for a projectile energy.           |
    +-------------------------------------------+-----------------------------------------------------------------------+
    | XYs3d                                     | Class to store the energy P(E'|E,mu).                                 |
    +-------------------------------------------+-----------------------------------------------------------------------+
    | LLNLAngularOfAngularEnergySubform         | Class to store the angular I = 1 data.                                |
    +-------------------------------------------+-----------------------------------------------------------------------+
    | LLNLAngularEnergyOfAngularEnergySubform   | Class to store the energy I = 3 data.                                 |
    +-------------------------------------------+-----------------------------------------------------------------------+
    | LLNLAngularEnergyForm                     | Form class to store LLNL I = 1 and 3 distribution data.               |
    +-------------------------------------------+-----------------------------------------------------------------------+
"""         

import math

from PoPs import IDs as IDsPoPsModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.processing import group as groupModule

from fudge.productData import averageProductMomentum as averageProductMomentumModule

from . import miscellaneous as miscellaneousModule

from . import base as baseModule
from . import angular as angularModule
from . import energy as energyModule
from . import angularEnergyMC as angularEnergyMCModule


def defaultAxes( energyUnit, probabilityLabel = 'P(mu,energy_out|energy_in)' ) :
    """
    Returns an :py:class:`Axes` instance for LLNL I = 3 data.

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
     Class to store the energy f(E') for a mu value.
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
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return( XYs1d )
    
class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    Class to store the energy P(E'|mu) for a projectile energy.
    """

    def  __add__( self, other ) :
        """
        This method adds *other* to *self* and returns the results.

        :param other:       Number or instance of :py:class:`XYs2d` to add to *self*.

        :return:            :py:class:`XYs2d` instance.
        """

        functionNd = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, outerDomainValue = self.outerDomainValue,
                label = self.label, interpolationQualifier = self.interpolationQualifier )

        if( isinstance( other, self.__class__ ) ) :
            x2s = list( set( [ function.outerDomainValue for function in self ] + [ function.outerDomainValue for function in other ] ) )
            for x2 in sorted( x2s ) :
                function = self.evaluate( x2 ) + other.evaluate( x2 )
                function.outerDomainValue = x2
                functionNd.append( function )
            raise NotImplementedError("Adding two LLNL_angularEnergy.XYs2d instances")
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
        return XYs1dModule.XYs1d( EpOfMu ).integrate()

    def averageMomentum( self ) :
        """
        The method calculates the average outgoing momentum.
        """

        MpOfMu = [ [ pdfOfMpAtMu.outerDomainValue, pdfOfMpAtMu.averageMomentum( ) ] for pdfOfMpAtMu in self ]
        return XYs1dModule.XYs1d( MpOfMu ).integrate()

    def normalize( self, insitu = True ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        return( multiD_XYsModule.XYs2d.normalize( self, insitu = insitu, dimension = 1 ) )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( XYs1d, ) )

class XYs3d( baseModule.Subform, multiD_XYsModule.XYs3d ) :
    """
    This class stores the data for LLNL I = 3 data. That is, the P(E'|E,mu) data.
    """

    def __init__( self, **kwargs ) :

        interpolationQualifier = kwargs.get('interpolationQualifier', xDataEnumsModule.InterpolationQualifier.unitBase)
        if interpolationQualifier == xDataEnumsModule.InterpolationQualifier.none:
            interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        baseModule.Subform.__init__( self )

    def check( self, info ) :
        """
        Check for incomplete angular distributions + any negative probabilities.
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
                    warnings.append(warning.NegativeProbability(
                        mu.rangeMin, PQU.PQU(energy_in.outerDomainValue, self.axes[-1].unit), mu=mu.outerDomainValue, obj=mu))

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
        This method returns the P(E'|mu) data at the projectile energy *domainValue*.

        :param domainValue:     Projectile energy to evaluate *self* at.
        :param epsilon:         If a tabulated projectile energy is within epsilon (relatively speaking) with *domainValue*, its P(E'|mu) is returned.

        :return:                A :py:class:`XYs2d` instance.
        """

        position, function1, function2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( domainValue )
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
                if( function1[0].outerDomainValue > -1.0 ) :  # KLUDGE1: This is a kludge but there is not much one can do.
                    function = function2.copy( )
                else :
                    mus = sorted( set( [ function1d.outerDomainValue for function1d in function1 ] + [ function1d.outerDomainValue for function1d in function2 ] ) )
                    function = XYs2d( axes = function1.axes )
                    for mu in mus :
                        function1d1 = function1.evaluate( mu, interpolationQualifier = "unitbase" )
                        function1d2 = function2.evaluate( mu, interpolationQualifier = "unitbase" )
                        if( not( isinstance( function1d1, XYs1d ) ) or not( isinstance( function1d2, XYs1d ) ) ) :
                            raise Exception( 'function1d1 and function1d2 must be an XYs1d instance' )
                        xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( domainValue, function1.outerDomainValue, function1d1.nf_pointwiseXY, 
                                                                                       function2.outerDomainValue, function1d2.nf_pointwiseXY, 1 )
                        xy = xy.thinDomain(1e-6)
                        xy = XYs1d( xy, outerDomainValue = mu )
                        function.append( xy )

        function.outerDomainValue = domainValue
        function.normalize( insitu = True )
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
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs3d` instance.
        """

        return( ( XYs2d, ) )

class LLNLAngularOfAngularEnergySubform( baseModule.Subform ) :
    """
    Sub-form class to store the angular I = 1 data.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | data      | This member stores the actual angular data as a 2d function.  |
    +-----------+---------------------------------------------------------------+

    :param data:    2d function which must be a :py:class:`angularModule.XYs2` instance.
    """

    moniker = 'LLNLAngularOfAngularEnergy'

    def __init__( self, data ) :

        if( not( isinstance( data, angularModule.XYs2d ) ) ) : raise TypeError( 'instance is not an angular.XYs2d' )
        baseModule.Subform.__init__( self )
        self.data = data
        data.setAncestor(self)

    @property
    def domainMin( self ) :
        """
        Returns the minimum projectile energy of *self*'s data.
        """

        return( self.data.domainMin )

    @property
    def domainMax( self ) :
        """
        Returns the maximum projectile energy of *self*'s data.
        """

        return( self.data.domainMax )

    @property
    def domainUnit( self ) :
        """
        Returns the unit of the projectile energy of *self*'s data.
        """

        return( self.data.domainUnit )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def copy( self ) :
        """
        Returns a copy of *self*.
        """

        return( LLNLAngularOfAngularEnergySubform( self.data.copy ) )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning
        from pqu import PQU as PQUModule

        warnings = []
        for idx,function in enumerate(self.data):

            integral = function.integrate()

            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.UnnormalizedDistribution( PQUModule.PQU( function.outerDomainValue, self.data.axes[-1].unit ), idx, integral, function ) )

            if function.rangeMin < 0.0:
                warnings.append(warning.NegativeProbability(
                    function.rangeMin, PQUModule.PQU(function.outerDomainValue, self.data.axes[-1].unit), obj=function))

        return warnings

    def integrate( self, energyIn, muOut ) :
        """
        This meethod integrates the distribution at projectile energy over the specified outgoing energy, mu and phi ranges.
        See function :py:func:`miscellaneousModule.domainLimits` for how to specify *energyOut*, *muOut* and *phiOut*.
        
        :param energyIn:            The energy of the projectile.
        :param muOut:               The outgoing mu range to integrate over.

        :return:                    A float representing the value of the integration.
        """

        return( self.data.integrate( energyIn, muOut ) )

    def normalize( self, insitu = False ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        if( not( insitu ) ) : raise ValueError( 'insitu = False is currently supported.' )

        LLNLAngularOfAngularEnergy2d = self
        LLNLAngularOfAngularEnergy2d.data.normalize( insitu = True, dimension = 1 )

        return( LLNLAngularOfAngularEnergy2d )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* as an :py:class:`angularEnergyMCModule.Angular` representation.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                :py:class:`angularEnergyMCModule.Angular` instances.
        """

        return( angularEnergyMCModule.Angular( self.data.to_xs_pdf_cdf1d( style, tempInfo, indent ) ) )

    def toXML_strList( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        data = angularModule.XYs2d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        result = cls( data )
        xPath.pop( )
        return( result )

class LLNLAngularEnergyOfAngularEnergySubform( baseModule.Subform ) :
    """
    Sub-form class to store the energy I = 3 data.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | data      | This member stores the actual energy data as a 3d function.   |
    +-----------+---------------------------------------------------------------+

    :param data:    3d function for P(E'|E,mu) which must be a :py:class:`XYs3d` instance.
    """

    moniker = 'LLNLAngularEnergyOfAngularEnergy'

    def __init__( self, data ) :

        if( not( isinstance( data, XYs3d ) ) ) : raise TypeError( 'instance is not an XYs3d' )
        baseModule.Subform.__init__( self )
        self.data = data

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def copy( self ) :
        """
        Returns a copy of *self*.
        """

        return( LLNLAngularEnergyOfAngularEnergySubform( self.data.copy( ) ) )

    def check( self, info ) :
        """
        Does a check of *self*'s data.

        :param info:        A dictionary with parameters used for determining if a difference is relevant.
        """

        from fudge import warning
        from pqu import PQU as PQUModule

        warnings = []
        for idx,XYs2d in enumerate(self.data):  # looping over incident energy

            energy_in_warnings = []
            for XYs1d in XYs2d:                 # looping over mu
                integral = XYs1d.integrate()

                if abs(integral - 1.0) > info['normTolerance']:
                    energy_in_warnings.append( warning.UnnormalizedDistributionAtMu( XYs1d.outerDomainValue, integral, obj=XYs1d ) )

                if( XYs1d.rangeMin < 0.0 ) :
                    energy_in_warnings.append(warning.NegativeProbability(
                        XYs1d.rangeMin, energy_in=XYs2d.outerDomainValue, mu=XYs1d.outerDomainValue, obj=XYs1d))

            if energy_in_warnings:
                warnings.append( warning.Context("Incident energy %s (index %d)" % ( PQUModule.PQU( XYs2d.outerDomainValue, self.data.axes[-1].unit ), idx ),
                        warningList=energy_in_warnings) )

        return warnings

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

        return( self.data.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )

    def normalize( self, insitu = False ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        if( not( insitu ) ) : raise ValueError( 'insitu = False is currently supported.' )

        LLNLAngularEnergyOfAngularEnergy3d = self
        LLNLAngularEnergyOfAngularEnergy3d.data.normalize( insitu = insitu )

        return( LLNLAngularEnergyOfAngularEnergy3d )

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

        return( self.data.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

    def toXML_strList( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        data = XYs3d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        result = cls( data )
        xPath.pop( )
        return( result )

class LLNLAngularEnergyForm( baseModule.Form ) :
    """
    This class is for legacy LLNL ENDL I = 1 and 3 data and is deprecated. Only use for the legacy ENDL data.
    +-----------------------+---------------------------------------------------+
    | Member                | Description                                       |
    +=======================+===================================================+
    | angularSubform        | The angular sub-form for P(mu|E).                 |
    +-----------------------+---------------------------------------------------+
    | angularEnergySubform  | The energy sub-form for P(E'|E,mu).               |
    +-----------------------+---------------------------------------------------+

    :param label:                   The label for this form.
    :param productFrame:            The frame the product data are specified in.
    :param angularSubform:          The angular sub-form for P(mu|E) which must be a :py:class:`LLNLAngularOfAngularEnergySubform` instance.
    :param angularEnergySubform:    The energy sub-form for P(E'|E,mu) which must be a :py:class:`LLNLAngularEnergyOfAngularEnergySubform` instance.
    """

    moniker = 'LLNLAngularEnergy'
    subformAttributes = ( 'angularSubform', 'angularEnergySubform' )

    def __init__( self, label, productFrame, angularSubform, angularEnergySubform ) :

        if( not( isinstance( angularSubform, LLNLAngularOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularOfAngularEnergySubform subform' )
        if( not( isinstance( angularEnergySubform, LLNLAngularEnergyOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularEnergyOfAngularEnergySubform subform' )
        baseModule.Form.__init__( self, label, productFrame, ( angularSubform, angularEnergySubform ) )

    @property
    def domainUnit( self ) :
        """
        Returns the unit of the projectile energy of *self*'s data.
        """

        return( self.angularSubform.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        def calculateDepositionMomentumAtMu( mu, parameters ) :

            f = ( mu - parameters.mu1 ) / ( parameters.mu2 - parameters.mu1 )
            P_mu = ( 1 - f ) * parameters.P1 + f * parameters.P2
            EpP = parameters.muEpPs.interpolateAtValue(mu, unitBase=True, extrapolation=xDataEnumsModule.Extrapolation.flat)
            Ep = EpP.integrateWithWeight_sqrt_x( )
            return( mu * P_mu * Ep )

        def calculateDepositionMomentum( muPs, muEpPs ) :

            class LLNLAngular_angularEnergy_MomentumParameters :

                def __init__( self, mu1, P1, mu2, P2, muEpPs ) :

                    self.mu1, self.P1, self.mu2, self.P2, self.muEpPs = mu1, P1, mu2, P2, muEpPs

                def nextMu( self, m2, P2 ) :

                    self.mu1, self.P1 = self.mu2, self.P2
                    self.mu2, self.P2 = mu2, P2

            parameters = LLNLAngular_angularEnergy_MomentumParameters( 0, 0, 0, 0, muEpPs )
            mu1, p = None, 0
            for mu2, P2 in muPs :
                parameters.nextMu( mu2, P2 )
                if( mu1 is not None ) :
                    p += miscellaneousModule.GnG_adaptiveQuadrature( calculateDepositionMomentumAtMu, mu1, mu2, parameters, 
                            miscellaneousModule.GaussQuadrature2, tolerance = 1e-4 )[0]
                mu1 = mu2
            return( p )

        class CalculateDepositionMomentumThicken :

            def __init__( self, angularSubform, angularEnergySubform, relativeTolerance, absoluteTolerance ) :

                self.angularSubform = angularSubform
                self.angularEnergySubform = angularEnergySubform
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                muPs = self.angularSubform.getAtEnergy( E )
                muEpPs = self.angularEnergySubform.interpolateAtV( E, unitBase = True )
                return( calculateDepositionMomentum( muPs, muEpPs ) )

        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = kwargs['momentumUnit']
        multiplicity = kwargs['multiplicity']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        product = kwargs['product']
        productMass = kwargs['productMass']

        angularSubform = self.angularSubform.data
        angularEnergySubform = self.angularEnergySubform.data

        depEnergy = miscellaneousModule.calculateDepositionEnergyFromAngular_angularEnergy( style.label, 
                angularSubform, angularEnergySubform, multiplicity, accuracy = energyAccuracy )

        if( product.pid == IDsPoPsModule.photon ) :
            depMomentum = miscellaneousModule.calculateDepositionEnergyFromAngular_angularEnergy( style.label, 
                    angularSubform, angularEnergySubform, multiplicity, True, accuracy = momentumAccuracy )
        else :
            depMomentum = []
            for indexE, muPs in enumerate( angularSubform ) :
                depMomentum.append( [ muPs.outerDomainValue, calculateDepositionMomentum( muPs, angularEnergySubform[indexE] ) ] )
            const = math.sqrt( 2. * productMass )
            for EMomenutem in depMomentum : EMomenutem[1] *= const * multiplicity.evaluate( EMomenutem[0] )

        axes = averageProductMomentumModule.defaultAxes( energyUnit = energyUnit, momentumUnit = momentumUnit )
        depMomentum = averageProductMomentumModule.XYs1d( data = depMomentum, axes = axes, label = style.label )

        return( [ depEnergy ], [ depMomentum ] )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.angularSubform.convertUnits( unitMap )
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

        if frame == xDataEnumsModule.Frame.centerOfMass:
            TypeError('Lab to center-of-mass translation not supported.')

        if energyIn < self.angularSubform.data[1].outerDomainValue:         # May need kludge (see KLUDGE1).
            if self.angularSubform.data[0][0][0] != -1:
                energyIn = self.angularSubform.data[1].outerDomainValue
        angular = self.angularSubform.data.evaluate(energyIn)
        energy = self.angularEnergySubform.data.evaluate(energyIn)

        energePrimes = []
        for xys1d in energy:
            energePrimes += xys1d.domainGrid
        energePrimes = sorted(set(energePrimes))

        energePrimeProbabilities = []
        for energePrime in energePrimes:
            muProbability = []
            for muIndex, muP in enumerate(angular):
                energyAtMu = energy[muIndex]
                energyProbabilities = energyAtMu.evaluate(energePrime)
                if energyProbabilities is None:
                    continue
                muProbability.append([muP[0], muP[1] * energyProbabilities])
            energePrimeProbabilities.append([energePrime, XYs1dModule.XYs1d(muProbability).integrate(domainMin=muMin, domainMax=muMax)])
        xys1d = energyModule.XYs1d(energePrimeProbabilities, axes=energyModule.defaultAxes(energyUnit=self.domainUnit))

        return xys1d

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Calls the **fixDomains** for the **angularSubform** and **angularEnergySubform** members.
        """

        numberOfFixes  = self.angularSubform.data.fixDomains(domainMin, domainMax, fixToDomain, tweakLower = True)
        numberOfFixes += self.angularEnergySubform.data.fixDomains(domainMin, domainMax, fixToDomain, tweakLower = True)

        return numberOfFixes

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

        position, function1, function2, frac, interpolation, interpolationQualifier2 = self.angularSubform.data.getBoundingSubFunctions( energyIn )
        if( position == '' ) :          # Near threshold, some mu domains do not span [-1, 1] and the interpolation qualifier is not unit-base. This is a kludge.
            if( function1.outerDomainValue != -1.0 ) : energyIn = function2.outerDomainValue
        angularEvaluate = self.angularSubform.integrate( energyIn, muOut )
        energyPartialIntegral = self.angularEnergySubform.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, 
                LegendreOrder = LegendreOrder )    # Does the phiOut integration.

        return( angularEvaluate * energyPartialIntegral )

    def normalize( self, insitu = False ) :
        """
        This methods returns a renormalized copy (or *self) of *self*.

        :param insitu:      If True, normalizes *self*, else copies *self* and return the normalized copy.

        :return:            Returns the normalzied function which can be *self* if *insitu* is True.
        """

        if( not( insitu ) ) : raise ValueError( 'insitu = False is currently supported.' )

        LLNLAngularEnergy3d = self

        LLNLAngularEnergy3d.angularSubform.normalize( insitu = True )
        LLNLAngularEnergy3d.angularEnergySubform.normalize( insitu = True )

        return( LLNLAngularEnergy3d )

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an :py:class:`angularEnergyMCModule.Form` instance representing *self*.

    
        :param style            The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        angular = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        _, angularEnergy = self.angularEnergySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( angularEnergyMCModule.Form( style.label, self.productFrame, angular, angularEnergy ) )

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
        if( verbosity > 2 ) : print( '%sGrouping %s' % ( indent, self.moniker ) )

        angularSubform = self.angularSubform.data
        angularEnergySubform = self.angularEnergySubform.data

        TM_1, TM_E = transferMatricesModule.ENDLEMuEpP_TransferMatrix( style, tempInfo, tempInfo['crossSection'], self.productFrame, 
            angularSubform, angularEnergySubform, tempInfo['multiplicity'],
            comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toXML_strList( self, indent = "", **kwargs ) : 
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]

        xmlString += self.angularSubform.toXML_strList( indent3, **kwargs )

        xmlString += self.angularEnergySubform.toXML_strList( indent3, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker 
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )
        for child in element :
            if( child.tag == LLNLAngularOfAngularEnergySubform.moniker ) :
                angularSubform = LLNLAngularOfAngularEnergySubform.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == LLNLAngularEnergyOfAngularEnergySubform.moniker ) :
                angularEnergySubform = LLNLAngularEnergyOfAngularEnergySubform.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise ValueError( 'Encountered unexpected child element "%s"' % child.tag )
        component = cls( element.get( 'label' ), element.get( 'productFrame' ), 
                angularSubform, angularEnergySubform )

        xPath.pop( )
        return( component )
