# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""             
This module contains classes and functions supporting uncorrelated double differential distribution.
                
    This module contains the following classes: 
            
    +-------------------+-----------------------------------------------------------------------+
    | Class             | Description                                                           |
    +===================+=======================================================================+
    | Subform           | Base class for :py:class:`AngularSubform` and                         |
    |                   | :py:class:`EnergySubform` classes.                                    |
    +-------------------+-----------------------------------------------------------------------+
    | AngularSubform    | This class represents the angular part of the distribution.           |
    +-------------------+-----------------------------------------------------------------------+
    | EnergySubform     | This class represents the energy part of the distribution.            |
    +-------------------+-----------------------------------------------------------------------+
    | Form              | Class representing the **GNDS** uncorrelated distribution. This       |
    |                   | distribution has 2 child nodes. One representing the angular part     |
    |                   | and one representing the energy part of the distribution.             |
    +-------------------+-----------------------------------------------------------------------+
"""

import math
from fudge.core.math import fudgemath

from PoPs import IDs as IDsPoPsModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import values as valuesModule

from . import angular as angularModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import angularEnergy as angularEnergyModule

from . import base as baseModule
from . import miscellaneous as miscellaneousModule


class Subform( ancestryModule.AncestryIO ) :
    """
    Base class for the :py:class:`AngularSubform` and :py:class:`EnergySubform` classes.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the actual data as a 2d function.      |
    +-----------+-----------------------------------------------------------+

    :param data:            2d function representing the angular or energy data.
    :param dataSubform:     Allowed class for *data*.
    """

    ancestryMembers = ( 'data', )

    def __init__( self, data, dataSubform ) :

        ancestryModule.AncestryIO.__init__( self )

        if( not isinstance( data, dataSubform ) ) : raise TypeError( "Needed %s distribution subform" % self.moniker )
        self.data = data
        self.data.setAncestor( self )

    @property
    def productFrame( self ) :
        """Returns the frame the product data are specified in."""

        return( self.ancestor.productFrame )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def copy( self ) :
        """
        Returns a copy of *self*.

        :return:        An instance that is the same class as *self*.
        """

        return self.__class__(self.data.copy())

    def fixDomains(self, energyMin, energyMax, domainsToFix):
        """
        Calls the *fixDomains* method on the *data* member.
        """

        return self.data.fixDomains(energyMin, energyMax, domainsToFix)

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :
        """
        This method returns a copy of *self* in a representation suitable for Monte Carlo transport.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
        """

        data = self.data.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( data is None ) : return( None )
        return self.__class__(data)

    def findEntity(self, entityName, attribute=None, value=None):
        """
        Overrides :py:meth:`ancestryModule.Ancestry.findEntity`.
        """

        if self.data.moniker == entityName:
            return self.data
        return ancestryModule.Ancestry.findEntity(self, entityName, attribute, value)

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        This method does nothing.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        pass

class AngularSubform( Subform ) :
    """
    Class representing the angular data for an uncorrelated distribution.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the actual data as a 2d function.      |
    +-----------+-----------------------------------------------------------+
    
    :param data:            2d function representing the angular data.
    """

    moniker = 'angular'

    def __init__( self, data ) :

        Subform.__init__( self, data, angularModule.Subform )

class EnergySubform( Subform ) :
    """
    Class representing the energy data for an uncorrelated distribution.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the actual data as a 2d function.      |
    +-----------+-----------------------------------------------------------+
    
    :param data:            2d function representing the angular data.
    """

    moniker = 'energy'

    def __init__( self, data ) :

        Subform.__init__( self, data, energyModule.Subform )

class Form( baseModule.Form ) :
    """
    Contains uncorrelated distributions for outgoing angle and outgoing energy. Use when the correlations
    between angle and energy are unknown.

    The following table list the primary members of this class:

    +-------------------+-----------------------------------------------------------+
    | Member            | Description                                               |
    +===================+===========================================================+
    | angularSubform    | Angular data.                                             |
    +-------------------+-----------------------------------------------------------+
    | energySubform     | Energy data.                                              |
    +-------------------+-----------------------------------------------------------+
    
    :param label:               The label for this form.
    :param productFrame:        The frame the product data are specified in.
    :param _angularSubform:     An :py:class:`AngularSubform` representing the energy data.
    :param _energySubform:      An :py:class:`EnergySubform` representing the energy data.
    """

    moniker = 'uncorrelated'
    subformAttributes = ( 'angularSubform', 'energySubform' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angularSubform, _energySubform ) :

        if( not isinstance( _angularSubform, AngularSubform ) ) : raise TypeError( "Needed angular distribution subform, got %s" %str(_angularSubform.__class__) )
        if( not isinstance( _energySubform, EnergySubform ) ) : raise TypeError( "Needed energy distribution subform, got %s" % str(_energySubform.__class__) )

        baseModule.Form.__init__( self, label, productFrame, ( _angularSubform, _energySubform ) )

    @property
    def domainUnit( self ) :
        """
        Returns the energy unit of the projectile.
        """

        return( self.angularSubform.data.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        if( isinstance( self.angularSubform.data, angularModule.Regions2d ) ) :
            raise Exception( 'Regions2d angular is currently not supported' )

        if( isinstance( self.energySubform.data, energyModule.Regions2d ) ) :
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

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """
        Calculates the outgoing particle's energy spectrum at projectile energy *energyIn* for frame *frame*,
        Currently, if frame='lab' and the product is a photon then its data are treated as being in the lab frame.
        The domain of angular integration can be set with the muMin and muMax keys in *kwargs*. Default angular integration
        is from -1 to 1.

        :param energy_in:           Energy of the projectile.
        :param frame:               The frame to calculate the energy spectrum in. 
        :param kwargs:              A dictionary that contains data to control the way this method acts.
        
        :return:                    XYs1d instance for the energy spectrum.
        """

        class AngualarAtEnergyCOM:

            def __init__(self, angular):

                self.angular = angular

            def probabilityCOM(self, energyPrimeCOM, muCOM):

                return self.angular.evaluate(muCOM)

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)

        if self.productFrame == frame:
            energySpectrum = self.energySubform.data.energySpectrumAtEnergy(energyIn)
            if muMin != -1.0 or muMax != 1.0:
                energySpectrum *= self.angularSubform.data.evaluate(energyIn).integrate(domainMin=muMin, domainMax=muMax)
            return energySpectrum

        if self.productFrame == xDataEnumsModule.Frame.lab:
            TypeError('Lab to center-of-mass translation not supported.')

        if self.product.pid != IDsPoPsModule.photon:
            energyProbabilityAtEnergy = self.energySubform.data.evaluate(energyIn)
            angularProbabilityAtEnergy = self.angularSubform.data.evaluate(energyIn)
            xys2d = miscellaneousModule.energyAngularSpectrumFromCOMSpectrumToLabAtEnergy(self, energyIn, energyProbabilityAtEnergy, 
                    AngualarAtEnergyCOM(angularProbabilityAtEnergy))
            data = [[xys1d.outerDomainValue, xys1d.integrate(domainMin=muMin, domainMax=muMax)] for xys1d in xys2d]
            return energyModule.XYs1d(data=data, axes=energyModule.defaultAxes(self.domainUnit))
        else:                          # Ignore centerOfMass corrections for photons, for now anyway.
            energySpectrum = self.energySubform.data.energySpectrumAtEnergy(energyIn)
            if muMin != -1.0 or muMax != 1.0:
                energySpectrum *= self.angularSubform.data.evaluate(energyIn).integrate(domainMin=muMin, domainMax=muMax)
            return energySpectrum

    def getSpectrumAtEnergy(self, energy):
        """
        This method is deprecated, please use energySpectrumAtEnergy instead.
        Returns the energy spectrum for self at projectile energy in the lab frame.
        """

        return self.energySpectrumAtEnergy(energy, xDataEnumsModule.Frame.lab)

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides :py:func:`ancestryModule.Ancestry.findEntity`.
        Uncorrelated contains both energy and angular form,
        need to check both to find desired entity
        """

        if self.energySubform.moniker == entityName: return self.energySubform
        elif self.angularSubform.moniker == entityName: return self.angularSubform
        return ancestryModule.Ancestry.findEntity( self, entityName, attribute, value )

    def fixDomains(self, energyMin, energyMax, domainsToFix):
        """
        Calls the *fixDomains* method on the *energy* and *angular* members.
        """

        numberOfFixes  = self.energySubform.fixDomains(energyMin, energyMax, domainsToFix)
        numberOfFixes += self.angularSubform.fixDomains(energyMin, energyMax, domainsToFix)

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

        angularEvaluate = 1.0
        angularData = self.angularSubform.data
        if( muOut is None ) :
            angularEvaluate = angularData.evaluate( energyIn ).LegendreCoefficient( LegendreOrder )
        else :
            angular1d = angularData.evaluate( energyIn )
            muMin, muMax = miscellaneousModule.domainLimits( muOut, angular1d.domainMin, angular1d.domainMax )
            try :
                if( muMax is None ) :
                    angularEvaluate = angular1d.evaluate( muMin )
                else :
                    angularEvaluate = angular1d.integrate( muMin, muMax )
            except :
                print( type( angularData ) )
                raise

        energyPartialIntegral = 0.0
        energyData = self.energySubform.data
        if( isinstance( energyData, ( energyModule.DiscreteGamma, energyModule.PrimaryGamma ) ) ) :
            energyPartialIntegral = energyData.integrate( energyIn, energyOut )
        else :
            energyData1d = energyData.evaluate( energyIn )
            energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, energyData1d.domainMin, energyData1d.domainMax )
            if( energyOutMax is None ) :
                energyPartialIntegral = energyData1d.evaluate( energyOutMin )
            else :
                if( isinstance( energyData1d, energyModule.Regions1d ) ) :
                    limits = { energyData1d.axes[1].label : ( energyOutMin, energyOutMax ) }
                    energyPartialIntegral = energyData1d.integrate(**limits)
                else :
                    energyPartialIntegral = energyData1d.integrate(energyOutMin, energyOutMax)

        phiEvaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( phiEvaluate * angularEvaluate * energyPartialIntegral )

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns an :py:class:`Form` instance representing *self* with (xs, pdf, cdf) data 
        for P(E') and P(mu).
    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                An instance of self.
        """     

        _angularSubform = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        _energySubform = self.energySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        if( ( _angularSubform is not None ) or ( _energySubform is not None ) ) :
            if( _angularSubform is None ) : _angularSubform = self.angularSubform.copy( )
            if( _energySubform is None ) : _energySubform = self.energySubform.copy( )
            _form = Form( style.label, self.productFrame, _angularSubform, _energySubform )
        else :
            _form = None
        return( _form )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']

        angularSubform = self.angularSubform.data
        energySubform = self.energySubform.data
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

        crossSection = tempInfo['crossSection']
        product = tempInfo['product']
        if( isinstance( energySubform, energyModule.DiscreteGamma ) ) :
            if( product.pid == IDsPoPsModule.photon ) :
                Ep = float( energySubform.value )
                TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( style, tempInfo, Ep, crossSection,
                        angularSubform, multiplicity = product.multiplicity.evaluated,
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            else :
                raise Exception( 'See Bret' )
        else :
            if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
                totalMass = energySubform.mass.getValueAs( massUnit )
                Q = tempInfo['reaction'].getQ( energyUnit, final = False )
                TM_1, TM_E = transferMatricesModule.NBodyPhaseSpace( style, tempInfo, crossSection, 
                        energySubform.numberOfProducts, totalMass, Q, tempInfo['multiplicity'], 
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            else :
                TM_1, TM_E = transferMatricesModule.uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, 
                        self.productFrame, angularSubform, energySubform, tempInfo['multiplicity'], 
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def spectrum( self, frame, **kwargs ) :
        """
        Returns an :py:class:`energyAngularModule.XYs3d` instance of *self*.

        :param frame:       The desired frame the distribution is to represent.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted
        """

        if( frame == self.productFrame ) : return( self.toPointwise_withLinearXYs( **kwargs ) )

        if frame == xDataEnumsModule.Frame.centerOfMass:
            raise Exception( 'Conversion of distribution from %s to %s not supported' % ( self.productFrame, frame ) )

        energiesIn = self.angularSubform.data.domainGrid + self.energySubform.data.domainGrid
        energiesIn = valuesModule.Values.sortAndThin( energiesIn, rel_tol = 1e-12 )
        xys3d = energyAngularModule.XYs3d( axes = energyAngularModule.defaultAxes( self.domainUnit ) )
        for energy in energiesIn : xys3d.append( self.spectrumAtEnergy( energy, frame ) )

        return( xys3d )

    def spectrumAtEnergy( self, energyIn, frame ) :
        """
        Returns a XYs2d instance of *self* evaluated at *energyIn*.

        :param energyIn:    The energy of the projectile that the distribution is evaluated at.
        :param frame:       The desired frame the distribution is to represent.
        """

        class AngualarAtEnergyCOM :

            def __init__( self, angular ) :

                self.angular = angular

            def probabilityCOM( self, energyPrimeCOM, muCOM ) :

                return( self.angular.evaluate( muCOM ) )

        angularSubform = self.angularSubform.data.toPointwise_withLinearXYs( )
        angularProbability = angularSubform.getAtEnergy( energyIn )   # FIXME why doesn't getSpectrumAtEnergy work here?

        energySubform = self.energySubform.data.toPointwise_withLinearXYs( )
        energyProbability = energySubform.getSpectrumAtEnergy( energyIn )

        if( frame == self.productFrame ) :
            probability = energyAngularModule.XYs2d( outerDomainValue = energyIn )
            for energyOut, P in energyProbability :
                efp = energyAngularModule.XYs1d( P * angularProbability, outerDomainValue = energyOut )
                probability.append( efp )

            return( probability )

        if frame == xDataEnumsModule.Frame.centerOfMass:
            raise Exception( 'Conversion of distribution from %s to %s not supported' % ( self.productFrame, frame ) )

        return( miscellaneousModule.energyAngularSpectrumFromCOMSpectrumToLabAtEnergy( self, energyIn, energyProbability, AngualarAtEnergyCOM( angularProbability ) ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    An :py:class:`angularEnergyModule.XYs3d` instance.
        """

        angularSubform = self.angularSubform.data.toPointwise_withLinearXYs( **kwargs )
        energySubform = self.energySubform.data.toPointwise_withLinearXYs( **kwargs )
        grid = angularSubform.domainGrid
        for e in energySubform.domainGrid :
            if( e not in grid ) : grid.append( e )
        grid.sort( )
        axes = angularEnergyModule.defaultAxes( energyUnit = energySubform.axes[-1].unit )
        f_E_mu_Ep = angularEnergyModule.XYs3d( axes = axes )
        for e in grid :
            f_mu_Ep = angularEnergyModule.XYs2d( outerDomainValue = e )
            af = angularSubform.getAtEnergy( e )   # FIXME why doesn't getSpectrumAtEnergy work here?
            ef = energySubform.getSpectrumAtEnergy( e )
            for mu, P in af :
                efp = P * ef
                efp.__class__ = angularEnergyModule.XYs1d # FIXME better way to assign class?
                efp.outerDomainValue = mu
                f_mu_Ep.append( efp )
            f_E_mu_Ep.append( f_mu_Ep )
        return( f_E_mu_Ep )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]

        xmlString += self.angularSubform.toXML_strList( indent2, **kwargs )
        xmlString += self.energySubform.toXML_strList( indent2, **kwargs )
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
            if( child.tag == angularModule.Form.moniker ) :
                _angularSubform = AngularSubform(angularModule.Form.parseNodeUsingClass(child[0], xPath, linkData, **kwargs))
            elif( child.tag == energyModule.Form.moniker ) :
                _energySubform = EnergySubform(energyModule.Form.parseNodeUsingClass(child[0], xPath, linkData, **kwargs))
            else :
                raise ValueError( 'unsupported tag = "%s"' % child.tag )
        uncorrelated = cls( element.get( 'label' ), element.get( 'productFrame' ), _angularSubform, _energySubform )
        xPath.pop( )
        return( uncorrelated )

def calculateAverageProductData( productFrame, angularSubform, energySubform, style, indent, **kwargs ) :
    """         
    This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
    Average means over all outgoing angle and energy.

    :param productFrame:    The frame the product data are specified in.
    :param angularSubform:  2d function of the angular data.
    :param energySubform:   2d function of the energy data.
    :param style:           The style instance which the calculated data will belong to.
    :param indent:          If this method does any printing, this is the amount of indentation of the printed line.
    :param kwargs:          A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
    """

    def fixLimits( EMin, EMax, Es ) :

        if( EMin not in Es ) : Es.append( EMin )
        if( EMax not in Es ) : Es.append( EMax )
        Es.sort( )
        while( Es[0] < EMin ) : del Es[0]
        while( Es[-1] > EMax ) : del Es[-1]
        return( valuesModule.Values.sortAndThin( Es, 1e-6 ) )

    def calculateAverageEnergy( self, Ein ) :

        averageEp = self.energySubform.averageEp( Ein )
        if self.productFrame == xDataEnumsModule.Frame.centerOfMass:
            A_Ein = self.massRatio * Ein
            if( not( self.angularSubform.isIsotropic( ) ) ) : averageEp += 2 * math.sqrt( A_Ein * averageEp ) * self.angularSubform.averageMu( Ein )
            averageEp += A_Ein
        averageEp *= self.multiplicity.evaluate( Ein )
        return( averageEp )

    def calculateAverageMomentum( self, Ein ) :

        pp, averageMu = 0., self.angularSubform.averageMu( Ein )
        if( averageMu != 0. ) : pp = energySubform.sqrtEp_AverageAtE( Ein ) * averageMu
        if self.productFrame == xDataEnumsModule.Frame.centerOfMass:
            pp += math.sqrt( self.massRatio * Ein )
        return( multiplicity.evaluate( Ein ) * math.sqrt( 2. * self.productMass ) * pp )

    def calculateAverageMomentumForPhoton( self, Ein ) :

        muAverage = angularSubform.averageMu( Ein )
        if( muAverage == 0. ) : return( 0. )
        return( multiplicity.evaluate( Ein ) * energySubform.averageEp( Ein ) * muAverage )

    class CalculateDepositionInfo :

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

    if( product.pid == IDsPoPsModule.photon ) : productFrame = xDataEnumsModule.Frame.lab    # All gamma data treated as in lab frame.
    if( isinstance( energySubform, energyModule.NBodyPhaseSpace ) ) :
        Q = kwargs['reaction'].getQ( energyUnit, final = False )
        energySubform = NBodyPhaseSpace( energySubform, massUnit, projectileMass, targetMass, productMass, Q )

    if kwargs['isInfiniteTargetMass']:
        aveEnergy = [ [ functional1d.outerDomainValue, functional1d.integrateWithWeight_x( ) ] for functional1d in energySubform ]
    else:
        Es = energySubform.getEnergyArray()
        if Es[0] is None: Es[0] = EMin
        if Es[-1] is None: Es[-1] = EMax
        if EMin < Es[0]: EMin = Es[0]
        if EMax > Es[-1]: EMax = Es[-1]
        Es = sorted(set(Es + multiplicity.domainGrid))
        Es = fixLimits(EMin, EMax, Es)
        massRatio = projectileMass * productMass / (projectileMass + targetMass)**2
        calculationData = CalculateDepositionInfo(productFrame, productMass, massRatio, angularSubform, energySubform, multiplicity)
        calculationData.setEvaluateAtX(calculateAverageEnergy)
        aveEnergy = [[E, calculationData.evaluateAtX(E)] for E in Es]
        absoluteTolerance = 1e-3 * energyAccuracy * max([Ep for E, Ep in aveEnergy])
        calculationData.setTolerances(energyAccuracy, absoluteTolerance)
        aveEnergy = fudgemath.thickenXYList(aveEnergy, calculationData)

    if isinstance(angularSubform, angularModule.Isotropic2d) and productFrame == xDataEnumsModule.Frame.lab:
        aveMomentum = [ [ aveEnergy[0][0], 0. ], [ aveEnergy[-1][0], 0. ] ]
    else :
        if( product.pid == IDsPoPsModule.photon ) :
            calculationData.setEvaluateAtX( calculateAverageMomentumForPhoton )
        else :
            calculationData.setEvaluateAtX( calculateAverageMomentum )
        Es = sorted( set( Es + angularSubform.getEnergyArray( EMin, EMax ) ) )
        Es = fixLimits( EMin, EMax, Es )
        aveMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in aveMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        aveMomentum = fudgemath.thickenXYList( aveMomentum, calculationData )

    return( aveEnergy, aveMomentum )
