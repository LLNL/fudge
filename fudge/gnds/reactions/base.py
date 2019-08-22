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

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule

import xData.ancestry as ancestryModule

import fudge
from fudge.core.utilities import fudgeExceptions

from .. import channels as channelsModule
from ..reactionData.doubleDifferentialCrossSection import doubleDifferentialCrossSection as doubleDifferentialCrossSectionModule
from ..reactionData import crossSection as crossSectionModule
from ..reactionData import availableEnergy as availableEnergyModule
from ..reactionData import availableMomentum as availableMomentumModule

__metaclass__ = type

class base_reaction( ancestryModule.ancestry ) :
    """Base class for all types of reaction."""

    ancestryMembers = ( 'crossSection', 'doubleDifferentialCrossSection', 'outputChannel' )

    def __init__( self, outputChannel, ENDF_MT, documentation = None, label = None, process = None ) :

        ancestryModule.ancestry.__init__( self )
        self.__label = label

        self.__doubleDifferentialCrossSection = doubleDifferentialCrossSectionModule.component( )
        self.__doubleDifferentialCrossSection.setAncestor( self )

        self.__crossSection = crossSectionModule.component( )
        self.__crossSection.setAncestor( self )

        if( not isinstance( outputChannel, channelsModule.channel ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'ouputChannel not instance of channel class.' )
        self.outputChannel = outputChannel
        self.outputChannel.setAncestor( self )

        self.ENDF_MT = int( ENDF_MT )
        self.__process = process
        self.documentation = {}
        if( not( documentation is None ) ) : self.addDocumentation( documentation )

    @property
    def doubleDifferentialCrossSection( self ) :

        return( self.__doubleDifferentialCrossSection )

    @property
    def crossSection( self ) :

        return( self.__crossSection )

    @property
    def label( self ) :
        """Returns the reaction's label."""

        if( self.__label is None ) : self.updateLabel( )
        return( self.__label )

    def updateLabel( self ) :
        """Sets the reaction's label from outputChannel products."""

        self.__label = str( self.outputChannel )
        if( self.process is not None ) : self.__label += ' [%s]' % self.process

    @property
    def process( self ) :

        return( self.__process )

    @process.setter
    def process( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'process must be a string' )
        self.__process = value

    def check( self, info ):
        """
        This method is usually not called directly. Use reactionSuite.check() instead.

        Checks cross section and outputChannel (and differential cross sections if available). Checks include:
          do Z/A balance? do Q-value and thresholds agree?
          Does cross section domain agree with each product distribution/multiplicity domain?
          Does energy balance?

        @:param info: dict
        @:return list of warnings
        """

        from fudge.gnds import warning
        from . import production as productionModule
        from fudge.gnds.reactionData.doubleDifferentialCrossSection.chargedParticleElastic.CoulombPlusNuclearElastic \
            import CoulombDepositionNotSupported
        warnings = []

        reactionSuite = self.getRootAncestor( )

        def particleZA( particleID ) :

            particle = reactionSuite.PoPs[particleID]
            if hasattr(particle, 'id') and particle.id in reactionSuite.PoPs.aliases:
                particle = reactionSuite.PoPs[ particle.pid ]
            return( chemicalElementMiscPoPsModule.ZA( particle ) )

        try:
# BRB6 hardwired
            info['Q'] = self.getQ('eV', final=False)
        except ValueError:
            pass
        cpcount = sum( [ ( particleZA( prod.id ) / 1000 ) > 0 for prod in self.outputChannel ] )
        info['CoulombOutputChannel'] = cpcount > 1

        differentialCrossSectionWarnings = self.doubleDifferentialCrossSection.check( info )
        if differentialCrossSectionWarnings:
            warnings.append( warning.context("doubleDifferentialCrossSection:", differentialCrossSectionWarnings) )

        crossSectionWarnings = self.crossSection.check( info )
        if crossSectionWarnings:
            warnings.append( warning.context("Cross section:", crossSectionWarnings) )

        if 'Q' in info: del info['Q']
        del info['CoulombOutputChannel']

        if info['crossSectionOnly']:
            return warnings             # otherwise continue to check outputs

        # compare calculated and listed Q-values:
        if not isinstance(self, productionModule.production): # can't reliably compute Q for production reactions
            try:
                Q = self.getQ('eV', final=False)
                Qcalc = info['availableEnergy_eV']
                for prod in self.outputChannel:
                    try:
                        Qcalc -= prod.getMass('eV/c**2') * prod.multiplicity.getConstant()
                    except Exception:   # multiplicity is not constant
                        if( prod.id == IDsPoPsModule.photon ) : continue
                        raise ValueError("Non-constant multiplicity")
                if abs(Q-Qcalc) > PQUModule.PQU(info['dQ']).getValueAs('eV'):
                    if self.outputChannel.process != channelsModule.processes.continuum:
                        warnings.append( warning.Q_mismatch( PQUModule.PQU(Qcalc,'eV'), PQUModule.PQU(Q,'eV'), self ) )
            except ValueError:
                pass    # this test only works if multiplicity and Q are both constant for all non-gamma products

        if not (isinstance( self.outputChannel, channelsModule.sumOfRemainingOutputChannels ) or
                    self.outputChannel.isFission( ) or isinstance( self, productionModule.production ) ) :
            # check that ZA balances:
            ZAsum = 0
            for product in self.outputChannel:
                if( product.id == IDsPoPsModule.photon ) : continue
                ZAsum += particleZA( product.id ) * product.multiplicity.getConstant()
            if ZAsum != info['compoundZA']:
                warnings.append( warning.ZAbalanceWarning( self ) )

        if self.outputChannel.isFission():
            from fudge.gnds.channelData.fissionEnergyReleased import fissionEnergyReleased
            if isinstance(self.outputChannel.Q.evaluated, fissionEnergyReleased):
                info['crossSectionDomain'] = self.crossSection.domainMin, self.crossSection.domainMax
                FERwarnings = self.outputChannel.Q.evaluated.check( info )
                if FERwarnings:
                    warnings.append( warning.context("fissionEnergyReleased:", FERwarnings) )
                del info['crossSectionDomain']

        # disabling for now: only complain if distributions are missing for transportables:
        """
        if (not any( [product.distributions.components for product in self.outputChannel] ) and not any(
                [dProd.distributions.components for prod in [p for p in self.outputChannel
                    if p.outputChannel is not None] for dProd in prod.outputChannel] ) ):
            # no distributions found for any reaction product or subsequent decay product
            warnings.append( warning.noDistributions( self ) )
            return warnings """

        info['crossSectionDomain'] = self.crossSection.domainMin, self.crossSection.domainMax
        info['isTwoBody'] = isinstance( self.outputChannel, channelsModule.twoBodyOutputChannel )

        for product in self.outputChannel:
            productWarnings = product.check( info )
            if productWarnings:
                warnings.append( warning.context("Product: %s" % product.label, productWarnings) )

        del info['crossSectionDomain']
        del info['isTwoBody']

        if info['checkEnergyBalance'] and not isinstance( self, productionModule.production ):
            # Calculate energy deposition data for all products, and for decay products.
            # Then recursively check the product list for energy balance at each step of the reaction/decay
            # At each step, if we have energy deposition for *every* product in the list, we can rigorously check
            # energy balance. Otherwise, can only check that we don't exceed available energy.
            try:
                self.calculateAverageProductData( style=info['averageProductDataStyle'], **info['averageProductDataArgs'] )
            except CoulombDepositionNotSupported:
                warnings.append( warning.SkippedCoulombElasticEnergyDeposition( self ) )
            except Exception as e:
                warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                if info['failOnException']: raise
                return warnings

            def checkProductsForEnergyBalance( products, Qs, fission = False, decay = False ):
                # sample usage: for the reaction n + F19 -> n + (F19_c -> F19 + gamma), this function
                # should be called twice, to test energy balance at each step of the reaction.
                # First call: products = [n, F19_c] and Qs = [Q_reaction],
                # second call: products = [n, F19, gamma] and Qs = [Q_reaction, Q_decay].
                edepWarnings = []

                averageProductDataLabel = info['averageProductDataStyle'].label
                energyDep = []
                for prod in products:
                    if averageProductDataLabel in prod.energyDeposition:
                        energyDep.append( [ prod.label, prod.energyDeposition[ averageProductDataLabel ] ] )
                if energyDep:
                    totalEDep = energyDep[0][1]
                    for idx in range(1,len(energyDep)):
                        if(     ( totalEDep.domainMin != energyDep[idx][1].domainMin ) or
                                ( totalEDep.domainMax != energyDep[idx][1].domainMax ) ) :
                            upperEps = 0
                            if totalEDep.domainMax != energyDep[idx][1].domainMax:
                                upperEps = 1e-8
                            try:
                                totalEDep, energyDep[idx][1] = totalEDep.mutualify(
                                        1e-8,upperEps,0, energyDep[idx][1], 1e-8,upperEps,0)
                            except Exception as e:
                                edepWarnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                                if info['failOnException']: raise
                                return warnings
                        totalEDep = totalEDep + energyDep[idx][1]
                else: totalEDep = []
                Qsum = sum(Qs)

                def energyDepositedPerProduct( energyDep, ein ):
                    """ if energy doesn't balance, determine how much is deposited in each product """
                    result = []
                    availableEnergy = ein + Qsum
                    if availableEnergy == 0: availableEnergy = sum( [edep.evaluate(ein) for p,edep in energyDep] )
                    for prod, edep in energyDep:
                        edep = edep.evaluate( ein )
                        if edep is None: result.append( (prod, 0) )
                        else: result.append( ( prod, 100.0 * edep/availableEnergy ) )
                    return sorted(result, key=lambda foo: foo[1])[::-1]

                # now we have total energy deposition for all particles, use to check energy balance.
                # a few special cases to consider:
                if fission:
                    # fission products aren't listed (so far, anyway), so about 85% of available energy should be missing:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if edep > abs((ein + Qsum) * info['fissionEnergyBalanceLimit']):
                            edepWarnings.append( warning.fissionEnergyImbalance( PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                elif len(products) == len(energyDep):
                    # have full energy dep. data for all products, so we can rigorously check energy balance:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( abs(edep - (ein+Qsum)) > abs((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and abs(edep - (ein+Qsum)) > PQUModule.PQU( info['dEnergyBalanceAbsolute'] )
                                    .getValueAs(totalEDep.axes[0].unit) ):
                            edepWarnings.append( warning.energyImbalance( PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                else:
                    # missing some products, so just check that outgoing energy doesn't exceed incoming:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( (edep - (ein+Qsum)) > ((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and (edep - (ein+Qsum)) > PQUModule.PQU( info['dEnergyBalanceAbsolute'] )
                                    .getValueAs(totalEDep.axes[0].unit) ):
                            edepWarnings.append( warning.energyImbalance( PQUModule.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )

                if edepWarnings:
                    context = "Energy balance"
                    if decay: context += " (after decay)"
                    context += " for products: " + ', '.join( [prod.id for prod in products] )
                    warnings.append( warning.context( context, edepWarnings ) )

                # now recursively check decay products, if any:
                for pidx, currentProd in enumerate(products):
                    if currentProd.outputChannel is not None:
                        checkProductsForEnergyBalance(
                                products[:pidx] + [p for p in currentProd.outputChannel] + products[pidx+1:],
                                Qs + [currentProd.outputChannel.Q.getConstant()],
                                decay = True, fission = False        # FIXME what about spontaneous fission decay?
                                )
                # end of helper function checkProductsForEnergyBalance

            try:
                Q = self.outputChannel.Q.getConstant()
                checkProductsForEnergyBalance( products = [p1 for p1 in self.outputChannel], Qs = [Q],
                        fission = self.outputChannel.isFission(), decay=False )
            except ValueError:
                pass    # FIXME this test currently disabled when Q non-constant

        return warnings

    @property
    def domainMin( self ) :

        return( self.crossSection.domainMin )

    @property
    def domainMax( self ) :

        return( self.crossSection.domainMax )

    @property
    def domainUnit( self ) :

        return( self.crossSection.domainUnit )

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.crossSection.domainUnitConversionFactor( unitTo ) )

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        from . import reaction as reactionModule

        self.doubleDifferentialCrossSection.convertUnits( unitMap )
        self.crossSection.convertUnits( unitMap )
        self.outputChannel.convertUnits( unitMap )

        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.convertUnits( unitMap )
            self.availableMomentum.convertUnits( unitMap )

    def getENDL_CS_ENDF_MT( self ) :
        """
        Returns the reaction's ENDL C, S and ENDF's MT values as integers in a python dictionary with keys 
        'C', 'S' and 'MT' (e.g., { 'C' : 11, 'S' : 1, 'MT' : 53 }).
        """

        from fudge.legacy.converting import endf_endl

        MT = self.ENDF_MT
        C, S = endf_endl.getCSFromMT( MT )
        return( { 'C' : C, 'S' : S, 'MT' : MT } )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001, heatAllPoints = False,
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False ) :

        if( len( self.doubleDifferentialCrossSection ) > 0 ) :
            print( "    Skipping doubleDifferentialCrossSection reaction " )
            return

        crossSection = self.crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
                heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero, addToSuite = True )
        return( crossSection )

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def thresholdQAs( self, unit, final = True ) :

        return( self.outputChannel.thresholdQAs( unit, final = final ) )

    def getDocumentation( self, name ) :

        return( self.documentation[name] )

    def getQ( self, unit, final = True ) :
        """Returns the Q-value for this reaction. Converted to float if possible, otherwise a string value is returned."""

        return( self.thresholdQAs( unit, final = final ) )

    def getReactionSuite( self ) :

        from .. import reactionSuite as reactionSuiteModule
        return( self.findClassInAncestry( reactionSuiteModule.reactionSuite ) )

    def cullStyles( self, styleList ) :

        from . import reaction as reactionModule

#        self.__doubleDifferentialCrossSection.cullStyles( styleList )
        self.__crossSection.cullStyles( styleList )
        self.outputChannel.cullStyles( styleList )
        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.cullStyles( styleList )
            self.availableMomentum.cullStyles( styleList )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average product data.

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        verbosity = kwargs['verbosity']
        indent2 = indent + kwargs['incrementalIndent']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.outputChannel.toString( simpleString=True ) ) )

        kwargs['reaction'] = self
        kwargs['EMin'] = self.domainMin
        kwargs['EMax'] = self.domainMax
        self.outputChannel.calculateAverageProductData( style, indent2, **kwargs )

        if( hasattr( self, 'availableEnergy' ) ) :
            axes = availableEnergyModule.defaultAxes( self.domainUnit )
            QToPointwiseLinear = self.outputChannel.QToPointwiseLinear( final = True )
            availableEnergy = availableEnergyModule.XYs1d( data = [ [ QToPointwiseLinear.domainMin, QToPointwiseLinear.domainMin ], [ self.domainMax, self.domainMax ] ], axes = axes )
            availableEnergy += QToPointwiseLinear
            self.availableEnergy.add( availableEnergyModule.XYs1d( data = availableEnergy, axes = axes, label = style.label ) )
        if( hasattr( self, 'availableMomentum' ) ) :
            massInE = kwargs['projectileMass']
            availableMomentum = availableMomentumModule.calculateMomentumPoints( style, massInE, self.domainMin, self.domainMax, self.domainUnit )
            self.availableMomentum.add( availableMomentum )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.outputChannel.toString( simpleString=True ) ) )

        self.outputChannel.processMC_cdf( style, tempInfo, indent2 )

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', isPhotoAtomic = False ) :

        self.crossSection.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent, isPhotoAtomic = isPhotoAtomic )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from . import reaction as reactionModule

        tempInfo['workFile'].append( 'r%s' % tempInfo['reactionIndex'] )
        tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) + self.toString( )

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.outputChannel.toString( simpleString=True ) ) )

        tempInfo['reaction'] = self
        norm = tempInfo['groupedFlux']
        tempInfo['groupedFlux'] = None

        crossSection = style.findFormMatchingDerivedStyle( self.crossSection )
# BRB FIXME The next line is a kludge, see note on crossSection.resonancesWithBackground.processMultiGroup.
        if( isinstance( crossSection, crossSectionModule.reference ) ) :
            crossSection = crossSection.crossSection
        if( isinstance( crossSection, crossSectionModule.resonancesWithBackground ) ) :
            crossSection = crossSection.ancestor['recon']
        if( not( isinstance( crossSection, crossSectionModule.XYs1d ) ) ) :
            crossSection = crossSection.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
        tempInfo['crossSection'] = crossSection
        tempInfo['multiGroupCrossSection'] = self.crossSection.processMultiGroup( style, tempInfo, indent2 )
        self.crossSection.remove( style.label )                             # Not normalized by tempInfo['groupedFlux'] so remove.

        tempInfo['groupedFlux'] = norm
        tempInfo['multiGroupCrossSectionNormed'] = self.crossSection.processMultiGroup( style, tempInfo, indent2 )     # Normalized by tempInfo['groupedFlux'].

        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.processMultiGroup( style, tempInfo, indent2 )
            self.availableMomentum.processMultiGroup( style, tempInfo, indent2 )
        self.outputChannel.processMultiGroup( style, tempInfo, indent2 )

        del tempInfo['workFile'][-1]

    def toString( self, indent = '' ) :

        return( self.outputChannel.toString( ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        from . import reaction as reactionModule

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        attributeString = ""
        if( self.process is not None ) : attributeString += ' process="%s"' % self.process
        attributeString += ' ENDF_MT="%s"' % self.ENDF_MT

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        fissionGenre = self.outputChannel.getFissionGenre( )
        if fissionGenre is not None: attributeString += ' fissionGenre="%s"' % fissionGenre

        xmlString[-1] += attributeString + '>'

        if self.documentation:
# BRB6 What is this
            xmlString.append( '%s<documentations>' % indent2 )
            for doc in self.documentation: xmlString += self.documentation[doc].toXMLList( indent3, **kwargs )
            xmlString[-1] += '</documentations>'

        xmlString += self.doubleDifferentialCrossSection.toXMLList( indent2, **kwargs )
        xmlString += self.crossSection.toXMLList( indent2, **kwargs )
        xmlString += self.outputChannel.toXMLList( indent2, **kwargs )

        if( isinstance( self, reactionModule.reaction ) ) :
            xmlString += self.availableEnergy.toXMLList( indent2, **kwargs )
            xmlString += self.availableMomentum.toXMLList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
        """Translate a <reaction> element from xml into a reaction class instance."""

        xPath.append( '%s[@label="%s"]' % ( element.tag, element.get( 'label' ) ) )

        crossSectionComponent = fudge.gnds.reactionData.crossSection.parseXMLNode(
                element.find( 'crossSection' ), xPath, linkData )

        outputChannel = fudge.gnds.channels.parseXMLNode( element.find('outputChannel'), xPath, linkData)
        outputChannel.fissionGenre = element.get( 'fissionGenre' )

        reaction = cls( outputChannel = outputChannel, label = element.get( 'label' ), ENDF_MT = int( element.get('ENDF_MT') ) )
        if element.get('process'):
            reaction.process = element.get('process')
        for crossSection in crossSectionComponent : reaction.crossSection.add( crossSection )

        if( element.find( 'documentations' ) ) :
            for doc in element.find( 'documentations' ) :
                reaction.addDocumentation( fudge.gnds.documentation.documentation.parseXMLNode( doc, xPath, linkData))

        if( element.find( 'doubleDifferentialCrossSection' ) ) :
            reaction.doubleDifferentialCrossSection.parseXMLNode( element.find( 'doubleDifferentialCrossSection' ), xPath, linkData )

        if( element.find( availableEnergyModule.component.moniker ) ) :
            reaction.availableEnergy = availableEnergyModule.parseXMLNode(
                element.find( availableEnergyModule.component.moniker ), xPath, linkData )

        if( element.find( availableMomentumModule.component.moniker ) ) :
            reaction.availableMomentum = availableMomentumModule.parseXMLNode(
                element.find( availableMomentumModule.component.moniker ), xPath, linkData )

        xPath.pop( )
        return( reaction )

def isGNDSReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, base_reaction)
