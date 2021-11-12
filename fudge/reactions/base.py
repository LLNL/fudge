# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideModule

import xData.ancestry as ancestryModule
from xData import standards as standardsModule

import fudge

from .. import outputChannel as outputChannelModule 
from ..reactionData.doubleDifferentialCrossSection import doubleDifferentialCrossSection as doubleDifferentialCrossSectionModule
from ..reactionData import crossSection as crossSectionModule
from ..reactionData import availableEnergy as availableEnergyModule
from ..reactionData import availableMomentum as availableMomentumModule

__metaclass__ = type

class FissionGenre :

    total = 'total'
    firstChance = 'firstChance'
    secondChance = 'secondChance'
    thirdChance = 'thirdChance'
    fourthChance = 'fourthChance'

    allowed = ( None, total, firstChance, secondChance, thirdChance, fourthChance )

class base_reaction( ancestryModule.ancestry ) :
    """Base class for all types of reaction."""

    ancestryMembers = ( 'crossSection', 'doubleDifferentialCrossSection', 'outputChannel' )

    def __init__( self, genre, ENDF_MT, fissionGenre = None, documentation = None, label = None ) :

        ancestryModule.ancestry.__init__( self )
        self.__label = label

        self.__doubleDifferentialCrossSection = doubleDifferentialCrossSectionModule.component( )
        self.__doubleDifferentialCrossSection.setAncestor( self )

        self.__crossSection = crossSectionModule.component( )
        self.__crossSection.setAncestor( self )

        self.__outputChannel = outputChannelModule.outputChannel( genre )
        self.__outputChannel.setAncestor( self )

        self.fissionGenre = fissionGenre
        self.ENDF_MT = int( ENDF_MT )

        self.documentation = {}
        if( not( documentation is None ) ) : self.addDocumentation( documentation )

    def __str__( self ) :

        return( self.label )

    @property
    def doubleDifferentialCrossSection( self ) :

        return( self.__doubleDifferentialCrossSection )

    @property
    def crossSection( self ) :

        return( self.__crossSection )

    @property
    def outputChannel( self ) :

        return( self.__outputChannel )

    @property
    def fissionGenre( self ) :

        return( self.__fissionGenre )

    @fissionGenre.setter
    def fissionGenre( self, value ) :

        if( value not in FissionGenre.allowed ) : raise Exception( 'Invalid fission genre "%s".' % value )
        self.__fissionGenre = value

    @property
    def label( self ) :
        """Returns the reaction's label."""

        if( self.__label is None ) : self.updateLabel( )
        return( self.__label )

    def updateLabel( self ) :
        """Sets the reaction's label from outputChannel products."""

        label = self.__outputChannel.toString( MT = self.ENDF_MT )
        if self.fissionGenre is not None:
            if '[' in label: raise Exception('Label already contains a process: "%s"' % label)
            if len(label) == 0:
                label = self.fissionGenre + ' fission'
            else:
                label += ' [%s fission]' % self.fissionGenre
        self.__label = label

    def findLinks( self, links ) :

        for ancestryMember in self.ancestryMembers :
            if( ancestryMember in ( 'doubleDifferentialCrossSection' ) ) : continue
            getattr( self, ancestryMember ).findLinks( links )

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def isFission( self ):

        return self.__fissionGenre is not None

    def isThermalNeutronScatteringLaw( self ) :

        for form in self.__doubleDifferentialCrossSection :
            if( form.isThermalNeutronScatteringLaw( ) ) : return( True )
        return( False )

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

        from fudge import warning
        from . import production as productionModule
        from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic.CoulombPlusNuclearElastic \
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
        cpcount = sum( [ ( particleZA( prod.id ) // 1000 ) > 0 for prod in self.__outputChannel ] )
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
                for prod in self.__outputChannel:
                    try:
                        Qcalc -= prod.getMass('eV/c**2') * prod.multiplicity.getConstant()
                    except Exception:   # multiplicity is not constant
                        if( prod.id == IDsPoPsModule.photon ) : continue
                        raise ValueError("Non-constant multiplicity")
                if abs(Q-Qcalc) > PQUModule.PQU(info['dQ']).getValueAs('eV'):
                    if self.__outputChannel.process != outputChannelModule.processes.continuum:
                        warnings.append( warning.Q_mismatch( PQUModule.PQU(Qcalc,'eV'), PQUModule.PQU(Q,'eV'), self ) )
            except ValueError:
                pass    # this test only works if multiplicity and Q are both constant for all non-gamma products

        if not (self.__outputChannel.genre == outputChannelModule.Genre.sumOfRemainingOutputChannels or
                    self.isFission( ) or isinstance( self, productionModule.production) ):
            # check that ZA balances:
            ZAsum = 0
            for product in self.__outputChannel:
                if( product.id == IDsPoPsModule.photon ) : continue
                ZAsum += particleZA( product.id ) * product.multiplicity.getConstant()
            if ZAsum != info['compoundZA']:
                warnings.append( warning.ZAbalanceWarning( self ) )

        # disabling for now: only complain if distributions are missing for transportables:
        """
        if (not any( [product.distributions.components for product in self.__outputChannel] ) and not any(
                [dProd.distributions.components for prod in [p for p in self.__outputChannel
                    if p.outputChannel is not None] for dProd in prod.outputChannel] ) ):
            # no distributions found for any reaction product or subsequent decay product
            warnings.append( warning.noDistributions( self ) )
            return warnings """

        info['crossSectionDomain'] = self.crossSection.domainMin, self.crossSection.domainMax
        info['isTwoBody'] = self.__outputChannel.genre == outputChannelModule.Genre.twoBody

        for product in self.__outputChannel:
            productWarnings = product.check( info )
            if productWarnings:
                warnings.append( warning.context("Product: %s" % product.label, productWarnings) )

        fissionFragmentWarnings = self.__outputChannel.fissionFragmentData.check( info )
        if fissionFragmentWarnings:
            warnings.append( warning.context("Fission fragment info:", fissionFragmentWarnings) )

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
                Q = self.__outputChannel.Q.getConstant()
                checkProductsForEnergyBalance( products = [p1 for p1 in self.__outputChannel], Qs = [Q],
                        fission = self.isFission(), decay=False )
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

        self.__doubleDifferentialCrossSection.convertUnits( unitMap )
        self.__crossSection.convertUnits( unitMap )
        self.__outputChannel.convertUnits( unitMap )

        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.convertUnits( unitMap )
            self.availableMomentum.convertUnits( unitMap )

    def amendForPatch( self, fromLabel, toLabel ) :

        from . import reaction as reactionModule

        self.__doubleDifferentialCrossSection.amendForPatch( fromLabel, toLabel )
        self.__crossSection.amendForPatch( fromLabel, toLabel )
        self.__outputChannel.amendForPatch( fromLabel, toLabel )

        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.amendForPatch( fromLabel, toLabel )
            self.availableMomentum.amendForPatch( fromLabel, toLabel )

    def diff( self, other, diffResults ) :

        self.__crossSection.diff( other.crossSection, diffResults )
        self.__outputChannel.diff( other.outputChannel, diffResults )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001, heatAllPoints = False,
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False, verbose = 0 ) :

        if( len( self.doubleDifferentialCrossSection ) > 0 ) :
            if( verbose > 0 ) : print( "    Skipping doubleDifferentialCrossSection reaction " )
            return

        crossSection = self.crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
                heatBelowThreshold, heatAllEDomain, setThresholdToZero = setThresholdToZero, addToSuite = True )
        return( crossSection )

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def thresholdQAs( self, unit, final = True ) :

        return( self.__outputChannel.thresholdQAs( unit, final = final ) )

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
        self.__outputChannel.cullStyles( styleList )
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

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.__outputChannel.toString( simpleString = True, MT = self.ENDF_MT ) ) )

        kwargs['reaction'] = self
        kwargs['EMin'] = self.domainMin
        kwargs['EMax'] = self.domainMax
        self.__outputChannel.calculateAverageProductData( style, indent2, **kwargs )

        if( hasattr( self, 'availableEnergy' ) ) :
            axes = availableEnergyModule.defaultAxes( self.domainUnit )
            QToPointwiseLinear = self.__outputChannel.QToPointwiseLinear( final = True )
            availableEnergy = availableEnergyModule.XYs1d( data = [ [ QToPointwiseLinear.domainMin, QToPointwiseLinear.domainMin ], [ self.domainMax, self.domainMax ] ], axes = axes )
            availableEnergy += QToPointwiseLinear
            self.availableEnergy.add( availableEnergyModule.XYs1d( data = availableEnergy, axes = axes, label = style.label ) )
        if( hasattr( self, 'availableMomentum' ) ) :
            massInE = kwargs['projectileMass']
            availableMomentum = availableMomentumModule.calculateMomentumPoints( style, massInE, self.domainMin, self.domainMax, self.domainUnit )
            self.availableMomentum.add( availableMomentum )

    def partialProductionIntegral( self, reaction_suite, productID, energyIn, energyOut = None, muOut = None, phiOut = None, 
            frame = standardsModule.frames.labToken, LegendreOrder = 0, **kwargs ) :

        if( isinstance( self.crossSection[0], crossSectionModule.CoulombPlusNuclearElastic ) ) :
            issueCounters = kwargs.get( 'issueCounters', None )
            if( issueCounters is not None ) : issueCounter = issueCounters.get( 'partialProductionIntegral:CoulombPlusNuclearElastic', 1, True )
            status = issueCounter.increment( )
            if( status ) : print( '    WARNING: partialProductionIntegral: skipping elastic Coulomb reaction %s.' % self )
            return( 0.0 )

        crossSection = self.crossSection.evaluate( energyIn )
        if( crossSection == 0.0 ) : return( 0.0 )
        _partialProductionIntegral = self.__outputChannel.partialProductionIntegral( reaction_suite, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( crossSection * _partialProductionIntegral )

    def processMC_cdf( self, style, tempInfo, indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.__outputChannel.toString( simpleString = True, MT = self.ENDF_MT ) ) )

        self.__outputChannel.processMC_cdf( style, tempInfo, indent2 )

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', isPhotoAtomic = False ) :

        self.crossSection.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent, isPhotoAtomic = isPhotoAtomic )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from . import reaction as reactionModule

        tempInfo['workFile'].append( 'r%s' % tempInfo['reactionIndex'] )
        tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) + self.toString( )

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        if( verbosity > 0 ) : print( '%s%s' % (indent, self.__outputChannel.toString( simpleString = True, MT = self.ENDF_MT ) ) )

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
        self.__outputChannel.processMultiGroup( style, tempInfo, indent2 )

        del tempInfo['workFile'][-1]

    def removeStyles( self, styleLabels ) :

        from . import reaction as reactionModule

        self.__doubleDifferentialCrossSection.removeStyles( styleLabels )
        self.__crossSection.removeStyles( styleLabels )
        self.__outputChannel.removeStyles( styleLabels )
        if( isinstance( self, reactionModule.reaction ) ) :
            self.availableEnergy.removeStyles( styleLabels )
            self.availableMomentum.removeStyles( styleLabels )

    def toString( self, indent = '' ) :

        return( self.__outputChannel.toString( indent = indent, MT = self.ENDF_MT ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        from . import reaction as reactionModule

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        attributeString = ""
        attributeString += ' ENDF_MT="%s"' % self.ENDF_MT
        if( self.fissionGenre is not None ) : attributeString += ' fissionGenre="%s"' % str( self.fissionGenre )

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        xmlString[-1] += attributeString + '>'

        if self.documentation:
# BRB6 What is this
            xmlString.append( '%s<documentations>' % indent2 )
            for doc in self.documentation: xmlString += self.documentation[doc].toXMLList( indent3, **kwargs )
            xmlString[-1] += '</documentations>'

        xmlString += self.__doubleDifferentialCrossSection.toXMLList( indent2, **kwargs )
        xmlString += self.__crossSection.toXMLList( indent2, **kwargs )
        xmlString += self.__outputChannel.toXMLList( indent2, **kwargs )

        if( isinstance( self, reactionModule.reaction ) ) :
            xmlString += self.availableEnergy.toXMLList( indent2, **kwargs )
            xmlString += self.availableMomentum.toXMLList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def parseNode( self, node, xPath, linkData, **kwargs ) :

        xPath.append( node.tag )

        if( node.find( 'documentations' ) ) :
            for doc in node.find( 'documentations' ) :
                self.addDocumentation(fudge.documentation.documentation.parseXMLNode(doc, xPath, linkData))

        self.doubleDifferentialCrossSection.parseXMLNode( node.find( doubleDifferentialCrossSectionModule.component.moniker ), xPath, linkData )
        self.crossSection.parseXMLNode( node.find( crossSectionModule.component.moniker ), xPath, linkData )

        if( node.find( outputChannelModule.outputChannel.moniker ) ) :
            self.outputChannel.parseXMLNode( node.find( outputChannelModule.outputChannel.moniker ), xPath, linkData )

        if( node.find( availableEnergyModule.component.moniker ) ) :
            self.availableEnergy = availableEnergyModule.component( )
            self.availableEnergy.setAncestor( self )
            self.availableEnergy.parseXMLNode( node.find( availableEnergyModule.component.moniker ), xPath, linkData )

        if( node.find( availableMomentumModule.component.moniker ) ) :
            self.availableMomentum = availableMomentumModule.component( )
            self.availableMomentum.setAncestor( self )
            self.availableMomentum.parseXMLNode( node.find( availableMomentumModule.component.moniker ), xPath, linkData )

        xPath.pop( )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
        """Translate a <reaction> element from xml into a reaction class instance."""

        xPath.append( '%s[@label="%s"]' % ( element.tag, element.get( 'label' ) ) )
        reaction = cls( outputChannelModule.Genre.NBody, label = element.get( 'label' ), ENDF_MT = int( element.get( 'ENDF_MT' ) ), 
                fissionGenre = element.get( 'fissionGenre' ) )
        xPath.pop( )

        reaction.parseNode( element, xPath, linkData )

        return( reaction )

def isGNDSReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, base_reaction)
