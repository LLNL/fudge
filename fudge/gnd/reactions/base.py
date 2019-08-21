# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry
import fudge

__metaclass__ = type

reactionToken = 'reaction'
productionToken = 'production'
summedReactionToken = 'summedReaction'
fissionComponentToken = 'fissionComponent'
partialGammaProductionToken = 'partialGammaProduction'

class base_reaction( ancestry ):
    """ base class for all types of reaction """

    def __init__( self, token, label, ENDF_MT, crossSection = None, documentation = None, attributes = {} ):

        ancestry.__init__( self, token, None, attribute = 'label' )
        self.setLabel( label )
        self.crossSection = None
        self.setCrossSection( crossSection )
        self.attributes = {}
        self.setAttribute( 'ENDF_MT', '%d' % ENDF_MT )
        for attribute in attributes : self.setAttribute( attribute, attributes[attribute] )
        self.documentation = {}
        if( not( documentation is None ) ) : addDocumentation( documentation )

    def __checkCrossSection__( self, info ):
        """Check for Z/A balance, correct Q/Thresholds, plus check within the reaction
        for problems with cross section. """

        from fudge.gnd import warning
        warnings = []

        info['Q'] = self.getQ('eV')
        crossSectionWarnings = self.crossSection.check( info )
        if crossSectionWarnings:
            warnings.append( warning.context("Cross section:", crossSectionWarnings) )
        del info['Q']

        return warnings

    def __checkOutputChannel__( self, info ):
        """Check for problems in reaction outputChannel, distributions and multiplicities,
        energy balance. """

        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []

        # compare calculated and listed Q-values:
        Q = self.getQ('eV')
        try:
            Qcalc = info['availableEnergy_eV']
            for prod in self.outputChannel:
                try:
                    Qcalc -= prod.getMass('eV/c**2') * prod.multiplicity.getConstant()
                except Exception:   # multiplicity is not constant
                    if prod.getName() == 'gamma': continue
                    raise ValueError, "Non-constant multiplicity"
            if abs(Q-Qcalc) > PQU(info['dQ']).getValueAs('eV'):
                warnings.append( warning.Q_mismatch( PQU(Qcalc,'eV'), PQU(Q,'eV'), self ) )
        except ValueError:
            pass    # this test only works if multiplicity is constant
 
        if not (self.outputChannel.genre is fudge.gnd.channels.sumOfRemainingOutputChannelsGenre
                or self.outputChannel.isFission()):
            # check that ZA balances:
            ZAsum = 0
            for product in self.outputChannel:
                if product.getName() == 'gamma': continue
                ZAsum += product.particle.getZ_A_SuffixAndZA()[3] * product.multiplicity.getConstant()
            if ZAsum != info['compoundZA']:
                warnings.append( warning.ZAbalanceWarning( self ) )

        if self.outputChannel.fissionEnergyReleased is not None:
            info['crossSectionDomain'] = self.crossSection.getDomain(asPQU=True)
            FERwarnings = self.outputChannel.fissionEnergyReleased.check( info )
            if FERwarnings:
                warnings.append( warning.context("fissionEnergyReleased:", FERwarnings) )
            del info['crossSectionDomain']

        # disabling for now: only complain if distributions are missing for transportables:
        """
        if (not any( [product.distributions.components for product in self.outputChannel] ) and not any(
                [dProd.distributions.components for prod in [p for p in self.outputChannel
                    if p.decayChannel is not None] for dProd in prod.decayChannel] ) ):
            # no distributions found for any reaction product or subsequent decay product
            warnings.append( warning.noDistributions( self ) )
            return warnings """

        info['crossSectionDomain'] = self.crossSection.getDomain()
        info['isTwoBody'] = (self.outputChannel.genre is fudge.gnd.channels.twoBodyGenre)

        for product in self.outputChannel:
            productWarnings = product.check( info )
            if productWarnings:
                warnings.append( warning.context("Product: %s" % product.getLabel(), productWarnings) )

        del info['crossSectionDomain']
        del info['isTwoBody']

        if info['checkEnergyBalance']:
            # get list of energy deposition data for all products. If a product has no energy dep. data,
            # but its decay product does, use that instead. Also determine whether reaction has *full*
            # energy deposition info: do all products (or all their decay prods) have energy dep info?
            try:
                self.calculateDepositionData( info['processInfo'], '' )
            except Exception, e:
                warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                return warnings

            energyDep = []
            energyDepCount = 0
            for product in self.outputChannel:
                if 'depositionEnergy' in product.data:
                    energyDep.append( [product.getLabel(), product.data['depositionEnergy']['pointwise'] ] )
                    energyDepCount += 1
                elif product.decayChannel is not None:
                    decayEnergyDep = [[dProd.getLabel(), dProd.data['depositionEnergy']['pointwise']]
                            for dProd in product.decayChannel if 'depositionEnergy' in dProd.data]
                    energyDep.extend( decayEnergyDep )
                    if len(decayEnergyDep) == len(product.decayChannel): energyDepCount += 1
            if energyDep:
                totalEDep = energyDep[0][1]
                for idx in range(1,len(energyDep)):
                    if totalEDep.getDomain() != energyDep[idx][1].getDomain():
                        try:
                            totalEDep, energyDep[idx][1] = totalEDep.mutualify(1e-8,0,0, energyDep[idx][1], 1e-8,0,0)
                        except Exception, e:
                            warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                            return warnings
                    totalEDep += energyDep[idx][1]
            else: totalEDep = []

            def energyDepositedPerProduct( energyDep, ein ):
                """ if energy doesn't balance, determine how much is deposited in each product """
                result = []
                for prod, edep in energyDep:
                    edep = edep.getValue( ein )
                    if edep is None: result.append( (prod, 0) )
                    else: result.append( ( prod, 100.0 * edep/(ein + Q) ) )
                return sorted(result, key=lambda foo: foo[1])[::-1]

            # now we have total energy deposition for all particles, use to check energy balance.
            # a few special cases to consider:
            if self.outputChannel.isFission():
                # fission products aren't listed (so far, anyway), so about 85% of available energy should be missing:
                for i,(ein,edep) in enumerate(totalEDep):
                    if edep > abs((ein + Q) * info['fissionEnergyBalanceLimit']):
                        warnings.append( warning.fissionEnergyImbalance( PQU(ein, totalEDep.axes[0].unit),
                            i, energyDepositedPerProduct(energyDep, ein), self ) )
            elif len(self.outputChannel) == energyDepCount:
                # have full energy dep. data for all products, so we can rigorously check energy balance:
                for i,(ein,edep) in enumerate(totalEDep):
                    if ( abs(edep - (ein+Q)) > abs((ein+Q) * info['dEnergyBalanceRelative'])
                            and abs(edep - (ein+Q)) > info['dEnergyBalanceAbsolute'] ):
                        warnings.append( warning.energyImbalance( PQU(ein, totalEDep.axes[0].unit),
                            i, energyDepositedPerProduct(energyDep, ein), self ) )
            else:
                # missing some products, so just check that outgoing energy doesn't exceed incoming:
                for i,(ein,edep) in enumerate(totalEDep):
                    if ( (edep - (ein+Q)) > ((ein+Q) * info['dEnergyBalanceRelative'])
                            and (edep - (ein+Q)) > info['dEnergyBalanceAbsolute'] ):
                        warnings.append( warning.energyImbalance( PQU(ein, totalEDep.axes[0].unit),
                            i, energyDepositedPerProduct(energyDep, ein), self ) )

        return warnings

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.crossSection.domainUnitConversionFactor( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.crossSection.getDomainUnit( ) )

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )

    def setAttribute( self, name, value ) :
        """Adds attribute name and its value to the list of attributes."""

        self.attributes[name] = value

    def getCrossSection( self ) :
        """Returns the reaction's cross section."""

        return( self.crossSection )

    def setCrossSection( self, crossSection ) :
        """Sets channel's cross section to crossSection."""

        if( self.crossSection is not None ) : self.crossSection.setParent( None )
        self.crossSection = crossSection
        if( crossSection is not None ) : crossSection.setParent( self )

    def getENDL_CS_ENDF_MT( self ) :
        """Returns the reaction's ENDL C, S and ENDF's MT values as integers in a python dictionary with keys 'C', 'S' and 'MT' 
        (e.g., { 'C' : 11, 'S' : 1, 'MT' : 53 })."""

        from fudge.legacy.converting import endf_endl

        MT = int( self.attributes['ENDF_MT'] )
        C, S = endf_endl.getCSFromMT( MT )
        return( { 'C' : C, 'S' : S, 'MT' : MT } )

    def getLabel( self ) :
        """Returns the reaction's label."""

        return( self.label )
    
    def setLabel( self, label ) :
        """Sets the reaction's label."""

        self.label = label

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def getDocumentation( self, name ) :

        return( self.documentation[doc.name] )

    def getQ( self, unit, final = True, groundStateQ = False ) :
        """Returns the Q-value for this reaction. Converted to float if possible, otherwis a string value is returned"""

        if 'constant' in self.Q.forms: return self.Q.getConstantAs( unit )
        else:
            raise Exception

    def setQ( self, Q ) :
        self.Q.setParent( self )
        self.Q = Q

    def toXMLList( self, flags, verbosityIndent = '', indent = '' ):

        if( flags['verbosity'] >= 10 ) :
            print '%s%s:' % ( verbosityIndent, self.moniker ),
            if hasattr(self, 'outputChannel'): print self.outputChannel.toString( simpleString = True )
            elif self.moniker == productionToken: print self
            elif self.moniker == summedReactionToken: print self.name

        attributeString = ""
        for attribute in self.attributes : attributeString += ' %s="%s"' % ( attribute, self.attributes[attribute] )

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        if hasattr(self, 'outputChannel'):  # for reaction, fissionComponent and partialGammaProduction
            xmlString[-1] += ' outputChannel="%s"' % self.outputChannel
            fissionGenre = self.outputChannel.getFissionGenre( )
            if fissionGenre is not None: attributeString += ' fissionGenre="%s"' % fissionGenre
        elif self.moniker == productionToken:
            xmlString[-1] += ' product="%s" Q="%s"' % (self.product, self.Q.getXMLAttribute() )
        elif self.moniker == summedReactionToken:
            xmlString[-1] += ' name="%s" Q="%s"' % (self.name, self.Q.getXMLAttribute() )

        xmlString[-1] += attributeString + '>'

        if self.documentation:
            indent2 = indent+'  '
            xmlString.append( '%s<documentations>' % indent2 )
            for doc in self.documentation: xmlString += self.documentation[doc].toXMLList( indent = indent2+'  ' )
            xmlString[-1] += '</documentations>'

        elif self.moniker == summedReactionToken:
            for summand in self.summands: xmlString.append( summand.toXML(indent+'  ') )

        xmlString += self.getCrossSection( ).toXMLList( indent = indent+'  ' )
        if hasattr(self, 'outputChannel'):
            xmlString += self.outputChannel.toXMLList( flags, verbosityIndent + '    ', indent = indent+'  ' )
            for key in self.data: xmlString += self.data[key].toXMLList( indent = indent+'  ' )
        elif hasattr(self, 'Q'):
            xmlString += self.Q.toXMLList( indent = indent+'  ' )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

def isGNDReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, base_reaction)
