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

from fudge.core.ancestry import ancestry
import fudge

__metaclass__ = type

reactionToken = 'reaction'
productionToken = 'production'
summedReactionToken = 'summedReaction'
fissionComponentToken = 'fissionComponent'
partialGammaProductionToken = 'partialGammaProduction'

class base_reaction( ancestry ):
    """Base class for all types of reaction."""

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
            # Calculate energy deposition data for all products, and for decay products.
            # Then recursively check the product list for energy balance at each step of the reaction/decay
            # At each step, if we have energy deposition for *every* product in the list, we can rigorously check
            # energy balance. Otherwise, can only check that we don't exceed available energy.
            try:
                self.calculateDepositionData( info['processInfo'], '' )
            except Exception, e:
                warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                return warnings

            def checkProductsForEnergyBalance( products, Qs, flags = {'fission':False, 'decay':False} ):
                # sample usage: for the reaction n + F19 -> n + (F19_c -> F19 + gamma), this function
                # should be called twice, to test energy balance at each step of the reaction.
                # First call: products = [n, F19_c] and Qs = [Q_reaction],
                # second call: products = [n, F19, gamma] and Qs = [Q_reaction, Q_decay].
                edepWarnings = []

                energyDep = [[prod.getLabel(), prod.data['depositionEnergy'].getNativeData()]
                        for prod in products if 'depositionEnergy' in prod.data]
                if energyDep:
                    totalEDep = energyDep[0][1]
                    for idx in range(1,len(energyDep)):
                        if totalEDep.getDomain() != energyDep[idx][1].getDomain():
                            try:
                                totalEDep, energyDep[idx][1] = totalEDep.mutualify(
                                        1e-8,0,0, energyDep[idx][1], 1e-8,0,0)
                            except Exception, e:
                                edepWarnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                                return warnings
                        totalEDep += energyDep[idx][1]
                else: totalEDep = []
                Qsum = sum(Qs)

                def energyDepositedPerProduct( energyDep, ein ):
                    """ if energy doesn't balance, determine how much is deposited in each product """
                    result = []
                    availableEnergy = ein + Qsum
                    if availableEnergy == 0: availableEnergy = sum( [edep.getValue(ein) for p,edep in energyDep] )
                    for prod, edep in energyDep:
                        edep = edep.getValue( ein )
                        if edep is None: result.append( (prod, 0) )
                        else: result.append( ( prod, 100.0 * edep/availableEnergy ) )
                    return sorted(result, key=lambda foo: foo[1])[::-1]

                # now we have total energy deposition for all particles, use to check energy balance.
                # a few special cases to consider:
                if flags['fission']:
                    # fission products aren't listed (so far, anyway), so about 85% of available energy should be missing:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if edep > abs((ein + Qsum) * info['fissionEnergyBalanceLimit']):
                            edepWarnings.append( warning.fissionEnergyImbalance( PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                elif len(products) == len(energyDep):
                    # have full energy dep. data for all products, so we can rigorously check energy balance:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( abs(edep - (ein+Qsum)) > abs((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and abs(edep - (ein+Qsum)) > info['dEnergyBalanceAbsolute'] ):
                            edepWarnings.append( warning.energyImbalance( PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                else:
                    # missing some products, so just check that outgoing energy doesn't exceed incoming:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( (edep - (ein+Qsum)) > ((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and (edep - (ein+Qsum)) > info['dEnergyBalanceAbsolute'] ):
                            edepWarnings.append( warning.energyImbalance( PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )

                if edepWarnings:
                    context = "Energy balance"
                    if flags['decay']: context += " (after decay)"
                    context += " for products: " + ', '.join( [prod.particle.name for prod in products] )
                    warnings.append( warning.context( context, edepWarnings ) )

                # now recursively check decay products, if any:
                for pidx, currentProd in enumerate(products):
                    if currentProd.decayChannel is not None:
                        flags['decay'] = True
                        checkProductsForEnergyBalance(
                                products[:pidx] + [p for p in currentProd.decayChannel] + products[pidx+1:],
                                Qs + [currentProd.decayChannel.getConstantQAs('eV')],
                                flags
                                )
                # end of helper function checkProductsForEnergyBalance

            checkProductsForEnergyBalance( products = [p1 for p1 in self.outputChannel], Qs = [Q],
                    flags = {'fission':self.outputChannel.isFission(), 'decay':False} )

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
        """Returns the Q-value for this reaction. Converted to float if possible, otherwise a string value is returned."""

        if( 'constant' in self.Q.forms ) :
            return( self.Q.getConstantAs( unit ) )
        else :      # BRB ?????? getQ needs work with redesign of reaction stuff.
            raise Exception

    def setQ( self, Q ) :

        self.Q.setParent( self )
        self.Q = Q

    def toXMLList( self, flags, indent = '' ):

        attributeString = ""
        for attribute in self.attributes : attributeString += ' %s="%s"' % ( attribute, self.attributes[attribute] )

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
# BRB ????????? self should always have an outputChannel
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
            xmlString += self.outputChannel.toXMLList( flags, indent = indent+'  ' )
            for key in self.data: xmlString += self.data[key].toXMLList( indent = indent+'  ' )
        elif hasattr(self, 'Q'):
            xmlString += self.Q.toXMLList( indent = indent+'  ' )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

def isGNDReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, base_reaction)
