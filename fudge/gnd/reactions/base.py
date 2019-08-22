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

import xData.ancestry as ancestryModule

from ..reactionData import crossSection as crossSectionModule

import fudge

__metaclass__ = type

class base_reaction( ancestryModule.ancestry ):
    """Base class for all types of reaction."""

    def __init__( self, label, ENDF_MT, documentation = None, date = None, EFL = None ) :

        ancestryModule.ancestry.__init__( self )
        self.setLabel( label )

        self.crossSection = crossSectionModule.component( )
        self.crossSection.setAncestor( self )

        self.date = date
        self.EFL = EFL
        self.ENDF_MT = int( ENDF_MT )
        self.__process = None
        self.documentation = {}
        if( not( documentation is None ) ) : self.addDocumentation( documentation )

    @property
    def process( self ) :

        return( self.__process )

    @process.setter
    def process( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'process must be a string' )
        self.__process = value

    def __checkCrossSection__( self, info ):
        """Check for Z/A balance, correct Q/Thresholds, plus check within the reaction
        for problems with cross section. """

        from fudge.gnd import warning
        warnings = []

        try:
            info['Q'] = self.getQ('eV')
        except ValueError:
            pass
        CoulombOutputChannel = False
        if hasattr(self,'outputChannel'):
            cpcount = sum( [prod.particle.getZ_A_SuffixAndZA()[0]>0 for prod in self.outputChannel] )
            if cpcount > 1:
                CoulombOutputChannel = True
        info['CoulombOutputChannel'] = CoulombOutputChannel

        crossSectionWarnings = self.crossSection.check( info )
        if crossSectionWarnings:
            warnings.append( warning.context("Cross section:", crossSectionWarnings) )

        if 'Q' in info: del info['Q']
        del info['CoulombOutputChannel']

        return warnings

    def __checkOutputChannel__( self, info ):
        """Check for problems in reaction outputChannel, distributions and multiplicities,
        energy balance. """

        from fudge.gnd import warning
        from pqu import PQU
        from ..productData.distributions.angular import CoulombDepositionNotSupported
        warnings = []

        # compare calculated and listed Q-values:
        try:
            Q = self.getQ('eV')
            Qcalc = info['availableEnergy_eV']
            for prod in self.outputChannel:
                try:
                    Qcalc -= prod.getMass('eV/c**2') * prod.multiplicity.getConstant()
                except Exception:   # multiplicity is not constant
                    if( prod.name == 'gamma' ) : continue
                    raise ValueError, "Non-constant multiplicity"
            if abs(Q-Qcalc) > PQU.PQU(info['dQ']).getValueAs('eV'):
                warnings.append( warning.Q_mismatch( PQU.PQU(Qcalc,'eV'), PQU.PQU(Q,'eV'), self ) )
        except ValueError:
            pass    # this test only works if multiplicity and Q are both constant
 
        if not (self.outputChannel.genre is fudge.gnd.channels.sumOfRemainingOutputChannelsGenre
                or self.outputChannel.isFission()):
            # check that ZA balances:
            ZAsum = 0
            for product in self.outputChannel:
                if( product.name == 'gamma' ) : continue
                ZAsum += product.particle.getZ_A_SuffixAndZA()[3] * product.multiplicity.getConstant()
            if ZAsum != info['compoundZA']:
                warnings.append( warning.ZAbalanceWarning( self ) )

        if self.outputChannel.fissionEnergyReleased is not None:
            info['crossSectionDomain'] = self.crossSection.domain(asPQU=True)
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

        info['crossSectionDomain'] = self.crossSection.domain()
        info['isTwoBody'] = (self.outputChannel.genre is fudge.gnd.channels.twoBodyGenre)

        for product in self.outputChannel:
            productWarnings = product.check( info )
            if productWarnings:
                warnings.append( warning.context("Product: %s" % product.label, productWarnings) )

        del info['crossSectionDomain']
        del info['isTwoBody']

        if info['checkEnergyBalance']:
            # Calculate energy deposition data for all products, and for decay products.
            # Then recursively check the product list for energy balance at each step of the reaction/decay
            # At each step, if we have energy deposition for *every* product in the list, we can rigorously check
            # energy balance. Otherwise, can only check that we don't exceed available energy.
            try:
                self.calculateDepositionData( info['processInfo'], info['tempInfo'], '' )
            except CoulombDepositionNotSupported, e:
                warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
            except Exception, e:
                warnings.append( warning.EnergyDepositionExceptionRaised( str(e), self ) )
                if info['failOnException']: raise
                return warnings

            def checkProductsForEnergyBalance( products, Qs, flags = {'fission':False, 'decay':False} ):
                # sample usage: for the reaction n + F19 -> n + (F19_c -> F19 + gamma), this function
                # should be called twice, to test energy balance at each step of the reaction.
                # First call: products = [n, F19_c] and Qs = [Q_reaction],
                # second call: products = [n, F19, gamma] and Qs = [Q_reaction, Q_decay].
                edepWarnings = []

                depositionStyle = info['processInfo'].style.label
                energyDep = []
                for prod in products:
                    if depositionStyle in prod.energyDeposition:
                        energyDep.append( [ prod.label, prod.energyDeposition[ depositionStyle ] ] )
                if energyDep:
                    totalEDep = energyDep[0][1]
                    for idx in range(1,len(energyDep)):
                        if totalEDep.domain() != energyDep[idx][1].domain():
                            try:
                                totalEDep, energyDep[idx][1] = totalEDep.mutualify(
                                        1e-8,0,0, energyDep[idx][1], 1e-8,0,0)
                            except Exception, e:
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
                if flags['fission']:
                    # fission products aren't listed (so far, anyway), so about 85% of available energy should be missing:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if edep > abs((ein + Qsum) * info['fissionEnergyBalanceLimit']):
                            edepWarnings.append( warning.fissionEnergyImbalance( PQU.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                elif len(products) == len(energyDep):
                    # have full energy dep. data for all products, so we can rigorously check energy balance:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( abs(edep - (ein+Qsum)) > abs((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and abs(edep - (ein+Qsum)) > info['dEnergyBalanceAbsolute'] ):
                            edepWarnings.append( warning.energyImbalance( PQU.PQU(ein, totalEDep.axes[0].unit),
                                i, ein+Qsum, energyDepositedPerProduct(energyDep, ein), self ) )
                else:
                    # missing some products, so just check that outgoing energy doesn't exceed incoming:
                    for i,(ein,edep) in enumerate(totalEDep):
                        if ( (edep - (ein+Qsum)) > ((ein+Qsum) * info['dEnergyBalanceRelative'])
                                and (edep - (ein+Qsum)) > info['dEnergyBalanceAbsolute'] ):
                            edepWarnings.append( warning.energyImbalance( PQU.PQU(ein, totalEDep.axes[0].unit),
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

            try:
                Q = self.getQ('eV')
                checkProductsForEnergyBalance( products = [p1 for p1 in self.outputChannel], Qs = [Q],
                        flags = {'fission':self.outputChannel.isFission(), 'decay':False} )
            except ValueError:
                pass    # FIXME this test currently disabled when Q non-constant

        return warnings

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.crossSection.domainUnitConversionFactor( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.crossSection.domainUnit( ) )

    def getENDL_CS_ENDF_MT( self ) :
        """Returns the reaction's ENDL C, S and ENDF's MT values as integers in a python dictionary with keys 'C', 'S' and 'MT' 
        (e.g., { 'C' : 11, 'S' : 1, 'MT' : 53 })."""

        from fudge.legacy.converting import endf_endl

        MT = self.ENDF_MT
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

        return( self.documentation[name] )

    def getQ( self, unit, final = True, groundStateQ = False ) :
        """Returns the Q-value for this reaction. Converted to float if possible, otherwise a string value is returned."""

        return( self.Q.getConstantAs( unit ) )

    def setQ( self, Q ) :

        self.Q.setAncestor( self )
        self.Q = Q

    def getReactionSuite( self ) :

        from fudge.gnd import reactionSuite as reactionSuiteModule
        return( self.findClassInAncestry( reactionSuiteModule.reactionSuite ) )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        attributeString = ""
        if( self.EFL is not None ) : attributeString += ' EFL="%s"' % self.EFL
        if( self.date is not None ) : attributeString += ' date="%s"' % self.date
        if( self.process is not None ) : attributeString += ' process="%s"' % self.process
        attributeString += ' ENDF_MT="%s"' % self.ENDF_MT

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
# BRB ????????? self should always have an outputChannel
        if hasattr(self, 'outputChannel'):  # for reaction, fissionComponent and partialGammaProduction
            xmlString[-1] += ' outputChannel="%s"' % self.outputChannel
            fissionGenre = self.outputChannel.getFissionGenre( )
            if fissionGenre is not None: attributeString += ' fissionGenre="%s"' % fissionGenre
        if( hasattr( self, 'summands' ) ) : xmlString[-1] += ' name="%s"' % self.name

        xmlString[-1] += attributeString + '>'

        if self.documentation:
            xmlString.append( '%s<documentations>' % indent2 )
            for doc in self.documentation: xmlString += self.documentation[doc].toXMLList( indent3, **kwargs )
            xmlString[-1] += '</documentations>'

        if( hasattr( self, 'summands' ) ) : 
            for summand in self.summands : xmlString.append( summand.toXML( indent2, **kwargs ) )

        xmlString += self.crossSection.toXMLList( indent2, **kwargs )
        if hasattr(self, 'outputChannel'):
            xmlString += self.outputChannel.toXMLList( indent2, **kwargs )
        elif hasattr(self, 'Q'):
            xmlString += self.Q.toXMLList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

def isGNDReaction( o ) :
    """Returns True if o is an instance of base_reaction or of a subclass thereof. """

    return isinstance(o, base_reaction)
