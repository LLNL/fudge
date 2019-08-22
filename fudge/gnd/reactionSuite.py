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

"""
reactionSuite.py contains the 'reactionSuite' class that in turn holds all reactions for a given target/projectile.
reactionSuite is the top-level class for the GND structure.
"""

import re
import datetime

from fudge.gnd import suites as suitesModule
from fudge.gnd.reactionData import crossSection as crossSectionModule
from xData import standards as standardsModule
from fudge.particles import nuclear
from alias import aliases
from .version import GND_VERSION
import styles as stylesModule

import xData.ancestry as ancestryModule
import xData.link as linkModule

import fudge

__metaclass__ = type

class reactionSuite( ancestryModule.ancestry ) :
    """This is the main class for a gnd projectile/target object. It contains
        * documentation
        * a list of all particles encountered in the file
        * resonance data
        * a list of reactions """

    moniker = 'reactionSuite'

    def __init__( self, projectile, target, GND_version = GND_VERSION, style = None,
            documentation = None, particleList = None, projectileFrame = standardsModule.frames.labToken, MAT = None ) :
        """
        Creates a new reactionSuite object, containing all reaction data for a projectile/target combination.

        :param projectile: incident particle (instance of xParticle.base)
        :param target: target particle (instance of xParticle.base)
        :param GND_version: GND format version (string)
        :param style: Indicates style of data in this reactionSuite (evaluated, linearized, processed, etc.).  Instance of styles.style
        :param documentation: Top-level documentation for the entire reactionSuite (instance of documentation.documentation)
        :param particleList: list of particles appearing inside (instance of xParticleList.xParticleList)
        :param projectileFrame: frame of the projectile. Allowed values are 'lab' or 'centerOfMass'
        :param MAT: integer ENDF MAT number (only needed for writing back to ENDF-6)
        :return:
        """

        ancestryModule.ancestry.__init__( self )

        if( GND_version not in ( GND_VERSION ) ) : raise Exception( "Unsupported GND structure '%s'!" % str(GND_version) )

        if( not( isinstance( projectile, fudge.gnd.xParticle.base ) ) ) : projectile = fudge.particles.nuclear.getZAOrNameAs_xParticle( projectile )
        self.projectile = projectile
        self.projectileFrame = projectileFrame
        self.MAT = MAT

        if( not( isinstance( target, fudge.gnd.xParticle.base ) ) ) : target = fudge.particles.nuclear.getZAOrNameAs_xParticle( target )
        self.target = target

        self.GND_version = GND_version
        self.particles = particleList
        if( self.particles is None ) : self.particles = fudge.gnd.xParticleList.xParticleList()
        self.particles.setAncestor( self )
        self.attributes = {}

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )
        if( style is not None ) : self.styles.add( style )

        self.aliases = aliases( )
        self.aliases.setAncestor( self )

        self.resonances = None

        self.__reactions = suitesModule.reactions()
        self.__reactions.setAncestor( self )
        self.__sums = suitesModule.sums()
        self.__sums.setAncestor( self )
        self.__productions = suitesModule.productions()
        self.__productions.setAncestor( self )
        self.__fissionComponents = suitesModule.fissionComponents()
        self.__fissionComponents.setAncestor( self )

        self.partialGammaProductions = []   # FIXME should be removed?


        self.documentation = {}
        if( documentation is not None ) : self.addDocumentation( documentation )

    def __iter__( self ) :

        for reaction in self.reactions : yield reaction
        for reaction in self.partialGammaProductions : yield reaction
        for sum_ in self.sums : yield sum_
        for reaction in self.fissionComponents : yield reaction
        for reaction in self.productions : yield reaction

    def getProductionReactions( self ) :

        return( self.productions )

    def __str__( self ) :
        """See method toString."""

        return( self.toString( ) )

    @property
    def styles(self):
        return self.__styles

    @property
    def reactions(self):
        return self.__reactions

    @property
    def sums(self):
        return self.__sums

    @property
    def productions(self):
        return self.__productions

    @property
    def fissionComponents(self):
        return self.__fissionComponents

    def addAlias( self, key, value, attributes = {} ) :

        self.aliases.add( key, value, attributes )

    def addNuclearMetaStableAlias( self, isotopeName, nuclearLevelName, metaStableIndex ) :
        """Adds a nuclear meta stable alias to self's aliases."""

        self.aliases.addNuclearMetaStable( isotopeName, nuclearLevelName, metaStableIndex )

    def getAliasesFor( self, value ) :
        """Returns a list of all aliases that have value value."""

        return( self.aliases.getAliasesFor( value ) )

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def _addReactionHelper( self, reactionList, reaction, sortValue=None ) :
        """Helper function: add to self.productions, self.summedReactions, etc"""
        # FIXME: this is only used by partialGammaProduction, should go away once that is migrated

        reaction.setAncestor( self, attribute = 'label' )
        if sortValue is None: sortValue = len(reactionList)
        reaction.sortValue = sortValue
        for idx in xrange( len( reactionList ) ) :
            if( reactionList[idx].sortValue > reaction.sortValue ) :
                reactionList.insert( idx, reaction )
                return
        reactionList.append( reaction )

    def addPartialGammaProduction( self, gammaProduction, sortValue=None ) :
        """Adds a gnd.partialGammaProduction.partialGammaProduction instance to the list."""

        self._addReactionHelper( self.partialGammaProductions, gammaProduction, sortValue )

    def addParticle( self, particle ) :
        """Add particle to the list of particles where particle should be an instance of gnd.xParticle.base."""

        self.particles.addParticle( particle )

    def addResonances( self, resonances ):
        """Add resonance parameter data (gnd.resonances.resonances instance) to the reactionSuite. """
        self.resonances = resonances
        self.resonances.setAncestor( self )

    def check( self, **kwargs ) :
        """
        Check all data in the reactionSuite, returning a gnd.warning.context object with list of warnings.
        
        Currently supported options:
            'branchingRatioSumTolerance' 1e-6        # branching ratios must sum to 1 (within this tolerance)
            'dQ'                         '1e-3 MeV'  # tolerance for comparing tabulated / calculated Q values
            'dThreshold'                 '1e-3 MeV'
            'crossSectionEnergyMax'      '20 MeV'    # warn if cross section ends before this energy
            'crossSectionOnly':          False
            'crossSectionMaxDiff':       1e-3        # for comparing crossSectionSum to summands
            'multiplicityMaxDiff':       1e-3        # for comparing multiplicitySum to summands
            'transportables'             ('n',)      # distribution required for these products
            'normTolerance':             1e-5        # for checking distribution normalization
            'checkEnergyBalance'         True
            'reconstructResonances'      True
            'dEnergyBalanceRelative'     1e-3
            'dEnergyBalanceAbsolute'     1.0
            'fissionEnergyBalanceLimit'  0.15        # at least 85% of available energy should go to fission products
            'failOnException'            False       # if True, crash instead of converting exceptions to warnings
            'cleanUpTempFiles'           True        # remove derived data that was produced during checking
 
        Currently unused options:
            'checkForEnDepData'         False
            'allowZeroE'                False
            'xCloseEps'                 None
            'maxAbsFloatValue'          None
            'checkEMin'                 True
            'checkMissing'              True
            'maxGammaMultiplicity'      100.
          """

        from fudge.gnd import warning
        from fudge.processing import processingInfo as processingInfoModule

        options = {
                'branchingRatioSumTolerance': 1e-6,
                'dQ': '1e-3 MeV',
                'dThreshold': '1e-3 MeV',
                'crossSectionEnergyMax': '20 MeV',
                'crossSectionOnly': False,
                'crossSectionMaxDiff': 1e-3,
                'multiplicityMaxDiff': 1e-3,
                'transportables': ('n',),
                'normTolerance': 1e-5,
                'checkEnergyBalance': True,
                'reconstructResonances': True,
                'dEnergyBalanceRelative': 1e-3,
                'dEnergyBalanceAbsolute': 1.0,
                'fissionEnergyBalanceLimit': 0.15,
                'failOnException': False,
                'cleanUpTempFiles': True,
                # currently unused options:
                'checkForEnDepData': False,
                'allowZeroE': False,
                'xCloseEps': None,
                'maxAbsFloatValue': None,
                'checkEMin': True,
                'checkMissing': True,
                'maxGammaMultiplicity': 100.,
                }
        for key in kwargs:
            if key in options: options[key] = kwargs[key]
            else: raise KeyError, "check() received unknown keyword argument '%s'" % key

        warnings = []

        if isinstance( self.target, fudge.gnd.xParticle.FissionProduct ):
            warnings.append( warning.FissionProductsNotImplemented() )
            return warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        evaluatedStyle = self.styles.getEvaluatedStyle()
        if evaluatedStyle is None:
            warnings.append( warning.NotImplemented( "Checking currently only supported for 'evaluated' style" ) )
            return warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        # assemble some useful info, to be handed down to children's check() functions:
        compoundZA = self.target.getZ_A_SuffixAndZA()[3] + self.projectile.getZ_A_SuffixAndZA()[3]
        CoulombChannel = (self.target.getZ_A_SuffixAndZA()[0] != 0 and self.projectile.getZ_A_SuffixAndZA()[0] != 0)
        elementalTarget=False
        if hasattr(self.target,'getMass'):
            kinematicFactor = (self.target.getMass('amu')+self.projectile.getMass('amu')) / self.target.getMass('amu')  # (M+m)/M
            availableEnergy_eV = self.target.getMass('eV/c**2') + self.projectile.getMass('eV/c**2')
        else:
            # For elemental targets, calculating these factors doesn't make sense since there is no defined target mass
            kinematicFactor=1.0
            elementalTarget=True
            availableEnergy_eV=None

        info = { 'reactionSuite': self, 'kinematicFactor': kinematicFactor, 'compoundZA': compoundZA, 
                'availableEnergy_eV': availableEnergy_eV, 'CoulombChannel': CoulombChannel, 'style': evaluatedStyle,
                'reconstructedStyle': None, 'depositionStyle': None }
        info.update( options )
        if elementalTarget:
            # For elemental targets, calculating energy balance isn't possible
            info['checkEnergyBalance']=False

        if self.resonances is not None:
            resonanceWarnings = self.resonances.check( info )
            if resonanceWarnings:
                warnings.append( warning.context('resonances', resonanceWarnings) )
            if ( options['reconstructResonances'] and self.resonances.reconstructCrossSection ):
                info['reconstructedStyle'] = self.styles.getTempStyleNameOfClass( stylesModule.crossSectionReconstructed )
                # convert resonance parameters to pointwise data, interpolable to .1 percent:
                try: self.reconstructResonances( styleName=info['reconstructedStyle'], accuracy=0.001, thin=False, verbose=False )
                except Exception as e:
                    warnings.append( warning.ExceptionRaised( "when reconstructing resonances: %s" % e ) )
                    if info['failOnException']: raise

        if info['checkEnergyBalance']:
            # setup options for calculating energy deposition
            info['depositionStyle'] = self.styles.getTempStyleNameOfClass( stylesModule.averageProductData )
            depositionStyle = stylesModule.averageProductData( label=info['depositionStyle'], date=str(datetime.date.today()) )
            self.styles.add( depositionStyle )
            depositionStyle.derivedStyles.add(evaluatedStyle)
            processInfo = processingInfoModule.processInfo2( style=depositionStyle, reactionSuite=self,
                    verbosity=1 )   # suppress messages
            info['processInfo'] = processInfo
            info['tempInfo'] = { 'reactionSuite' : self }

        if self.projectile.name == 'n':
            # test Wick's limit: 0-degree elastic xsc >= ( total xsc * k/4pi )^2
            try:
                elastic = self.getReaction('elastic')
                total = self.getReaction('total')
            except Exception as e:
                warnings.append( warning.testSkipped("Wick's limit", e) )
            else:
                try :
                    if( info['reconstructedStyle'] in elastic.crossSection ) :
                        elastic_xsc = elastic.crossSection[info['reconstructedStyle']]
                        total_xsc = total.crossSection[info['reconstructedStyle']]
                    else:
                        elastic_xsc = elastic.crossSection.toPointwise_withLinearXYs()
                        total_xsc = total.crossSection.toPointwise_withLinearXYs()
                    elastic_distribution = elastic.outputChannel.getProductWithName('n').distribution[evaluatedStyle.label].angularSubform
                    if isinstance( elastic_distribution, fudge.gnd.productData.distributions.angular.pointwise ):
                        linearized = elastic_distribution.toPointwise_withLinearXYs()
                        forward_scattering = crossSectionModule.pointwise( axes=crossSectionModule.pointwise.defaultAxes() )
                        for energy_in in linearized.getEnergyArray('eV'):
                            forward_scattering.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                    elif isinstance( elastic_distribution, fudge.gnd.productData.distributions.angular.piecewise ):
                        forward_scattering = crossSectionModule.piecewise( axes=crossSectionModule.pointwise.defaultAxes() )
                        for region in elastic_distribution:
                            ptw = crossSectionModule.pointwise()
                            linearized = region.toPointwise_withLinearXYs()
                            for energy_in in linearized.getEnergyArray('eV'):
                                ptw.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                            forward_scattering.append( ptw )
                        forward_scattering = forward_scattering.toPointwise_withLinearXYs(1e-8,1e-8)

                    mutualDomain = zip( *[dat.domain() for dat in elastic_xsc, total_xsc, forward_scattering] )
                    mutualDomain = (max(mutualDomain[0]), min(mutualDomain[1]))
                    egrid = [e for e in forward_scattering.domainGrid() if mutualDomain[0] <= e <= mutualDomain[1]]

                    # get probability at mu=1.0:
                    forward_scattering = [forward_scattering.evaluate(e) for e in egrid]
                    elastic_xsc = [elastic_xsc.evaluate(e) for e in egrid]
                    total_xsc = [total_xsc.evaluate(e) for e in egrid]

                    wlcons = 3.05607e-8 * kinematicFactor**2 # ( sqrt(2 * neutronMass) / (4 * pi * hbar) * (M+m)/M )^2 in 1/(eV*b)
                    for i1 in range(len(egrid)):
                        if forward_scattering[i1] * elastic_xsc[i1] < wlcons * egrid[i1] * total_xsc[i1]**2:
                            ratio = (forward_scattering[i1] * elastic_xsc[i1]) / (wlcons * egrid[i1] * total_xsc[i1]**2 )
                            warnings.append( warning.WicksLimitError( 1-ratio, egrid[i1] ) )
                except Exception as e:
                    warnings.append( warning.ExceptionRaised( "when checking Wick's limit: %s" % e ) )
                    if info['failOnException']: raise

        particleWarnings = self.particles.check( info )
        if particleWarnings: warnings.append( warning.context('particles', particleWarnings) )

        for reaction in self :
            reactionWarnings = reaction.check( info )
            if reactionWarnings: warnings.append( warning.context('%s label %s: %s'
                % (reaction.moniker, reaction.getLabel(), reaction), reactionWarnings ) )

        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info

        if options['cleanUpTempFiles']:
            if info['reconstructedStyle'] is not None: self.removeStyle( info['reconstructedStyle'] )
            if info['depositionStyle'] is not None: self.removeStyle( info['depositionStyle'] )

        return result

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. reactionSuite contains several different lists,
        so may need to descend into those to find desired entity
        """
        if entityName in ('reaction','summedReaction','fissionComponent','production','reactionSum'):
            for entity in getattr( self, entityName+'s' ):
                if getattr( entity, attribute, None ) == value:
                    return entity
        return( ancestryModule.ancestry.findEntity( self, entityName, attribute, value ) )

    def hasAlias( self, key ) :

        return( key in self.aliases )

    def hasParticle( self, name ) :

        return( self.particles.hasParticle( name ) )

    def getProjectileFrame( self ) :

        return( self.projectileFrame )

    def getAlias( self, key ) :

        return( self.aliases[key] )

    def getDocumentation( self, name ) :

        return( self.documentation[name] )

    def getParticle( self, name ) :

        return( self.particles.getParticle( name ) )

    def getMassRatio( self ) :

        if not hasattr( self, '__massRatio' ):
            M = self.target.getMass('amu')
            m = self.projectile.getMass('amu')
            self.__massRatio = (M / (M+m))
        return self.__massRatio

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )

    def getIsotopeName( self, *args ):
        """ return name of compound nucleus formed by nuc1+nuc2+...
        if a nucleus in 'args' starts with '-', subtract instead """
        def parse(name):
            if 'gamma' in name: return 1,0,0,1
            sign,symbol,A,mult = re.search(r"([-]?)([a-zA-Z]+)(_natural|[0-9]*)(?:\[multiplicity:\')?([0-9]+)?", name).groups()
            if not mult: mult = 1
            if symbol=='n' and not A: A = 1
            if A!='_natural': A=int(A)
            return int(sign + '1'), nuclear.elementZFromSymbol(symbol), A, int(mult)
        retZ, retA = 0,0
        try:
            for nucleus in args:
                sign, Z, A, mult = parse(nucleus)
                retZ += sign*mult*Z
                if '_natural' in (A,retA): retA='_natural'
                else: retA += sign*mult*A
            return '%s%s' % (nuclear.elementSymbolFromZ(retZ), retA)
        except:
            print "      WARNING: couldn't extract isotope name from product list!"

    def getReaction( self, channel ):
        """
        Search list of reactions for a specified channel.
        The 'channel' argument should be either a reaction type ('elastic','capture','fission', etc) or
        the list of outgoing particles ('n + Pu239' for example).
        Raises 'KeyError' if specified channel can't be found.
        """

        # translate special channel names:
        if channel=='elastic': channel = channel_tr = '%s + %s' % (self.projectile,self.target)
        elif channel=='capture': channel_tr='z,gamma'
        else: channel_tr = channel

        # check if 'channel' == one of the fudge reactions:
        chStrings = []
        for reaction in self :
            if( str( reaction ) == channel ) : return( reaction )
            chStrings.append( str( reaction ) )
        
        # make list containing a set() of products for each reaction. Ignore energy-dependent multiplicity
        chSets = [set(aa.split(' + ')) for aa in [a.replace('(','').replace(')','').replace('->','+').replace(
            "[multiplicity:'energyDependent']",'') for a in chStrings] ]
        chSetsNoGamma = [s.copy() for s in chSets]; [a.discard('gamma') for a in chSetsNoGamma]
        if set(channel.split(' + ')) in chSets: return self.reactions[ chSets.index( set(channel.split(' + ')) ) ]
        if set(channel.split(' + ')) in chSetsNoGamma: return self.reactions[ chSetsNoGamma.index( set(channel.split(' + ')) ) ]

        if 'fission' in channel.lower():
            channel_fiss = channel.lower()
            def matchlist(*args):
                return any( [a in channel_fiss for a in args] )
            if channel_fiss=='fission' or 'total' in channel_fiss: genre='total'
            elif matchlist('first','1st'): genre='firstChance'
            elif matchlist('second','2nd'): genre='secondChance'
            elif matchlist('third','3rd'): genre='thirdChance'
            elif matchlist('fourth','4th'): genre='fourthChance'
            else:
                print "Can't determine fission genre from '%s'" % channel_fiss
                genre = None
            retVal  = [ r1 for r1 in self.reactions         if r1.outputChannel.fissionGenre == genre ]
            retVal += [ r1 for r1 in self.fissionComponents if r1.outputChannel.fissionGenre == genre ]
            if len(retVal)==1:
                return retVal[0]
        else:
            # check if channel is in form '(z,2na)', 'n,4n' or similar:
            patt = re.match('^[(]?[znpagdt],([a-zA-Z0-9]+)[)]?$', channel_tr)
            if patt:
                thisChannelSet = set()
                match = re.findall('([1-9]?)(gamma|He3|[npagdt]?)[+]?', patt.groups()[0] )
                for mul, prod in match:
                    if not prod : continue
                    prod = { 'g' : 'gamma' , 'gamma' : 'gamma', 'n' : 'n', 'p' : 'H1', 'd' : 'H2', 't' : 'H3', 'He3' : 'He3', 'a' : 'He4' }[prod]
                    if mul: prod += "[multiplicity:'%s']" % mul
                    if prod in thisChannelSet:
                        raise KeyError, "Please specify multiplicity explicitly ('z,2n' instead of 'z,nn')"
                    thisChannelSet.add(prod)
                if not thisChannelSet:
                    raise KeyError, "Channel '%s' could not be found!" % channel
                # also add final nucleus to the set:
                proj, target = str( self.projectile ), str( self.target )
                thisChannelSet.add( self.getIsotopeName( *([proj, target] + ['-'+a for a in thisChannelSet]) ) )
                #thisChannelSet.discard('gamma') # don't use gammas for comparison
                if thisChannelSet in chSets: return self.reactions[ chSets.index( thisChannelSet ) ]
                if thisChannelSet in chSetsNoGamma: return self.reactions[ chSetsNoGamma.index( thisChannelSet ) ]

        raise KeyError, "Channel '%s' could not be found!" % channel

    def getTemperature( self, style ) :

        return( self.styles[style].temperature )

    def calculateDepositionData( self, processInfo, verbosityIndent = '' ) :
        """
        Calculate average energy and momentum data for all products of all reactions.
        Resulting data information is stored within each product.
        """

        if( processInfo.verbosity >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        processInfo.checkStyle( stylesModule.averageProductData )
        self.styles.add( processInfo.style )
        tempInfo = { 'reactionSuite' : self }
        for reaction in self.reactions :
            reaction.calculateDepositionData( processInfo, tempInfo, verbosityIndent + processInfo.verbosityIndentStep )

    def inputParticlesToReactionString( self, prefix = "", suffix = "" ) :

        return( "%s%s + %s%s" % ( prefix, str( self.projectile ), str( self.target ), suffix ) )

    def process( self, processInfo, verbosityIndent = '' ) :

        from fudge.core.utilities import times

        if( processInfo.verbosity >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        t0 = times.times( )
        try :
            logFile = processInfo.dict['logFile']
        except :
            logFile = None

        tempInfo = { 'reactionSuite' : self }
        tempInfo['incidentEnergyUnit'] = self.reactions[0].crossSection.getIncidentEnergyUnit( )
        tempInfo['massUnit'] = tempInfo['incidentEnergyUnit'] + '/c**2'
        tempInfo['masses'] = { 'Projectile' : self.projectile.getMass( tempInfo['massUnit'] ) }
        tempInfo['masses']['Target'] = self.target.getMass( tempInfo['massUnit'] )
        tempInfo['masses']['Product'] = None
        tempInfo['masses']['Residual'] = None

        self.reconstructResonances( styleName = 'reconstructed', accuracy = 1e-3, verbose = False )
        processInfo['workFile'] = []
        for channel in self :
            channel.process( processInfo, tempInfo, verbosityIndent + processInfo.verbosityIndentStep )
        for channel in self.fissionComponents :
            channel.process( processInfo, tempInfo, verbosityIndent + processInfo.verbosityIndentStep )
        if( 'LLNL_Pn' in processInfo.dict['styles'] ) :
            subStyle = fudge.gnd.miscellaneous.subStyleLLNL_Pn( processInfo['particles'] )
            style = fudge.gnd.miscellaneous.style( 'LLNL_Pn', subStyle = subStyle )
            style.setAttribute( 'flux', processInfo['flux']['name'] )
            self.styles.add( style )
        if( logFile is not None ) : logFile.write( str( t0 ) + '\n' )

    def reconstructResonances( self, styleName, accuracy = None, thin = True, verbose = False ):
        """Turn resonance parameters into pointwise cross sections, then merge the results with
        tabulated pointwise cross sections. Resulting pointwise cross sections are stored
        alongside the original 'resonancesWithBackground' data in the reactionSuite.

        options:
            accuracy (double) - target accuracy during reconstruction. For example, 0.001
            thin (boolean) - enable/disable thinning after resonance reconstruction.
                Disabling thinning makes it easier to check for consistency of summed cross sections.
            verbose (boolean) - turn on/off verbosity"""

        from . import sums as sumsModule

        if( self.resonances is None ) : return
        if not self.resonances.reconstructCrossSection:
            return # nothing to do
        from fudge.processing.resonances import reconstructResonances

        xsecs = reconstructResonances.reconstructResonances(self, tolerance = accuracy, verbose = verbose)
        epsilon = 1e-8  # for joining multiple regions together

        evalStyle, = [style for style in self.styles if isinstance(style,stylesModule.evaluated)]
        newStyle = stylesModule.crossSectionReconstructed( label = styleName, date = str( datetime.date.today() ) )

        # for each reaction, add tabulated pointwise data (ENDF MF=3) to reconstructed resonances:
        possibleChannels = { 'elastic' : False, 'capture' : True, 'fission' : True, 'total' : False, 'nonelastic' : False }
        for reaction in self :
            if isinstance( reaction, sumsModule.multiplicitySum ): continue
            if not isinstance( reaction.crossSection[ evalStyle.label ], crossSectionModule.resonancesWithBackground ):
                continue
            # which reconstructed cross section corresponds to this reaction?
            RRxsec = None
            if str( reaction ) in xsecs:
                RRxsec = xsecs[ str( reaction ) ]
            else :
                for possibleChannel in possibleChannels :
                    if( possibleChannels[possibleChannel] ) :
                        if( possibleChannel in str( reaction ) ) : RRxsec = xsecs[possibleChannel]
                    if( RRxsec is None ) :
                        try :
                            if( reaction is self.getReaction( possibleChannel ) ) : RRxsec = xsecs[possibleChannel]
                        except :
                            pass
                    if( RRxsec is not None ) : break
            if( RRxsec is None ) :
                if verbose:
                    print( "Warning: couldn't find appropriate reconstructed cross section to add to reaction %s" % reaction )
                continue

            background = reaction.crossSection[evalStyle.label].tabulatedData
            background = background.toPointwise_withLinearXYs( epsilon, epsilon )
            RRxsec = RRxsec.toPointwise_withLinearXYs( epsilon, epsilon )

            background, RRxsec = background.mutualify(0,0,0, RRxsec, -epsilon,epsilon,True)
            RRxsec = background + RRxsec    # result is a crossSection.pointwise instance
            if RRxsec.rangeMin() < 0:
                # turn any negative xsc to 0
                RRxsec = RRxsec.clip( rangeMin=0 )
                if verbose:
                    print( "Warning: negative cross section encountered for %s; changed to 0 b" % reaction )
            if thin:
                RRxsec = RRxsec.thin( accuracy or .01 )
            RRxsec.label = newStyle.label
            reaction.crossSection.add( RRxsec )
        self.styles.add( newStyle )

    def removeStyle( self, style ) :
        """
        Remove the given style everywhere it appears in this reactionSuite
        :param style: may be a style instance or its name
        :return:
        """

        if isinstance(style, stylesModule.style): style = style.label

        def removeStyleFromComponent( styleName, component ):
            if styleName in component:
                component.remove( styleName )

        for reaction in self :
            # currently checks in cross section and deposition data, should also check multiplicity and distributions
            if hasattr(reaction, 'crossSection'):
                removeStyleFromComponent( style, reaction.crossSection )

            if not hasattr(reaction, 'outputChannel'): continue
            for product in reaction.outputChannel:
                removeStyleFromComponent( style, product.energyDeposition )
                removeStyleFromComponent( style, product.momentumDeposition )
                if product.decayChannel is None: continue
                for dproduct in product.decayChannel:
                    removeStyleFromComponent( style, dproduct.energyDeposition )
                    removeStyleFromComponent( style, dproduct.momentumDeposition )
        self.styles.remove( style )

    def saveToOpenedFile( self, fOut, **kwargs ) :

        xmlString = self.toXMLList( **kwargs )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, **kwargs ) :
        """Save the reactionSuite in GND/xml structure to specified file.
        To suppress extra messages, change the 'verbosity' flag:
            >>>self.saveToFile( "output.xml", flags={'verbosity':0} )
        """

        with open( fileName, "w" ) as fout :
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n' )
            self.saveToOpenedFile( fout, **kwargs )

    def toXMLList( self, indent = "", **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        attributeString = ''
        for attribute in self.attributes : attributeString += ' %s="%s"' % ( attribute, self.attributes[attribute] )
        xmlString = [ '%s<%s projectile="%s" target="%s" version="%s"%s xmlns:xlink="http://www.w3.org/1999/xlink" projectileFrame="%s">'
            % ( indent, self.moniker, self.projectile.name, self.target.name, self.GND_version, attributeString, self.projectileFrame ) ]

        xmlString += self.styles.toXMLList( indent2, **kwargs )

        xmlString.append( '%s<documentations>' % indent2 )
        for doc in self.documentation : xmlString += self.documentation[doc].toXMLList( indent3, **kwargs )
        xmlString[-1] += '</documentations>'

        xmlString += self.aliases.toXMLList( indent2, **kwargs )

        xmlString += self.particles.toXMLList( indent2, **kwargs )
        
        if self.resonances is not None:
            xmlString += self.resonances.toXMLList( indent2, **kwargs )

        xmlString += self.reactions.toXMLList( indent2, **kwargs )
        xmlString += self.sums.toXMLList( indent2, **kwargs )
        xmlString += self.fissionComponents.toXMLList( indent2, **kwargs )
        xmlString += self.productions.toXMLList( indent2, **kwargs )
        # FIXME: remove partialGammaProductions
        for reaction in self.partialGammaProductions : xmlString += reaction.toXMLList( indent2, **kwargs )

        xmlString.append( '%s</%s>' % ( indent, self.moniker ) )
        return( xmlString )

    def setAttribute( self, name, value ) :
        """Adds attribute name and its value to the list of attributes."""

        self.attributes[name] = value

    def toString( self, indent = '' ) :
        """Returns a string representation of an reactionSuite."""

        s = indent + self.inputParticlesToReactionString( suffix = ' --> ' )
        indent2 = len( s ) * ' '
        indent3 = ''
        for channel in self.reactions :
            s += channel.toString( indent = indent3 )
            indent3 = indent2
        return( s )

def readXML( gndFile ):
    """
    Read a GND/xml file and create a new reactionSuite instance from the result.

    :param gndFile: path to a GND file, as a string.
    :return: reactionSuite instance containing all data from the file.
    """
    from xml.etree import cElementTree
    rsElement = cElementTree.parse( gndFile ).getroot()
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode
    rsElement = xmlNode( rsElement, xmlNode.etree )
    return parseXMLNode( rsElement )

def parseXMLNode( rsElement ):
    """Translates a <reactionSuite> xml node into a reactionSuite instance. Users should use the 'readXML' function instead."""

    xPath = ['reactionSuite']   # Keep track of location in the tree, in case errors come up.
    try :
        GND_version = rsElement.get( 'version' )
        if GND_version is None:
            GND_version = rsElement.get( 'format' )

        particles = fudge.gnd.xParticleList.parseXMLNode( rsElement.find( 'particles' ) )
        projectile = particles.getParticle( rsElement.get( 'projectile' ) )
        target = particles.getParticle( rsElement.get('target') )
        projectileFrame = rsElement.get( 'projectileFrame' )
#        formatVersion = rsElement.get( 'formatVersion' ) #FIXME: obsolete? DAB 11/24/2015
        rs = reactionSuite( projectile, target, GND_version = GND_version, particleList = particles,
                            projectileFrame = projectileFrame )

        linkData = { 'particles' : particles, 'reactionSuite' : rs, 'unresolvedLinks' : [], 'format' : GND_version }
        rs.styles.parseXMLNode( rsElement.find('styles'), xPath, linkData )
        for doc in rsElement.find( 'documentations' ) :
            rs.addDocumentation( fudge.gnd.documentation.documentation.parseXMLNode( doc, xPath, linkData ) )

        for child in rsElement :
            if( child.tag in ( 'styles', 'documentations', 'particles' ) ) :
                continue    # already read above
            elif( child.tag == 'aliases' ) :
                for alias in child :
                    attributes = {}
                    for key in alias.keys( ) :
                        if( key in [ 'key', 'value' ] ) : continue
                        attributes[key] = alias.get( key )
                    rs.addAlias( alias.get( 'key' ), alias.get( 'value' ), attributes = attributes )
            elif( child.tag == 'resonances' ) :
                rs.addResonances( fudge.gnd.resonances.resonances.parseXMLNode( child, xPath, linkData ) )
            elif( child.tag == 'reactions' ) :
                rs.reactions.parseXMLNode( child, xPath, linkData )
            elif( child.tag == 'sums' ) :
                rs.sums.parseXMLNode( child, xPath, linkData )
            elif( child.tag == 'productions' ) :
                rs.productions.parseXMLNode( child, xPath, linkData)
            elif( child.tag == 'fissionComponents' ) :
                rs.fissionComponents.parseXMLNode( child, xPath, linkData )
            # FIXME: remove partialGammaProduction
            elif( child.tag == 'partialGammaProduction' ) :
                rs.addPartialGammaProduction( fudge.gnd.partialGammaProduction.parseXMLNode( child, xPath, linkData ) )
            else :
                print( "Warning: encountered unexpected element '%s' in reactionSuite!" % child.tag )
    except Exception :
        print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
        raise

    # Fix links.
    for quant in linkData['unresolvedLinks'] :
        if isinstance( quant, linkModule.link ) :
            quant.link = quant.follow( rs )
        else :
            raise 'hell'

    return( rs )
