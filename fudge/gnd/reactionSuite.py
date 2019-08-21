# <<BEGIN-copyright>>
# <<END-copyright>>

"""
reactionSuite.py contains the 'reactionSuite' class that in turn holds all reactions for a given target/projectile.
reactionSuite is the top-level class for the GND format
"""

formatMajorVersion = 1
formatMinorVersion = 2
formatVersion = "gnd version %s.%s" % ( formatMajorVersion, formatMinorVersion )

import re
from fudge.core.ancestry import ancestry
from fudge.core.utilities import fudgeExceptions, fudgeZA, brb
from pqu import physicalQuantityWithUncertainty
from fudge.legacy.converting import endfFormats, endf_endl
from alias import aliases

from fudge.processing import processingInfo
from fudge.processing.deterministic import gndToLLNLSnCOut
import fudge
from fudge.core.math.xData import XYs

__metaclass__ = type

monikerReactionSuite = 'reactionSuite'

class reactionSuite( ancestry ) :
    """This is the main class for a gnd projectile/target object. It contains
        * documentation
        * a list of all particles encountered in the file
        * resonance data
        * a list of reactions """

    def __init__( self, projectile, target, style = None, documentation = None, particleList = None ) :
        """Creates a new reactionSuite object. 'projectile' and 'target' should both be
        gnd.xParticle.xParticle or gnd.xParticle.nuclearLevel instances."""

        ancestry.__init__( self, monikerReactionSuite, None )

        if( not( isinstance( projectile, fudge.gnd.xParticle.xParticle ) ) ) : projectile = fudge.particles.nuclear.getZAOrNameAs_xParticle( projectile )
        self.projectile = projectile

        if( not( isinstance( target, fudge.gnd.xParticle.xParticle ) ) ) :
            if( not( isinstance( target, fudge.gnd.xParticle.nuclearLevel ) ) ) : target = fudge.particles.nuclear.getZAOrNameAs_xParticle( target )
        self.target = target

        self.format = formatVersion
        self.particles = particleList
        if( self.particles is None ) : self.particles = fudge.gnd.xParticleList.xParticleList()
        self.particles.setParent( self )
        self.attributes = {}
        self.aliases = aliases( )
        self.reactions = []
        self.partialGammaProductions = []
        self.summedReactions = []
        self.productions = []
        self.fissionComponents = []
        self.styles = {}
        if( style is not None ) : self.addStyle( style )
        self.documentation = {}
        if( documentation is not None ) : self.addDocumentation( documentation )
        self.referredData = fudge.gnd.referredData.referredData( parent = self )

    def __getitem__( self, i ) :
        """Returns the (i-1)^th reaction channel."""

        return( self.reactions[i] )

    def __len__( self ) :
        """Returns the number of reaction channels in the reactionSuite."""

        return( len( self.reactions ) )
    
    def _getBaseAndDerivedReactions( self ):
        """ get list of all possible reactions (including discrete, summed and production reactions ) """

        return (self.reactions + self.partialGammaProductions + self.summedReactions +
                self.fissionComponents + self.productions)

    def __str__( self ) :
        """See method toString."""

        return( self.toString( ) )

    def addAlias( self, key, value ) :

        self.aliases[key] = value

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def _addReactionHelper( self, reactionList, reaction, sortValue=None ) :
        """Helper function: add to self.reactions, self.productions, self.summedReactions, etc"""

        reaction.setParent( self )
        if sortValue is None: sortValue = len(reactionList)
        reaction.sortValue = sortValue
        for i in xrange( len( reactionList ) ) :
            if( reactionList[i].sortValue > reaction.sortValue ) :
                reactionList.insert( i, reaction )
                return
        reactionList.append( reaction )

    def addReaction( self, reaction, sortValue=None ) :
        """Adds a reaction (gnd.reaction.reaction instance) to the list of reactions."""

        self._addReactionHelper( self.reactions, reaction, sortValue )

    def addPartialGammaProduction( self, gammaProduction, sortValue=None ) :
        """Adds a gnd.partialGammaProduction.partialGammaProduction instance to the list."""

        self._addReactionHelper( self.partialGammaProductions, gammaProduction, sortValue )

    def addSummedReaction( self, summedReaction, sortValue=None ) :
        """Adds a summed reaction (gnd.summedReaction.summedReaction instance) to the list."""

        self._addReactionHelper( self.summedReactions, summedReaction, sortValue )

    def addFissionComponent( self, fissionComponent, sortValue=None ) :
        """Adds a fission component reaction (gnd.fissionComponent.fissionComponent instance) to the list."""

        self._addReactionHelper( self.fissionComponents, fissionComponent, sortValue )

    def addProductionReaction( self, productionReaction, sortValue=None ):
        """Adds a production reaction (gnd.production.production instance) to the list."""

        self._addReactionHelper( self.productions, productionReaction, sortValue )

    def addParticle( self, particle ) :
        """Add either a particle or an excited level to the list of particles.
        'particle' should be either gnd.xParticle.xParticle or gnd.xParticle.nuclearLevel """

        self.particles.addParticle( particle )

    def addResonances( self, resonances ):
        """Add resonance parameter data (gnd.resonances.resonances instance) to the reactionSuite. """
        self.resonances = resonances
        self.resonances.parent = self
    
    def addReferredData( self, referredData ) :
        """Add a referredData section (used when data from two or more reactions are given as energy-dependent
        multiples of a 'reference' data). Currently only used when translating ENDF MF=9 data, and may be
        deprecated in the future.
        See gnd.referredData module for more info. """

        return( self.referredData.appendDatum( referredData ) )

    def addStyle( self, style ) :
        """Add document style, including information on library, version and type of data.
        Data type is usually 'evaluated', but 'processed' and 'experimental' are also possible data types.
        See gnd.miscellaneous.style for more info. """

        self.styles[style.name] = style

    def check( self, **kwargs ) :
        """
        Check all data in the reactionSuite, returning a list of warnings. 
        
        Currently supported options:
            'dQ'                        '1e-3 MeV'
            'dThreshold'                '1e-3 MeV'
            'crossSectionEnergyMax'     '20 MeV'
            'crossSectionOnly':         False
            'transportables'            ('n',)      # distribution required for these products
            'normTolerance':            1e-5
            'checkEnergyBalance'        True
            'reconstructResonances'     True
            'dEnergyBalanceRelative'    1e-3
            'dEnergyBalanceAbsolute'    1.0
            'fissionEnergyBalanceLimit' 0.15        # at least 85% should go to fission products
 
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

        options = {
                'dQ': '1e-3 MeV',
                'dThreshold': '1e-3 MeV',
                'crossSectionEnergyMax': '20 MeV',
                'crossSectionOnly': False,
                'transportables': ('n',),    # distribution required for these products
                'normTolerance': 1e-5,
                'checkEnergyBalance': True,
                'reconstructResonances': True,
                'dEnergyBalanceRelative': 1e-3,
                'dEnergyBalanceAbsolute': 1.0,
                'fissionEnergyBalanceLimit': 0.15,  # at least 85% should go to fission products
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

        # assemble some useful info, to be handed down to children's check() functions:
        kinematicFactor = (self.target.getMass('amu')+self.projectile.getMass('amu')) / self.target.getMass('amu')  # (M+m)/M
        compoundZA = self.target.getZ_A_SuffixAndZA()[3] + self.projectile.getZ_A_SuffixAndZA()[3]
        availableEnergy_eV = self.target.getMass('eV/c**2') + self.projectile.getMass('eV/c**2')
        info = { 'reactionSuite': self, 'kinematicFactor': kinematicFactor, 'compoundZA': compoundZA, 
                'availableEnergy_eV': availableEnergy_eV, }
        info.update( options )

        warnings = []

        if ( options['reconstructResonances'] and hasattr(self, 'resonances')
                and self.resonances.reconstructCrossSection ):
            # convert resonance parameters to pointwise data, interpolable to .1 percent:
            self.reconstructResonances( verbose=False )

        if info['checkEnergyBalance']:
            # setup options for calculating energy deposition
            info['processInfo'] = processingInfo.processInfo( self, verbosity=1 )   # suppress messages

        if self.projectile.name == 'n':
            # test Wick's limit: 0-degree elastic xsc >= ( total xsc * k/4pi )^2
            try:
                elastic_distribution = ( self.getReaction('elastic').outputChannel.getProductWithName('n')
                    .distributions.toPointwiseLinear() )
                egrid = elastic_distribution.getEnergyArray('eV')
                forward_scattering = [e.getValue(1.0) for e in elastic_distribution]    # probability at mu=1.0
                elastic_xsc = self.getReaction('elastic').crossSection.toPointwiseLinear()
                elastic_xsc = [elastic_xsc.getValue(e) for e in egrid]
                total_xsc = self.getReaction('total').crossSection.toPointwiseLinear()
                total_xsc = [total_xsc.getValue(e) for e in egrid]

                wlcons = 3.05607e-8 * kinematicFactor**2 # ( sqrt(2 * neutronMass) / (4 * pi * hbar) * (M+m)/M )^2 in 1/(eV*b)
                for i in range(len(egrid)):
                    if forward_scattering[i] * elastic_xsc[i] < wlcons * egrid[i] * total_xsc[i]**2:
                        ratio = (forward_scattering[i] * elastic_xsc[i]) / (wlcons * egrid[i] * total_xsc[i]**2 )
                        warnings.append( warning.WicksLimitError( 1-ratio, egrid[i] ) )
            except Exception as e:
                warnings.append( warning.ExceptionRaised( "when checking Wick's limit: %s" % e ) )

        particleWarnings = self.particles.check( info )
        if particleWarnings: warnings.append( warning.context('particles', particleWarnings) )

        for reaction in self._getBaseAndDerivedReactions():
            reactionWarnings = reaction.check( info )
            if reactionWarnings: warnings.append( warning.context('%s label %s: %s'
                % (reaction.moniker, reaction.getLabel(), reaction), reactionWarnings ) )

        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info
        return result

    def hasParticle( self, name ) :

        return( self.particles.hasParticle( name ) )

    def getDocumentation( self, name ) :

        return( self.documentation[doc.name] )

    def getParticle( self, name ) :

        return( self.particles.getParticle( name ) )

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )
    
    def getIsotopeName( self, *args ):
        """ return name of compound nucleus formed by nuc1+nuc2+...
        if a nucleus in 'args' starts with '-', subtract instead """
        def parse(name):
            if 'gamma' in name: return 1,0,0,1
            sign,Z,A,mult = re.search(r"([-]?)([a-zA-Z]+)(_natural|[0-9]*)(?:\[multiplicity:\')?([0-9]+)?", name).groups()
            if not mult: mult = 1
            if Z=='n' and not A: A = 1
            if A!='_natural': A=int(A)
            return int(sign + '1'), fudgeZA.SymbolToZ(Z), A, int(mult)
        retZ, retA = 0,0
        try:
            for nucleus in args:
                sign, Z, A, mult = parse(nucleus)
                retZ += sign*mult*Z
                if '_natural' in (A,retA): retA='_natural'
                else: retA += sign*mult*A
            return '%s%s' % (fudgeZA.ZToSymbol(retZ), retA)
        except:
            print "      WARNING: couldn't extract isotope name from product list!"

    def getReaction( self, channel ):
        """Search list of reactions for a specified channel.
        The 'channel' argument should be either a reaction type ('elastic','capture','fission', etc) or
        the list of outgoing particles ('n + Pu239' for example).
        Raises 'KeyError' if specified channel can't be found. """

        # translate special channel names:
        if channel=='elastic': channel='%s + %s' % (self.projectile,self.target)
        elif channel=='capture': channel_tr='z,gamma'
        else: channel_tr = channel

        # check if 'channel' == one of the fudge reactions:
        allChannels = self._getBaseAndDerivedReactions()
        chStrings = [str(a) for a in allChannels]
        if channel in chStrings: return allChannels[ chStrings.index( channel ) ]
        
        # make list containing a set() of products for each reaction. Ignore energy-dependent multiplicity
        chSets = [set(aa.split(' + ')) for aa in [a.replace('(','').replace(')','').replace('->','+').replace(
            "[multiplicity:'energyDependent']",'') for a in chStrings] ]
        chSetsNoGamma = [s.copy() for s in chSets]; [a.discard('gamma') for a in chSetsNoGamma]
        if set(channel.split(' + ')) in chSets: return self.reactions[ chSets.index( set(channel.split(' + ')) ) ]
        if set(channel.split(' + ')) in chSetsNoGamma: return self.reactions[ chSetsNoGamma.index( set(channel.split(' + ')) ) ]
        
        # if channel is in form '(z,2na)', 'n,4n' or similar:
        patt = re.match('^[(]?[zn],([a-zA-Z0-9]+)[)]?$', channel_tr)
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
        
        elif 'fission' in channel.lower():
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
            retVal = [a for a in self.reactions + self.fissionComponents if a.outputChannel.fissionGenre==genre]
            if len(retVal)==1:
                return retVal[0]
        raise KeyError, "Channel '%s' could not be found!" % channel

    def getReferredData( self, key ) :

        return( self.referredData[key] )

    def getTemperature( self ) :

        return( self.attributes['temperature'] )

    def calculateDepositionData( self, processInfo, verbosityIndent = '' ) :
        """Calculate energy and momentum deposition for all reactions and all products.
        Resulting deposition information is stored within each product. """

        if( processInfo['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        for channel in self : channel.calculateDepositionData( processInfo, verbosityIndent + '    ' )

    def inputParticlesToReactionString( self, prefix = "", suffix = "" ) :

        return( "%s%s + %s%s" % ( prefix, str( self.projectile ), str( self.target ), suffix ) )

    def process( self, processInfo, verbosityIndent = '' ) :

        from fudge.core.utilities import times
        t0 = times.times( )
        try :
            logFile = processInfo.dict['logFile']
        except :
            logFile = None
        tempInfo = processingInfo.tempInfo( )
        tempInfo['reactionSuite'] = self
        if( processInfo['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        self.reconstructResonances( accuracy = 1e-3, verbose = False )
        for channel in self : channel.process( processInfo, tempInfo, verbosityIndent + '    ' )
        if( 'LLNL_MC' in processInfo.dict['styles'] ) :
            style = fudge.gnd.miscellaneous.style( 'LLNL_MC' )
            self.addStyle( style )
        if( 'LLNL_Pn' in processInfo.dict['styles'] ) :
            subStyle = fudge.gnd.miscellaneous.subStyleLLNL_Pn( processInfo['particles'] )
            style = fudge.gnd.miscellaneous.style( 'LLNL_Pn', subStyle = subStyle )
            style.setAttribute( 'flux', processInfo['flux']['name'] )
            self.addStyle( style )
        if( logFile is not None ) : logFile.write( str( t0 ) + '\n' )
    
    def reconstructResonances( self, accuracy = None, verbose = False ):
        """Turn resonance parameters into pointwise cross sections, then merge the results with
        tabulated pointwise cross sections. Resulting pointwise cross sections are stored
        alongside the original 'resonancesWithBackground' data in the reactionSuite. """

        if( not( hasattr( self, 'resonances' ) ) ) : return
        if not self.resonances.reconstructCrossSection:
            return # nothing to do
        from fudge.processing.resonances import fudgeReconstructResonances
        from fudge.core.math.xData import axes

        xsecs = fudgeReconstructResonances.reconstructResonances(self, accuracy, verbose)
        epsilon = 1e-8  # for joining multiple regions together

        # for each channel, add tabulated pointwise data (ENDF MF=3) to reconstructed resonances:
        needReconstruction = [reac for reac in self._getBaseAndDerivedReactions()
                if reac.crossSection.nativeData == fudge.gnd.tokens.resonancesWithBackgroundFormToken]
        for channel in needReconstruction:
            # which reconstructed cross section corresponds to this channel?
            if str(channel) in xsecs: RRxsec = xsecs[ str(channel) ]
            elif channel is self.getReaction('elastic'): RRxsec = xsecs[ 'elastic' ]
            elif channel is self.getReaction('capture') or 'capture' in str(channel):
                RRxsec = xsecs[ 'capture' ]
            elif channel is self.getReaction('total'):
                RRxsec = xsecs[ 'total' ]
            elif 'fission' in str(channel):
                RRxsec = xsecs[ 'fission' ]
            else:
                raise Exception, "Couldn't find appropriate reconstructed cross section to add to channel %s"%channel

            background = channel.crossSection.forms[ channel.crossSection.nativeData ].tabulatedData
            background = background.toPointwiseLinear( epsilon, epsilon )
            RRxsec = RRxsec.toPointwiseLinear( epsilon, epsilon )
            # before adding to another pointwise region, y-value of upper point must be 0:
            RRdomain, pwDomain = RRxsec.getDomain(), background.getDomain()
            if RRdomain[0] < pwDomain[0] or RRdomain[1] > pwDomain[1]:
                raise Exception("Resonance region domain must fit inside background domain!")
            if pwDomain[0] < RRdomain[0]:
                RRxsec = RRxsec.dullEdges( lowerEps=-epsilon, positiveXOnly=True )
            if pwDomain[1] > RRdomain[1]:
                RRxsec = RRxsec.dullEdges( upperEps=epsilon, positiveXOnly=True )
            RRxsec = background + RRxsec    # result is a crossSection.pointwise instance
            if RRxsec.yMin() < 0:
                # turn any negative xsc to 0
                RRxsec = RRxsec.applyFunction( lambda y,dummy: 0 if y<0 else y, None )
                if verbose:
                    print( "Warning: negative cross section encountered for %s; changed to 0 b" % channel )
            RRxsec = RRxsec.thin( accuracy or .01 )
            channel.crossSection.addForm( RRxsec )

    def removeReaction( self, reactionName ):
        del self.reactions[ self.reactions.index( self.getReaction( reactionName ) ) ]

    def removeDerivedData( self ):
        """ Clean up an evaluation by removing everything except the original 'nativeData'
        that was included by the evaluator. """
        for reaction in self._getBaseAndDerivedReactions():
            xscNative = reaction.crossSection.nativeData
            reaction.crossSection.forms = {xscNative: reaction.crossSection.forms[ xscNative ]}
            for product in getattr( reaction, 'outputChannel', [] ):
                product.data = {}
                if product.decayChannel is not None:
                    for decayProd in product.decayChannel:
                        decayProd.data = {}

    def saveToOpenedFile( self, fOut, flags = None, verbosityIndent = '' ) :
        if not flags:
            flags = processingInfo.tempInfo()
            flags['verbosity'] = 31

        xmlString = self.toXMLList( flags, verbosityIndent = verbosityIndent )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, flags = None, verbosityIndent = '' ) :
        """Save the reactionSuite in GND/xml format to specified file.
        To suppress extra messages, change the 'verbosity' flag:
            >>>self.saveToFile( "output.xml", flags={'verbosity':0} ) 
        """

        with open(fileName,"w") as fout:
            self.saveToOpenedFile( fout, flags, verbosityIndent )

    def toLLNLSnCOut( self, toOtherData, first = False, last = False, verbosityIndent = '' ) :

        tempInfo = processingInfo.tempInfo( )
        tempInfo['reactionSuite'] = self
        if( toOtherData['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        return( gndToLLNLSnCOut.toLLNLSnCOut( self, toOtherData, tempInfo, first, last, verbosityIndent ) )

    def toXMLList( self, flags, verbosityIndent = '' ) :

        attributeString = ""
        for attribute in self.attributes : attributeString += ' %s="%s"' % ( attribute, self.attributes[attribute] )
        xmlString = [ '<?xml version="1.0" encoding="UTF-8"?>' ]
        xmlString.append( ( '<%s projectile="%s" target="%s" format="%s"%s xmlns:xlink="http://www.w3.org/1999/xlink">' \
            % ( self.moniker, self.projectile.getName( ), self.target.getName( ), self.format, attributeString ) ) )

        xmlString.append( '  <styles>' )
        for style in self.styles : xmlString += self.styles[style].toXMLList( indent = '    ' )
        xmlString[-1] += '</styles>'

        xmlString.append( '  <documentations>' )
        for doc in self.documentation : xmlString += self.documentation[doc].toXMLList( indent = '    ' )
        xmlString[-1] += '</documentations>'

        xmlString += self.aliases.toXMLList( indent = '  ' )
        xmlString.append( '  <particles>' )
        particles = sorted( self.particles.values( ) )
        for particle in particles: xmlString += particle.toXMLList( indent = '    ' )
        xmlString[-1] += '</particles>'
        
        if hasattr(self,'resonances'):
            xmlString += self.resonances.toXMLList( flags, verbosityIndent + '    ', indent = '  ' )

        if( flags['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        for reaction in self._getBaseAndDerivedReactions():
            xmlString += reaction.toXMLList( flags, verbosityIndent + '    ', indent = '  ' )
        xmlString += self.referredData.toXMLList( indent='  ' )

        xmlString.append( '</%s>' % self.moniker)
        return( xmlString )

    def toENDF6( self, flags, verbosityIndent = '', covarianceSuite = None ) :

        from fudge.legacy.converting import gndToENDF6

        if( flags['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
        verbosityIndent2 = verbosityIndent + ' ' * ( len( self.inputParticlesToReactionString( suffix = " -->" ) ) + 1 )
        projectile, target = self.projectile, self.target
        projectileZA = projectile.getZ_A_SuffixAndZA( )[-1]
        targetZA, MAT = endf_endl.ZAAndMATFromParticleName( target.getName( ) )
        targetZ, targetA = target.getZ_A_SuffixAndZA( )[:2]
        targetInfo = processingInfo.tempInfo( )
        targetInfo['reactionSuite'] = self
        targetInfo['ZA'] = targetZA
        # need neutron mass in eV/c**2, but it may not be in the particle list:
        if( 'n' in self.particles ) :
            targetInfo['neutronMass'] = self.getParticle( 'n' ).getMass( 'eV/c**2' )
        else:
            from fudge.structure import masses
            neutronAmu = masses.getMassFromZA( 1 )
            targetInfo['neutronMass'] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( neutronAmu, 'amu' ).getValueAs('eV/c**2')
        targetInfo['mass'] = target.getMass( 'eV/c**2' ) / targetInfo['neutronMass']
        try :
            targetInfo['LIS'] = target['levelIndex']
        except :
            targetInfo['LIS'] = 0
        targetInfo['LISO'] = 0
        for key in self.aliases.aliases :
            if( self.aliases[key] == target.getName( ) ) :
                targetInfo['LISO'] = int( key.split( '_m' )[1] )
                break
        MAT += targetInfo['LISO']
        targetInfo['delayedRates'] = []
        targetInfo['totalDelayedNubar'] = None
        targetInfo['MTs'], targetInfo['MF8'], targetInfo['LRs'] = {}, {}, {}
        endfMFList = { 1 : { 451 : [] }, 2 : {}, 3 : {}, 4 : {}, 5 : {}, 6 : {}, 8 : {}, 9 : {}, 10 : {}, 12 : {}, 13 : {}, 14 : {}, 15 : {},
                31 : {}, 32 : {}, 33 : {}, 34 : {}, 35 : {}, 40 : {} }
        if hasattr(self,'resonances'):      # Add resonances, independent of reaction channels
            self.resonances.toENDF6( endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent2 )
        for channel in self._getBaseAndDerivedReactions() :
            channel.toENDF6( endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent2 )
        gndToENDF6.upDateENDFMF8Data( endfMFList, targetInfo )
    
        # gamma decay data:
        for particle in self.particles.values():
            for key in particle.levels:
                level = particle.levels[key]
                if level.gammas: # non-empty gamma information
                    for baseMT in [ 50, 600, 650, 700, 750, 800 ] :
                        residualZA = endf_endl.ENDF_MTZAEquation( projectileZA, targetZA, baseMT )[0][-1]
                        if( fudgeZA.ZAToGNDName( residualZA ) == particle.name ) : break
                    level.toENDF6( baseMT, endfMFList, flags, targetInfo )

        MFs = sorted( endfMFList.keys( ) )
        directory = []
        endfList = []

        totalNubar = None
        totalDelayedNubar = targetInfo['totalDelayedNubar']
        if( 455 in endfMFList[5] ) :
            MF5MT455s = endfMFList[5][455]

            endfMFList[1][455]  = [ endfFormats.endfHeadLine( targetZA, targetInfo['mass'], 0, 2, 0, 0 ) ] # Currently, only LDG = 0, LNU = 2 is supported.
            endfMFList[1][455] += [ endfFormats.endfHeadLine( 0, 0, 0, 0, len( targetInfo['delayedRates'] ), 0 ) ]
            endfMFList[1][455] += endfFormats.endfDataList( targetInfo['delayedRates'] )
            n = len( totalDelayedNubar )
            endfMFList[1][455] += [ endfFormats.endfHeadLine( 0, 0, 0, 0, 1, n ) ]
            endfMFList[1][455] += endfFormats.endfInterpolationList( [ n, 2 ] )
            data = []
            for xy in totalDelayedNubar.copyDataToXYs( ) : data += xy
            endfMFList[1][455] += endfFormats.endfDataList( data ) + [ endfFormats.endfSENDLineNumber( ) ]

            MF5MT455List = [ endfFormats.endfHeadLine( targetZA, targetInfo['mass'], 0, 0, len( MF5MT455s ), 0 ) ]
            for MF5MT455 in MF5MT455s : MF5MT455List += MF5MT455
            if( len( MF5MT455s ) == 0 ) :
                del endfMFList[5][455]
            else :
                endfMFList[5][455] = MF5MT455List + [ endfFormats.endfSENDLineNumber( ) ]
        if(   'promptNubar' in targetInfo.dict ) :
            endfMFList[1][456]  = [ endfFormats.endfHeadLine( targetZA, targetInfo['mass'], 0, 2, 0, 0 ) ]       # Currently, LNU = 2 is supported.
            promptNubar = targetInfo['promptNubar']
            n = len( promptNubar )
            endfMFList[1][456] += [ endfFormats.endfHeadLine( 0, 0, 0, 0, 1, n ) ]
            endfMFList[1][456] += endfFormats.endfInterpolationList( [ n, 2 ] )
            data = []
            for xy in promptNubar.copyDataToXYs( ) : data += xy
            endfMFList[1][456] += endfFormats.endfDataList( data ) + [ endfFormats.endfSENDLineNumber( ) ]
            totalNubar = promptNubar
            try :
                if( not( totalDelayedNubar is None ) ) : totalNubar = totalNubar + totalDelayedNubar
            except :    # The following is a kludge for some "bad" data.
                if( ( totalNubar.domainMax( unitTo = 'MeV' ) == 30. ) and 
                        ( totalDelayedNubar.domainMax( unitTo = 'MeV' ) == 20. ) ) :
                    totalDelayedNubar[-1] = [ totalNubar.domainMax( ), totalDelayedNubar.getValue( totalDelayedNubar.domainMax( ) ) ]
                totalNubar = totalNubar + totalDelayedNubar
        elif( 'totalNubar' in targetInfo.dict ) :
            totalNubar = targetInfo['totalNubar']
        if( not( totalNubar is None ) ) :
            if( hasattr( totalNubar, 'toENDF6' ) ) :
                totalNubar.toENDF6( 452, endfMFList, flags, targetInfo )
            elif( isinstance( totalNubar, XYs.XYs ) ) :
                n, data = len( totalNubar ), []
                endfMFList[1][452]  = [ endfFormats.endfHeadLine( targetZA, targetInfo['mass'], 0, 2, 0, 0 ) ]
                endfMFList[1][452] += [ endfFormats.endfHeadLine( 0, 0, 0, 0, 1, n ) ]
                endfMFList[1][452] += endfFormats.endfInterpolationList( [ n, 2 ] )
                for xy in totalNubar.copyDataToXYs( ) : data += xy
                endfMFList[1][452] += endfFormats.endfDataList( data ) + [ endfFormats.endfSENDLineNumber( ) ]
            else :
                raise Exception( 'unsupported totalNubar type = %s' % brb.getType( totalNubar ) )

        # add covariances if available:
        if covarianceSuite:
            covarianceSuite.toENDF6( endfMFList, flags, targetInfo )

        for MF in MFs :                 # Build directory information.
            MFData = endfMFList[MF]
            MTs = sorted( MFData.keys( ) )
            for MT in MTs :
                if( ( MF == 1 ) and ( MT == 451 ) ) : continue
                data = MFData[MT]
                directory.append( "%33d%11d%11d%11d" % ( MF, MT, len( data ) - 1, 0 ) )
        endfDoc = self.documentation.get( 'endfDoc' )
        if( endfDoc is None ) :
            docHeader2 = [  ' %2d-%-2s-%3d LLNL       EVAL-OCT03 Unknown' % ( targetZ, fudge.particles.nuclear.elementSymbolFromZ( targetZ ), targetA ),
                            '                      DIST-DEC99                       19990101   ',
                            '----ENDL              MATERIAL %4d' % MAT,
                            '-----INCIDENT %s DATA' % 
                                { 1 : 'NEUTRON', 1001 : 'PROTON', 1002 : 'DEUTERON', 1002 : 'TRITON', 2003 : 'HELION', 2004 : 'ALPHA' }[projectileZA], 
                            '------ENDF-6 FORMAT' ]
            endfDoc = [ 'LLNL ENDL file translated to ENDF6 by FUDGE.', '' ' ************************ C O N T E N T S ***********************' ]
        else :
            docHeader2 = []
            endfDoc = endfDoc.getLines( )
        directory.insert( 0, "%33d%11d%11d%11d" % ( 1, 451, len( directory ) + len( endfDoc ) + 5, 0 ) )  # Add 4 lines of metadata at start, plus current line of directory
        directory.append( 99999 )

        # update the documentation, including metadata on first 4 lines:
        try: self.getReaction('fission'); LFI = True
        except KeyError: LFI = False
        LRP = -1
        if( getattr( self, 'resonances', None ) ) :
            if ( self.resonances.scatteringRadius ): LRP = 0
            elif ( self.resonances.reconstructCrossSection ): LRP = 1
            elif ( self.resonances.unresolved and not self.resonances.resolved
                    and self.resonances.unresolved.tabulatedWidths.forSelfShieldingOnly ): LRP = 1
            else: LRP = 2
        EMAX = max( [ reaction.crossSection.domainMax( unitTo = 'eV' ) for reaction in self ] )
        version = self.styles['evaluated'].getAttribute( 'version' )
        if( version == 'ENDL' ) :           # additional ENDF meta-data. If the library is unknown, use NLIB=-1
            NVER, LREL, NMOD = 1, 1, 1
            NLIB = -1
            temperature = self.attributes['temperature'].getValueAs( 'k * K' )  # k is the Boltzmann constant.
        else :
            NVER, LREL, NMOD = map( int, self.styles['evaluated'].getAttribute( 'version' ).split( '.' ) )    # version stored as '7.2.1'
            NLIB = { "ENDF/B" : 0,      "ENDF/A" : 1,   "JEFF" : 2,     "EFF" : 3,      "ENDF/B (HE)" : 4,  "CENDL" : 5,                    "JENDL" : 6,
                     "SG-23" : 21,      "INDL/V" : 31,  "INDL/A" : 32,  "FENDL" : 33,   "IRDF" : 34,        "BROND (IAEA version)" : 35,    "INGDB-90" : 36,
                     "FENDL/A" : 37,    "BROND" : 41 }.get( self.styles['evaluated'].attributes['library'], -1 )
            temperature = self.attributes['temperature'].getValueAs( 'K' )
        NFOR = 6    # ENDF-6 format
        ITYPE = 0   # other ITYPE sublibraries not yet supported
        NSUB = projectileZA * 10 + ITYPE
        LDRV = 0
        STA = 0
        if( isinstance( self.target, fudge.gnd.xParticle.nuclearLevel ) or self.target.attributes.get( 'unstable' ) ) : STA = 1
        if( targetInfo['LISO'] ) : STA = 1
        docHeader = [ endfFormats.endfHeadLine( targetZA, targetInfo['mass'], LRP, LFI, NLIB, NMOD ),
                endfFormats.endfHeadLine( self.target.getLevelAsFloat( 'eV' ), STA, self.target.getLevelIndex(), targetInfo['LISO'], 0, NFOR ),
                endfFormats.endfHeadLine( self.projectile.getMass( 'eV/c**2' ) / targetInfo['neutronMass'], EMAX, LREL, 0, NSUB, NVER ),
                endfFormats.endfHeadLine( temperature, 0, LDRV, 0, len( endfDoc ), len( directory[:-1] ) ) ]
        new_doc = fudge.gnd.documentation.documentation( 'endf', '\n'.join( docHeader + docHeader2 + endfDoc ) )
        endfMFList[1][451] += endfFormats.toEndfStringList( new_doc )
        
        endfList = [ "%66s%s" % ( " ", endfFormats.endfFENDLine( 1 )[66:75] ) ]
        endfMFList[1][451] += directory
        for MF in MFs :
            MFData = endfMFList[MF]
            MTs = sorted( MFData.keys( ) )
            for MT in MTs :
                data = MFData[MT]
                for i, datum in enumerate( data ) :
                    if( datum == 99999 ) :
                        endfList.append( endfFormats.endfSENDLine( MAT, MF ) )
                    else :
                        endfList.append( '%-66s%4d%2d%3d%5d' % ( datum, MAT, MF, MT, i + 1 ) )
            if( len( MTs ) > 0 ) : endfList.append( endfFormats.endfFENDLine( MAT ) )
        endfList.append( endfFormats.endfMENDLine( ) )
        endfList.append( endfFormats.endfTENDLine( ) )
        endfList.append( '' )
        return( '\n'.join( endfList ) )

    def setAttribute( self, name, value ) :
        """Adds attribute name and its value to the list of attributes."""

        self.attributes[name] = value

    def toString( self, indent = '' ) :
        """Returns a string representation of an reactionSuite."""

        s = indent + self.inputParticlesToReactionString( suffix = ' --> ' )
        indent2 = len( s ) * ' '
        indent3 = ''
        for channel in self :
            s += channel.toString( indent = indent3 )
            indent3 = indent2
        return( s )

def readXML( gndFile ):
    """Read a GND/xml file into fudge.
    Argument 'gndFile' is the path to a GND file.
    Returns a reactionSuite instance containing all data from the file.
    """
    from xml.etree import cElementTree
    rsElement = cElementTree.parse( gndFile ).getroot()
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode
    rsElement = xmlNode( rsElement, xmlNode.etree )
    return parseXMLNode( rsElement )

def parseXMLNode( rsElement ):
    """Translates a <reactionSuite> xml node into a reactionSuite instance.
    Users should use the 'readXML' function instead. """
    format = rsElement.get('format')
    if format not in ('gnd version 1.2',):
        raise Exception("GND file is in an unsupported format ('%s')!" % format)
    styles = [fudge.gnd.miscellaneous.style.parseXMLNode( style ) for style in rsElement.find('styles')]
    documentations = [fudge.gnd.documentation.parseXMLNode(doc) for doc in rsElement.find('documentations')]
    particles = fudge.gnd.xParticleList.parseXMLNode( rsElement.find('particles') )
    projectile = particles.getParticle( rsElement.get('projectile') )
    target = particles.getParticle( rsElement.get('target') )

    rs = reactionSuite( projectile, target, None, None, particles )
    for style in styles: rs.addStyle( style )
    for doc in documentations: rs.addDocumentation( doc )
    rs.version = rsElement.get('version')
    rs.attributes['temperature'] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(rsElement.get('temperature'))

    linkData = {'particles':particles, 'reactionSuite':rs, 'unresolvedLinks':[], 'format':format}

    for child in rsElement:
        if child.tag in ('styles','documentations','particles'):
            continue    # already read above
        elif child.tag == 'aliases':
            for alias in child: rs.addAlias( alias.get('key'), alias.get('value') )
        elif child.tag == 'resonances':
            rs.addResonances( fudge.gnd.resonances.resonances.parseXMLNode( child ) )
        elif child.tag == 'reaction':
            rs.addReaction( fudge.gnd.reaction.parseXMLNode( child, linkData ) )
        elif child.tag == 'partialGammaProduction':
            rs.addPartialGammaProduction( fudge.gnd.partialGammaProduction.parseXMLNode( child, linkData ) )
        elif child.tag == 'summedReaction':
            rs.addSummedReaction( fudge.gnd.summedReaction.parseXMLNode( child, linkData ) )
        elif child.tag == 'fissionComponent':
            fc = fudge.gnd.reaction.parseXMLNode( child, linkData )
            fc.__class__ = fudge.gnd.fissionComponent.fissionComponent
            fc.moniker = fudge.gnd.reactions.base.fissionComponentToken
            rs.addFissionComponent( fc )
        elif child.tag == 'production':
            rs.addProductionReaction( fudge.gnd.production.parseXMLNode( child, linkData ) )
        elif child.tag == 'referredData':
            rs.referredData = fudge.gnd.referredData.parseXMLNode( child, linkData )
            rs.referredData.parent = rs
        else:
            print("Warning: encountered unexpected element '%s' in reactionSuite!" % child.tag)

    # fix links:
    from fudge.gnd.link import follow
    for path, quant in linkData['unresolvedLinks']:
        if type(quant) is fudge.gnd.link.link: quant.link = follow(path, rs)
        else: quant.setReference( follow(path, rs) )

    return rs
