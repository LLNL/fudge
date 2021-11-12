# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
reactionSuite.py contains the 'reactionSuite' class that in turn holds all reactions for a given target/projectile.
reactionSuite is the top-level class for the GNDS structure.
"""

import re
import os

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs import database as PoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.groups import chemicalElement as chemicalElementPoPsModule
from PoPs.decays import misc as decaysMiscPoPsModule
from PoPs.families import unorthodox as unorthodoxPoPsModule

from xData import formatVersion as formatVersionModule
from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import standards as standardsModule
from xData import XYs as XYs1dModule
from xData.Documentation import documentation as documentationModule

from fudge.core.utilities import guessInteraction as guessInteractionModule
from fudge.processing import group as groupModule
from fudge.processing import nuclearPlusCoulombInterference as nuclearPlusCoulombInterferenceModule
from fudge.processing.deterministic import multiGroupUpScatter as multiGroupUpScatterModule

from fudge import resonances as resonancesModule
from fudge import suites as suitesModule
from fudge import sums as sumsModule
from fudge import styles as stylesModule
from fudge import warning as warningModule
from fudge import product as productModule
from . import institution as institutionModule

from .reactionData import crossSection as crossSectionModule
from .channelData import Q as QModule
from .productData import energyDeposition as energyDepositionModule
from .productData.distributions import angular as angularModule

from brownies.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule

import fudge

__metaclass__ = type

class Interaction :
    """Module (class) that defines the all interaction tokens."""

    nuclear = 'nuclear'
    atomic = 'atomic'
    TNSL = 'TNSL'
    LLNL_TNSL = 'LLNL_TNSL'

    allowed = ( nuclear, atomic, TNSL, LLNL_TNSL )

class nullDevice :
    """For internal use. Used in methods to set logFile when one is not entered."""

    def write( self, **kwargs ) : pass

class reactionSuite( ancestryModule.ancestry ) :
    """
    This is the main class for a gnds projectile/target object. It contains
        * optional external files (e.g. for covariances or binary data stored in another file)
        * a 'styles' section describing all evaluated and processed data styles, documentation, etc.
        * a 'PoPs' database with a list of all particles encountered in the file
        * optional resonance data
        * a list of reactions
        * additional lists: summed reactions (e.g. total), production reactions, 'orphan products', etc.
    """

    moniker = 'reactionSuite'
    ancestryMembers = ( '[externalFiles', '[styles', 'PoPs', 'resonances', '[reactions', '[orphanProducts', 'sums', '[productions', '[incompleteReactions', '[fissionComponents', 'applicationData' )
    childNodeOrder = {
            formatVersionModule.version_1_10 : ( suitesModule.externalFiles.moniker,                stylesModule.styles.moniker, 
                                                 'documentations',                                  PoPsModule.database.moniker,
                                                 resonancesModule.resonances.resonances.moniker,    suitesModule.reactions.moniker, 
                                                 suitesModule.orphanProducts.moniker,               sumsModule.sums.moniker, 
                                                 suitesModule.fissionComponents.moniker,            suitesModule.productions.moniker,
                                                 suitesModule.incompleteReactions.moniker,          suitesModule.applicationData.moniker ),
            formatVersionModule.version_2_0_LLNL_4 : (  suitesModule.externalFiles.moniker,                 stylesModule.styles.moniker,
                                                        PoPsModule.database.moniker,
                                                        resonancesModule.resonances.resonances.moniker,     suitesModule.reactions.moniker, 
                                                        suitesModule.orphanProducts.moniker,                sumsModule.sums.moniker, 
                                                        suitesModule.fissionComponents.moniker,             suitesModule.productions.moniker,
                                                        suitesModule.incompleteReactions.moniker,           suitesModule.applicationData.moniker ) }

    def __init__( self, projectile, target, evaluation, interaction = None, formatVersion = formatVersionModule.default, style = None,
            projectileFrame = standardsModule.frames.labToken, MAT = None, PoPs = None, sourcePath = None ) :
        """
        Creates a new reactionSuite object, containing all reaction data for a projectile/target combination.

        :param projectile:      Projectile ID (string)
        :param target:          Target ID (string)
        :param evaluation:      The evaluation for this protare (string)
        :param interaction:     The GNDS interaction for the protare (string)
        :param formatVersion:   GNDS format version (string)
        :param style:           Indicates style of data in this reactionSuite (evaluated, processed, etc.).  Instance of styles.style
        :param projectileFrame: Frame of the projectile. Allowed values are 'lab' or 'centerOfMass'
        :param MAT:             Integer ENDF MAT number (only needed for writing back to ENDF-6)
        :param PoPs:            A PoPs database that is copied to self's PoPs database
        :return:
        """

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( projectile, str ) ) ) : raise TypeError( 'projectile ID not a string' )
        self.__projectile = projectile

        if( not( isinstance( target, str ) ) ) : raise TypeError( 'target ID not a string' )
        self.__target = target

        if( not( isinstance( evaluation, str ) ) ) : raise TypeError( 'evaluation not a string' )
        self.__evaluation = evaluation

        self.interaction = interaction
        self.__guessedInteraction = False

        if( formatVersion not in formatVersionModule.allowed ) : raise Exception( "Unsupported GNDS PoPs structure '%s'!" % str( formatVersion ) )
        self.formatVersion = formatVersion

        self.__projectileFrame = projectileFrame

        self.MAT = MAT

        self.__sourcePath = sourcePath

        self.__externalFiles = suitesModule.externalFiles()
        self.__externalFiles.setAncestor( self )

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )
        if( style is not None ) : self.styles.add( style )

        if( PoPs is not None ) :
            self.__PoPs = PoPs.copy( )
        else :
            self.__PoPs = PoPsModule.database( 'protare_internal', '1.0' )
        self.__PoPs.setAncestor( self )

        self.__resonances = None

        self.__reactions = suitesModule.reactions()
        self.__reactions.setAncestor( self )

        self.__orphanProducts = suitesModule.orphanProducts( )
        self.__orphanProducts.setAncestor( self )

        self.__sums = sumsModule.sums()
        self.__sums.setAncestor( self )

        self.__productions = suitesModule.productions()
        self.__productions.setAncestor( self )

        self.__incompleteReactions = suitesModule.incompleteReactions()
        self.__incompleteReactions.setAncestor( self )

        self.__fissionComponents = suitesModule.fissionComponents()
        self.__fissionComponents.setAncestor( self )

        self.__applicationData = suitesModule.applicationData()
        self.__applicationData.setAncestor( self )

        self._externalLinks = []    # keep track of links to external files

    def __iter__( self ) :

        for reaction in self.reactions : yield reaction
        for orphanProduct in self.orphanProducts : yield orphanProduct
        for _sum in self.sums.crossSectionSums : yield _sum
        for reaction in self.fissionComponents : yield reaction
        for reaction in self.productions : yield reaction
        for reaction in self.incompleteReactions : yield reaction

    def __str__( self ) :
        """See method toString."""

        string = self.inputParticlesToReactionString( ) + ( ' (%s)' %  self.evaluation )
        
        return( string )

    @property
    def compound_nucleus( self ) :
        """
        Compute the compound nucleus id corresponding to this reaction.

        Note, we don't check to see if the compound nucleus is a valid concept, so
        this will fail on electro-atomic, photo-atomic and natural element target data
        """

        ZAproj = chemicalElementMiscPoPsModule.ZA( self.PoPs[ self.projectile ] )
        ZAtarg = chemicalElementMiscPoPsModule.ZA( self.PoPs[ self.target ] )
        Zcompound, Acompound = divmod( ZAproj + ZAtarg, 1000 )
        return chemicalElementMiscPoPsModule.idFromZAndA(Zcompound, Acompound)

    @property
    def format( self ) :

        return( self.formatVersion )

    @property
    def projectile( self ) :
        """Returns self's projectile."""

        return( self.__projectile )

    @property
    def target( self ) :
        """Returns self's target."""

        return( self.__target )

    @property
    def projectileFrame( self ) :

        return( self.__projectileFrame )

    @property
    def evaluation( self ) :
        """Returns self's evaluation."""

        return( self.__evaluation )

    @property
    def interaction( self ) :
        """Returns self's interaction."""

        return( self.__interaction )

    @interaction.setter
    def interaction( self, value ) :

        if( value is None ) :
            value = Interaction.nuclear
            print( 'Need to specify interaction when calling reactionSuite.__init__. Setting it to "%s". Please update your code as in the future this will execute a raise.' % value )

        if( value not in Interaction.allowed ) : raise TypeError( 'Invalid interaction = "%s"' % value )
        self.__interaction = value

    @property
    def guessedInteraction(self) :
        '''Returns the guessedInteraction value. This is for pre GNDS 2.0 use only and should not be used with GNDS 2.0 or higher
        (i.e., it should always be False for GNDS 2.0 or higher.'''

        return self.__guessedInteraction

    @guessedInteraction.setter
    def guessedInteraction(self, value) :
        '''Set self's guessedInteraction member. This is for pre GNDS 2.0 use only and should not be used with GNDS 2.0 or higher.'''

        self.__guessedInteraction = value

    @property
    def sourcePath( self ) :
        """Returns the sourcePath member which is the path to the reactionSuite file for self if self is from a file."""

        return( self.__sourcePath )

    @property
    def externalFiles( self ) :
        """Returns the externalFiles instance."""

        return self.__externalFiles

    @property
    def styles( self ) :
        """Returns the PoPs instance."""

        return self.__styles
    
    @property
    def PoPs( self ) :
        """Returns the PoPs instance."""

        return( self.__PoPs )

    @property
    def resonances(self):
        return self.__resonances

    @resonances.setter
    def resonances( self, _resonances ) :
        """Set resonances to the specified resonances instance."""

        if not isinstance(_resonances, resonancesModule.resonances.resonances):
            raise TypeError("Must be resonances instance, not '%s'" % type(_resonances))
        self.__resonances = _resonances
        self.__resonances.setAncestor( self )

    @property
    def reactions( self ) :
        """Returns the reactions instance."""

        return self.__reactions

    @property
    def orphanProducts( self ) :
        """Returns the orphanProducts instance."""

        return self.__orphanProducts

    @property
    def sums( self ):
        """Returns the sums instance."""

        return self.__sums

    @property
    def productions( self ) :
        """Returns the productions instance."""

        return self.__productions

    @property
    def incompleteReactions( self ) :
        """Returns the incompleteReactions instance."""

        return self.__incompleteReactions

    @property
    def fissionComponents( self ) :
        """Returns the fissionComponents instance."""

        return self.__fissionComponents

    @property
    def applicationData( self ) :
        """Returns the applicationData instance."""

        return self.__applicationData

    @property
    def domainMin( self ) :

        return( self.reactions[0].crossSection.domainMin )

    @property
    def domainMax( self ) :

        return( self.reactions[0].crossSection.domainMax )

    @property
    def domainUnit( self ) :

        return( self.reactions[0].crossSection.domainUnit )

    def isThermalNeutronScatteringLaw( self ) :

        for reaction in self.reactions :
            if( reaction.isThermalNeutronScatteringLaw( ) ) : return( True )
        return( False )

    def photonBranchingData( self ) :

        branchingData = {}
        for chemicalElement in self.PoPs.chemicalElements :
            for isotope in chemicalElement :
                branchingData[isotope.symbol] = decaysMiscPoPsModule.photonBranchingData( self.PoPs, isotope.symbol )
        return( branchingData )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        self.styles.convertUnits( unitMap )
        if self.resonances is not None: self.resonances.convertUnits( unitMap )
        self.reactions.convertUnits( unitMap )
        self.orphanProducts.convertUnits( unitMap )
        self.sums.convertUnits( unitMap )
        self.productions.convertUnits( unitMap )
        self.PoPs.convertUnits( unitMap )
        self.fissionComponents.convertUnits( unitMap )

    def check( self, **kwargs ) :
        """
        Check all data in the reactionSuite, returning a gnds.warning.context object with list of warnings.

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
            'dEnergyBalanceAbsolute'     '1 eV'
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

        from fudge import warning

        options = {
                'branchingRatioSumTolerance': 1e-6,
                'dQ': '1e-3 MeV',
                'dThreshold': '1e-3 MeV',
                'crossSectionEnergyMin': '1e-5 eV',
                'crossSectionEnergyMax': '20 MeV',
                'crossSectionOnly': False,
                'crossSectionMaxDiff': 1e-3,
                'multiplicityMaxDiff': 1e-3,
                'transportables': ('n',),
                'normTolerance': 1e-5,
                'checkEnergyBalance': True,
                'reconstructResonances': True,
                'dEnergyBalanceRelative': 1e-3,
                'dEnergyBalanceAbsolute': '1 eV',
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
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError("check() received unknown keyword argument '%s'" % key)

        warnings = []

        if isinstance( self.target, PoPsModule.unorthodoxModule.particle ):
            warnings.append( warning.UnorthodoxParticleNotImplemented( self.target))
            return warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        evaluatedStyle = self.styles.getEvaluatedStyle()
        if evaluatedStyle is None:
            warnings.append( warning.NotImplemented( "Checking currently only supported for 'evaluated' style" ) )
            return warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)

        # assemble some useful info, to be handed down to children's check() functions:
        incidentEnergyUnit = 'eV'
        massUnit = 'eV/c**2'
        projectile = self.PoPs[self.projectile]
        projectileZ, dummy, projectileZA, dummy = chemicalElementMiscPoPsModule.ZAInfo( projectile )
        target = self.PoPs[self.target]
        if hasattr(target,'id') and target.id in self.PoPs.aliases:
            target = self.PoPs[ target.pid ]
        targetZ, dummy, targetZA, dummy = chemicalElementMiscPoPsModule.ZAInfo( target )
        compoundZA = targetZA + projectileZA
        CoulombChannel = ( targetZ != 0 ) and ( projectileZ != 0 )
        elementalTarget=False
        if hasattr( target, 'getMass' ):
            projectileMass_amu = projectile.getMass( 'amu' )
            targetMass_amu = target.getMass( 'amu' )
            kinematicFactor = ( targetMass_amu + projectileMass_amu ) / targetMass_amu    # (M+m)/M

            projectileMass = projectile.getMass( massUnit )
            targetMass = target.getMass( massUnit )
            availableEnergy_eV = projectileMass + targetMass
        else:
            # For elemental targets, calculating these factors doesn't make sense since there is no defined target mass
            kinematicFactor=1.0
            elementalTarget=True
            availableEnergy_eV=None
        if isinstance(target,chemicalElementPoPsModule.chemicalElement): elementalTarget=True

        info = { 'reactionSuite': self, 'kinematicFactor': kinematicFactor, 'compoundZA': compoundZA,
                'availableEnergy_eV': availableEnergy_eV, 'CoulombChannel': CoulombChannel, 'style': evaluatedStyle,
                'reconstructedStyleName': None }
        info.update( options )
        if elementalTarget:
            # For elemental targets, calculating energy balance isn't possible
            info['checkEnergyBalance']=False

        externalFileWarnings = []
        for externalFile in self.externalFiles:
            from LUPY import checksums
            fullPath = os.path.join(
                os.path.split(self.sourcePath)[0], externalFile.path
            )
            if not os.path.exists(fullPath):
                warnings.append(warning.missingExternalFile(externalFile.path))
                continue

            if externalFile.checksum is not None:
                algorithm = externalFile.algorithm
                checksum = checksums.checkers[algorithm].from_file(fullPath)
                if checksum != externalFile.checksum:
                    warnings.append(
                        warning.wrongExternalFileChecksum(externalFile.path,
                                                          externalFile.checksum,
                                                          checksum))

        if externalFileWarnings:
            warnings.append( warning.context('externalFiles', externalFileWarnings))

        if self.resonances is not None:
            resonanceWarnings = self.resonances.check( info )
            if resonanceWarnings:
                warnings.append( warning.context('resonances', resonanceWarnings) )
            if ( options['reconstructResonances'] and self.resonances.reconstructCrossSection ):
                info['reconstructedStyleName'] = self.styles.getTempStyleNameOfClass( stylesModule.crossSectionReconstructed )
                # convert resonance parameters to pointwise data, interpolable to .1 percent:
                reconstructedStyle = stylesModule.crossSectionReconstructed( info['reconstructedStyleName'], derivedFrom=evaluatedStyle.label )
                try: self.reconstructResonances( reconstructedStyle, accuracy=0.001, thin=False, verbose=False )
                except Exception as e:
                    warnings.append( warning.ExceptionRaised( "when reconstructing resonances: %s" % e ) )
                    if info['failOnException']: raise

        if info['checkEnergyBalance']:
            # setup options for calculating average product energy and momentum
            averageProductDataStyle = stylesModule.averageProductData(
                    label = self.styles.getTempStyleNameOfClass( stylesModule.averageProductData ),
                    derivedFrom=evaluatedStyle.label )

            self.styles.add( averageProductDataStyle )
            info['averageProductDataStyle'] = averageProductDataStyle
            info['averageProductDataArgs'] = { 'verbosity':1,   # additional required arguments
                        'incrementalIndent':'  ', 'energyAccuracy':1e-6, 'momentumAccuracy':1e-6,
                        'incidentEnergyUnit': incidentEnergyUnit, 'massUnit': massUnit,
                        'momentumUnit': incidentEnergyUnit + '/c',
                        'projectileMass': projectileMass, 'targetMass': targetMass,
                        'reactionSuite':self, 'photonBranchingData': self.photonBranchingData(),
                        'isInfiniteTargetMass' : isinstance( target, chemicalElementPoPsModule.chemicalElement ) }

        if self.projectile == 'n' and self.target != 'n':
            # test Wick's limit: 0-degree elastic xsc >= ( total xsc * k/4pi )^2
            elastic = self.getReaction('elastic')
            total = self.getReaction('total')
            if total is None:
                warnings.append( warning.testSkipped("Wick's limit", "can't find reaction 'total'" ) )
            else:
                try :
                    if( info['reconstructedStyleName'] in elastic.crossSection ) :
                        elastic_xsc = elastic.crossSection[info['reconstructedStyleName']]
                        total_xsc = total.crossSection[info['reconstructedStyleName']]
                    else:
                        elastic_xsc = elastic.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8,upperEps=1e-8)
                        total_xsc = total.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8,upperEps=1e-8)
                    elastic_distribution = elastic.outputChannel.getProductWithName('n').distribution[evaluatedStyle.label].angularSubform
                    if isinstance( elastic_distribution, fudge.productData.distributions.angular.XYs2d ):
                        linearized = elastic_distribution.toPointwise_withLinearXYs()
                        forward_scattering = crossSectionModule.XYs1d( data=[], axes=crossSectionModule.XYs1d.defaultAxes() )
                        for energy_in in linearized.getEnergyArray():
                            forward_scattering.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                    elif isinstance( elastic_distribution, fudge.productData.distributions.angular.regions2d ):
                        forward_scattering = crossSectionModule.regions1d( axes=crossSectionModule.XYs1d.defaultAxes() )
                        for region in elastic_distribution:
                            ptw = crossSectionModule.XYs1d( data=[] )
                            linearized = region.toPointwise_withLinearXYs()
                            for energy_in in linearized.getEnergyArray():
                                ptw.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                            forward_scattering.append( ptw )
                        forward_scattering = forward_scattering.toPointwise_withLinearXYs( accuracy = 1e-5,
                                lowerEps = 1e-8, upperEps = 1e-8 )

                    mutualDomain = list( zip(
                        *[dat.domain() for dat in (elastic_xsc, total_xsc, forward_scattering)]) )
                    mutualDomain = (max(mutualDomain[0]), min(mutualDomain[1]))
                    egrid = [e for e in forward_scattering.domainGrid if mutualDomain[0] <= e <= mutualDomain[1]]

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

        particleWarnings = self.PoPs.check( **info )
        if particleWarnings: warnings.append( warning.context('PoPs', particleWarnings) )

        info['isInfiniteTargetMass'] = isinstance( target, chemicalElementPoPsModule.chemicalElement )
        for reaction in self :
            reactionWarnings = reaction.check( info )
            if reactionWarnings: warnings.append( warning.context('%s label %s'
                % (reaction.moniker, reaction.label), reactionWarnings ) )

        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info

        if options['cleanUpTempFiles']:
            if info['reconstructedStyleName'] is not None: self.removeStyle( info['reconstructedStyleName'] )

        return result

    def diff( self, other ) :

        diffResults1 = warningModule.diffResults( self, other )

        if( not( isinstance( other, reactionSuite ) ) ) : raise ValueError( 'Other not a reactionSuite instance.' )
        if( self.domainUnit != other.domainUnit ) : raise ValueError( 'Energy units not the same.' )

        for style in self.styles :
            if( not( isinstance( style, ( stylesModule.evaluated, stylesModule.crossSectionReconstructed ) ) ) ) :
                raise TypeError( 'Style "%s" not allowed when diffing.' % style.moniker )

        for style in other.styles :
            if( not( isinstance( style, ( stylesModule.evaluated, stylesModule.crossSectionReconstructed ) ) ) ) :
                raise TypeError( 'Style "%s" not allowed when diffing.' % style.moniker )

        self.reactions.diff( other.reactions, diffResults1 )

        return( diffResults1 )

    def elasticReactionLabel( self ) :
        """
        Returns the label for the elastic channel. When constructing the label, if target is an alias, use the alias' pid; otherwise, use the target id.
        """

        return( '%s + %s' % ( self.projectile, self.targetID( protareID = False ) ) )

    def targetID( self, protareID = False ) :
        """
        Returns the ID for the target. If the target is a nuclear meta-stable, then the target ID as listed by the 'target' attribute for the
        reactionSuite (i.e., protare) will be the meta-stable's ID while everywhere else it will be the nuclide excitation ID. Otherwise,
        the target ID as listed by the 'target' attribute for the reactionSuite is the same as it is everywhere else. For example, for the
        protare 'n + Am242', the target will be listed as 'Am242' everywhere, while for the protare 'n + Am242_m1', the target will be listed
        as 'Am242_m1' by the 'target' attribute for the reactionSuite and 'Am242_e2' everywhere else. Note, for the protare 'n + Am242_m1'
        the PoPs.aliases instance must have the alias 'Am242_m1' for 'Am242_e2'. If protareID is True, the target ID as listed by the 'target'
        attribute for the reactionSuite will be returned, otherwise, the target ID list everywhere else will be returned.
        """

        if( protareID ) : return( self.target )
        return( self.PoPs.final( self.target ).id )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. reactionSuite contains several different lists,
        so may need to descend into those to find desired entity.
        """

        if entityName in ('reaction','summedReaction','fissionComponent','production','reactionSum'):
            for entity in getattr( self, entityName+'s' ):
                if getattr( entity, attribute, None ) == value:
                    return entity
        return( ancestryModule.ancestry.findEntity( self, entityName, attribute, value ) )

    def findLinks( self ) :

        links = []
        for ancestryMember in self.ancestryMembers :
            if( hasattr( self, ancestryMember ) ) :
                item = getattr( self, ancestryMember )
                if( hasattr( item, 'findLinks' ) ) : getattr( item, 'findLinks' )( links )
        return( links )

    def fixLinks( self ) :

        ENDFconversions = []
        if( 'LLNL' in self.applicationData ) :
            ENDFconversions = [ LLNL_data for LLNL_data in self.applicationData['LLNL'] if isinstance( LLNL_data, ENDFconversionFlagsModule.ENDFconversionFlags ) ]
            if( len( ENDFconversions ) > 1 ) : raise Exception( 'Too many ENDFconversionFlags in applicationData["LLNL"] section.' )
            if( len( ENDFconversions ) == 1 ) : ENDFconversions = ENDFconversions[0]

        links = self.findLinks( )
        for instance, link, xlink in links :
            link2 = link
            if( xlink[0] == '/' ) :                     # Cannot look for relative xlink as it is relative to link which is not known.
                try :
                    link2 = self.followXPath( instance.path )
                except ancestryModule.XPathNotFound :
                    link2 = None
                except :
                    raise

            if( link2 is None ) :
                if( isinstance( instance, ENDFconversionFlagsModule.conversion ) ) :
                    ENDFconversions.flags.pop( ENDFconversions.flags.index( instance ) )
                elif( isinstance( instance, sumsModule.add ) ) :
                    summands = instance.ancestor
                    summands.pop( summands.index( instance ) )
            elif( link is not link2 ) :
                instance.link = link2

    def hasParticle( self, name ) :

        return( name in self.PoPs )

    def heatCrossSections( self, style, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001, heatAllPoints = False,
            doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False, verbose = 0 ) :

        EMin = self.domainMin
        for reaction in self.reactions :
            reaction.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero, verbose = verbose )

        for orphanProduct in self.orphanProducts :
            orphanProduct.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero, verbose = verbose )

        for production in self.productions :
            production.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero, verbose = verbose )

        for _sum in self.sums.crossSectionSums :
            _sum.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero, verbose = verbose )

        for fissionComponent in self.fissionComponents :
            fissionComponent.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero, verbose = verbose )

    def getParticle( self, name ) :

        return( self.PoPs[ name ] )

    def getMassRatio( self ) :

        if not hasattr( self, '__massRatio' ):
            M = self.PoPs[self.target].mass[0].float( 'amu' )
            m = self.PoPs[self.projectile].mass[0].float( 'amu' )
            self.__massRatio = (M / (M+m))
        return self.__massRatio

# BRB6 This does not belong here.
    def getIsotopeName( self, *args ):
        """
        Return name of compound nucleus formed by nuc1+nuc2+...
        if a nucleus in 'args' starts with '-', subtract instead.
        """

# CALEB: This logic need to use PoPs stuff.
        def parse(name):
            if( IDsPoPsModule.photon in name ) : return 1,0,0,1
            sign, mult, symbol, A = re.search(r"([-]?)([0-9]+)?([a-zA-Z]+)(_natural|[0-9]*)", name).groups()
            if not mult: mult = 1
            if( ( symbol == IDsPoPsModule.neutron ) and not A ) : A = 1
            if A!='_natural': A=int(A)
            if( name == IDsPoPsModule.neutron ) :
                Z = 0
            else :
                Z = chemicalElementMiscPoPsModule.ZFromSymbol[symbol]
            return( int( sign + '1' ), Z, A, int( mult ) )

        retZ, retA = 0,0
        try:
            for nucleus in args:
                sign, Z, A, mult = parse(nucleus)
                retZ += sign*mult*Z
                if '_natural' in (A,retA):
                    retA='_natural'
                else:
                    retA += sign*mult*A
            return( '%s%s' % ( chemicalElementMiscPoPsModule.symbolFromZ[retZ], retA ) )
        except:
            print ("      WARNING: couldn't extract isotope name from product list!")

# BRB6 many hardwired in this method.
    def getReaction( self, channel ):
        """
        Search list of reactions for a specified channel.
        The 'channel' argument should be either a reaction type ('elastic','capture','fission', etc) or
        the list of outgoing particles ('n + Pu239' for example).

        If 'channel' is an int, then we assume it's an ENDF MT

        Raises 'KeyError' if specified channel can't be found.
        """
        import numpy

        # If channel is an int, it might be an ENDF MT, so check those first
        if isinstance(channel, (int, numpy.integer)): # must be an ENDF MT
            for reactionList in self.reactions, self.sums.crossSectionSums, self.fissionComponents, self.productions:
                for reaction in reactionList:
                    if reaction.ENDF_MT == channel: return reaction
            return None

        # translate special channel names:
        if channel == 'elastic': channel = channel_tr = self.elasticReactionLabel( )
        elif channel=='capture': channel_tr='z,gamma'
        else: channel_tr = channel

        # check if 'channel' == one of the fudge reactions:
        chStrings, reacs = [], []
        for reaction in self :
            if( str( reaction ) == channel ) : return( reaction )
            chStrings.append( str( reaction ) )
            reacs.append( reaction )

        # make list containing a set() of products for each reaction.
        chStrings = [tmp.replace('(','').replace(')','').replace('->','+') for tmp in chStrings]
        chSets = [set(tmp.split(' + ')) for tmp in chStrings]
        chSetsNoPhoton = [tmp.copy() for tmp in chSets]
        for tmp in chSetsNoPhoton: tmp.discard('photon')
        if set(channel.split(' + ')) in chSets:
            return reacs[ chSets.index( set(channel.split(' + ')) ) ]
        if set(channel.split(' + ')) in chSetsNoPhoton:
            return reacs[ chSetsNoPhoton.index( set(channel.split(' + ')) ) ]

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
                print("Can't determine fission genre from '%s'" % channel_fiss)
                genre = None
            retVal  = [ r1 for r1 in self.reactions         if r1.fissionGenre == genre ]
            retVal += [ r1 for r1 in self.fissionComponents if r1.fissionGenre == genre ]
            if len(retVal)==1:
                return retVal[0]
        else:
            # check if channel is in form '(z,2na)', 'n,4n' or similar:
            patt = re.match('^[(]?[znpdthag],([a-zA-Z0-9]+)[)]?$', channel_tr)
            if patt:
                thisChannelSet = set()
                match = re.findall('([1-9]?)(gamma|He3|[npagdt]?)[+]?', patt.groups()[0] )
                for mul, prod in match:
                    if not prod : continue
                    prod = { 'g' : IDsPoPsModule.photon , 'gamma' : IDsPoPsModule.photon, 'n' : 'n', 'p' : 'H1',
                             'd' : 'H2', 't' : 'H3', 'h' : 'He3', 'a' : 'He4' }[prod]
                    if mul: prod = mul + prod
                    if prod in thisChannelSet:
                        raise KeyError("Please specify multiplicity explicitly ('z,2n' instead of 'z,nn')")
                    thisChannelSet.add(prod)
                if not thisChannelSet:
                    return None
                # also add residual to the set:
                proj, target = self.projectile, self.target
                thisChannelSet.add( self.getIsotopeName( *([proj, target] + ['-'+a for a in thisChannelSet]) ) )
                if thisChannelSet in chSets:
                    return reacs[ chSets.index( thisChannelSet ) ]
                if thisChannelSet in chSetsNoPhoton:
                    return reacs[ chSetsNoPhoton.index( thisChannelSet ) ]

        return None

    def calculateAverageProductData( self, style, indent = '', additionalReactions = [], **kwargs ) :
        """
        Calculate average energy and momentum data for all products of all reactions.
        Resulting data are stored within each product. Example usage is:

        from fudge import reactionSuite as reactionSuiteModule
        from fudge import styles as stylesSuiteModule
        reactionSuite = reactionSuiteModule.readXML( "16.xml" )

        style = stylesModule.averageProductData( label = 'productData', 'eval' )
        reactionSuite.calculateAverageProductData( style )

        :param style: The style to use.
        :param indent: string; The amount to indent and verbose output.
        :param kwargs: string; All other parameters.
        :return:
        """

        if( not( isinstance( style, stylesModule.averageProductData ) ) ) : raise TypeError( 'Invalid style' )

        verbosity = kwargs.get( 'verbosity', 0 )
        kwargs['verbosity'] = verbosity

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        kwargs['incrementalIndent'] = incrementalIndent
        indent2 = indent + incrementalIndent

        logFile = kwargs.get( 'logFile', nullDevice( ) )
        kwargs['logFile'] = logFile

        energyAccuracy = kwargs.get( 'energyAccuracy', 1e-5 )
        kwargs['energyAccuracy'] = energyAccuracy
        momentumAccuracy = kwargs.get( 'momentumAccuracy', 1e-3 )
        kwargs['momentumAccuracy'] = momentumAccuracy

        kwargs['incidentEnergyUnit'] = self.domainUnit
        kwargs['momentumUnit'] = kwargs['incidentEnergyUnit'] + '/c'
        kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'

        kwargs['projectileMass'] = self.PoPs[self.projectile].getMass( kwargs['massUnit'] )
        target = self.PoPs[self.target]
        kwargs['isInfiniteTargetMass'] = isinstance( target, chemicalElementPoPsModule.chemicalElement )
        if( kwargs['isInfiniteTargetMass'] ) :
            kwargs['targetMass'] = None
        else :
            if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
            kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        kwargs['reactionSuite'] = self

        for reaction in self.reactions : reaction.calculateAverageProductData( style, indent = indent2, **kwargs )
        for reaction in self.orphanProducts : reaction.calculateAverageProductData( style, indent = indent2, **kwargs )
        for reaction in additionalReactions : reaction.calculateAverageProductData( style, indent = indent2, **kwargs )

    def buildCrossSectionSum( self, sumLabel, ENDF_MT, styleLabel, reactionKeys, Q ) :
        """
        For the list of indices and/or labels in reactionKeys, calcualate the summed crosstions and create a crossSectionSum
        for it. If insert is True, the crossSectionSum is added to self's sums/crossSectionSums node. If N keys match the
        same reaction in this, that reaction will be added N times.

        :param sumLabel: Label for the constructed crossSectionSum node
            :type sumLabel: str
        :param ENDF_MT: The ENDF MT for the crossSectionSum node
            :type ENDF_MT: int
        :param styleLabel: The label for the function1d added to the Q and crossSection nodes
            :type styleLabel: str
        :param reactionKeys: All reactions with an index or label matching an item of this list are added to cross section sum.
            :type reactionKeys: iterator, list or tuple
        :param Q: The Q value for the cross section sum
            :type Q: float

        :return: A crossSectionSum instance representing the requested reaction cross section sums
        :rtype: crossSectionSum instance
        """

        crossSectionSum = sumsModule.crossSectionSum( sumLabel, ENDF_MT )

        sum = XYs1dModule.XYs1d( axes = crossSectionModule.defaultAxes( self.domainUnit ) )
        for key in reactionKeys :
            reaction = self.reactions[key]
            sum += reaction.crossSection.toPointwise_withLinearXYs( lowerEps = 1e-6, upperEps = 1e-6 )
            crossSectionSum.summands.append( sumsModule.add( link = reaction.crossSection ) )

        crossSectionSum.crossSection.add( crossSectionModule.XYs1d( data = sum, label = styleLabel, axes = crossSectionModule.defaultAxes( self.domainUnit ) ) )

        crossSectionSum.Q.add( QModule.constant1d( Q, domainMin = sum.domainMin, domainMax = sum.domainMax, axes = QModule.defaultAxes( self.domainUnit ), label = styleLabel ) )

        return( crossSectionSum )

    def inputParticlesToReactionString( self, prefix = "", suffix = "" ) :

        return( "%s%s + %s%s" % ( prefix, self.projectile, self.target, suffix ) )

    def cullStyles( self, styleName, removeDocumentation, removeApplicationData ) :
        """
        Removes all forms from each component except for the requested style or, if it is not present
        its nearest derived from style.
        """

        style = self.__styles[styleName]
        if( isinstance( style, stylesModule.heatedMultiGroup ) ) :
            styleList = [ style ]
        if( isinstance( style, stylesModule.griddedCrossSection ) ) :
            styleList = []
            while( style is not None ) :
                styleList.append( style )
                style = style.derivedFromStyle
        else :
            raise TypeError( 'style name "%s" is for type "%s" which is not supported for culling.' % ( styleName, style.moniker ) )

        self.__reactions.cullStyles( styleList )
        self.__orphanProducts.cullStyles( styleList )
        self.__sums.cullStyles( styleList )
        self.__productions.cullStyles( styleList )
        self.__incompleteReactions.cullStyles( styleList )
        self.__fissionComponents.cullStyles( styleList )
# FIXME What about self.resonances

    def CoulombPlusNuclearMuCutoffs( self ) :

        muCutoffs = []
        for style in self.styles :
            if( isinstance( style, stylesModule.CoulombPlusNuclearElasticMuCutoff ) ) : muCutoffs.append( style.muCutoff )

        return( list( set( muCutoffs ) ) )

    def hasFission( self ) :

        for reaction in self.__reactions:
            if reaction.isFission(): return True

        return False

    def loadCovariances(self):
        """
        Load all external files of type 'covarianceSuite', and resolve links between self and covarianceSuites
        @return: list of loaded covarianceSuites
        """
        from LUPY import GNDSType as GNDSTypeModule
        from fudge.covariances import covarianceSuite as covarianceSuiteModule
        covariances = []
        for external in self.externalFiles:
            # TODO submit proposal to include file type attribute in the externalFile GNDS node so we don't need to run GNDSTypeModule.type
            gndsType, metadata = GNDSTypeModule.type(external.realpath())
            if gndsType == covarianceSuiteModule.covarianceSuite.moniker:
                covariances.append(covarianceSuiteModule.readXML(external.realpath(), reactionSuite=self))
        return covariances

    def partialProductionIntegral( self, productID, energyIn, energyOut = None, muOut = None, phiOut = None, frame = standardsModule.frames.labToken,
                LegendreOrder = 0, **kwargs ) :

        partialProductionIntegralSum = 0.0

        for reaction in self.reactions :
            partialProductionIntegralSum += reaction.partialProductionIntegral( self, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )
        for orphanProduct in self.orphanProducts :
            partialProductionIntegralSum += orphanProduct.partialProductionIntegral( self, productID, energyIn, energyOut = energyOut, muOut = muOut, 
                    phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder, **kwargs )

        return( partialProductionIntegralSum )

    def processCoulombPlusNuclearMuCutoff( self, style, excludeRutherfordScattering = False ) :

        nuclearPlusCoulombInterference = None

        for reaction in self.reactions :
            RutherfordScatteringExcluded = reaction.processCoulombPlusNuclearMuCutoff( style, excludeRutherfordScattering = excludeRutherfordScattering )

            if( RutherfordScatteringExcluded is not None ) :
                crossSection, function2d = RutherfordScatteringExcluded
                firstProduct = reaction.outputChannel[0]
                secondProduct = reaction.outputChannel[1]

                nuclearPlusCoulombInterference = nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference( )

                try :
                    LLNL = self.applicationData['LLNL']
                except :
                    LLNL = institutionModule.institution( 'LLNL' )
                    self.applicationData.add( LLNL )
                LLNL.append( nuclearPlusCoulombInterference )

                nuclearPlusCoulombInterference.reaction.crossSection.add( crossSection )

                Q = reaction.outputChannel.Q[0].copy( )
                Q.label = style.label
                nuclearPlusCoulombInterference.reaction.outputChannel.Q.add( Q )

                product1 = productModule.product( firstProduct.id, firstProduct.label )
                product1.multiplicity.add( firstProduct.multiplicity[0].copy( ) )
                product1.multiplicity[0].label = style.label
                angularTwoBody1 = angularModule.twoBodyForm( crossSection.label, standardsModule.frames.centerOfMassToken, function2d )
                product1.distribution.add( angularTwoBody1 )
                nuclearPlusCoulombInterference.reaction.outputChannel.products.add( product1 )

                product2 = productModule.product( secondProduct.id, secondProduct.label )
                product2.multiplicity.add( secondProduct.multiplicity[0].copy( ) )
                product2.multiplicity[0].label = style.label
                recoil = angularModule.recoil( link = product1.distribution[style.label], relative = True )
                angularTwoBody2 = angularModule.twoBodyForm( crossSection.label, standardsModule.frames.centerOfMassToken, recoil )
                product2.distribution.add( angularTwoBody2 )
                nuclearPlusCoulombInterference.reaction.outputChannel.products.add( product2 )

                nuclearPlusCoulombInterference.reaction.outputChannel.process = 'no Rutherford'
                nuclearPlusCoulombInterference.reaction.updateLabel( )

        return( nuclearPlusCoulombInterference )

    def processThermalNeutronScatteringLaw( self, style, indent = '', verbosity = 0 ) :

        kwargs = { 'reactionSuite' : self }
        kwargs['product'] = self.PoPs[IDsPoPsModule.neutron]
        kwargs['incidentEnergyUnit'] = self.domainUnit
        kwargs['temperatureUnit'] = kwargs['incidentEnergyUnit'] + '/k'
        kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'
        kwargs['momentumUnit'] = kwargs['incidentEnergyUnit'] + '/c'
        kwargs['neutronMass'] = self.PoPs[self.projectile].getMass( kwargs['massUnit'] )
        kwargs['productMass'] = kwargs['neutronMass']
        kwargs['accuracy'] = 0.001
        kwargs['epsilon'] = 1e-6
        kwargs['verbosity'] = verbosity

        energyMax = self.styles.getEvaluatedStyle().projectileEnergyDomain.max
        if( energyMax > 0.0 ) : kwargs['energyMax'] = energyMax

        for reaction in self.reactions : reaction.processThermalNeutronScatteringLaw( style, kwargs )

    def processMC_cdf( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', additionalReactions = [] ) :

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        tempInfo = { 'reactionSuite' : self }
        tempInfo['verbosity'] = verbosity
        tempInfo['incrementalIndent'] = incrementalIndent
        tempInfo['brokenLinks'] = []

        for reaction in self.reactions : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )
        for reaction in self.orphanProducts : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )
        for reaction in additionalReactions : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )

        for original, newReference in tempInfo['brokenLinks']:
            newReference.link = original.referenceInstance.ancestor[ newReference.label ]

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', additionalReactions = [] ) :
        """
        Generate a union grid for all cross sections, spanning the same domain as elastic scattering.
        """

        def mutualDomain( style, reaction, thresholds ) :

            crossSection = style.findFormMatchingDerivedStyle( reaction.crossSection )
            if( isinstance( crossSection, crossSectionModule.regions1d ) ) :
                crossSection = crossSection.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-6 )
            crossSection = crossSection.domainSlice( domainMin, domainMax )                 # truncate to match elastic domain
            if( ( crossSection.domainMin > domainMin ) and ( crossSection[0][1] != 0 ) ) :
                print('    WARNING, domains not mutual, setting first y-value to 0.')
                crossSection[0] = [ crossSection[0][0], 0 ]
            if( ( crossSection.domainMax < domainMax ) and ( crossSection[-1][1] != 0 ) ) :
                print('    WARNING, domains not mutual, setting last y-value to 0.')
                crossSection[-1] = [ crossSection[-1][0], 0 ]
            thresholds.add( crossSection[0][0] )
            return( crossSection )

        def photoAtomic( style, reaction, thresholds ) :

            crossSection = style.findFormMatchingDerivedStyle( reaction.crossSection )
            if( isinstance( crossSection, crossSectionModule.regions1d ) ) :
                values = []
                x1, y1 = crossSection[0][0]
                for iRegion, region in enumerate( crossSection ) :
                    domainGrid = region.domainGrid
                    values += domainGrid
            else :
                x1, y1 = crossSection[0]
                values = crossSection.domainGrid
            if( ( x1 > domainMin ) and ( y1 != 0 ) ) :
                xEpsilon = 1e-6 * x1 + x1
                if( xEpsilon < values[1] ) : values.insert( 1, xEpsilon )

            return( set( values ) )

        evaluationStyle = style.findDerivedFromStyle( stylesModule.evaluated )
        projectileEnergyDomain = evaluationStyle.projectileEnergyDomain
        domainMin = projectileEnergyDomain.min
        domainMax = projectileEnergyDomain.max

        axes = style.findFormMatchingDerivedStyle( self.reactions[0].crossSection ).axes.copy( )

        isPhotoAtomic = False
        if( not( self.isThermalNeutronScatteringLaw( ) ) ) :
            target = self.PoPs[self.target]
            isChemicalElement = isinstance( target, chemicalElementPoPsModule.chemicalElement )
            isUnorthodox = isinstance( target, unorthodoxPoPsModule.particle )
            isPhotoAtomic = ( self.projectile == IDsPoPsModule.photon ) and ( isChemicalElement or isUnorthodox )

        thresholds = set( )
        total = crossSectionModule.XYs1d( data = [], axes = axes )
        if( isPhotoAtomic ) :
            values = set( )
            for reaction in self.reactions : values.update( set( photoAtomic( style, reaction, thresholds ) ) )
            values = sorted( values )
        else :
            for reaction in self.reactions : total += mutualDomain( style, reaction, thresholds )
            for reaction in self.orphanProducts : total += mutualDomain( style, reaction, thresholds )
            values = total.domainGrid

        tooCloseXsIndicies = []
        x1 = -1
        for i1, x2 in enumerate( values ) :
            if( abs( x2 - x1 ) < 1e-12 * x2 ) : tooCloseXsIndicies.append( i1 )
            x1 = x2
        tooCloseXsIndicies.reverse( )
        for index in tooCloseXsIndicies :                   # Remove close points but any threshold points.
            if( values[index] in thresholds ) : index -= 1
            del values[index]

        values = valuesModule.values( values )
        grid = axesModule.grid( total.axes[1].label, 1, total.axes[1].unit, axesModule.pointsGridToken, values )
        style.grid = grid

        for reaction in self.reactions :
            reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent, isPhotoAtomic = isPhotoAtomic )

        if( not( isPhotoAtomic ) ) :
            for reaction in self.orphanProducts :
                reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for production in self.productions :
                production.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for sum in self.sums.crossSectionSums :
                sum.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for fissionComponent in self.fissionComponents :
                fissionComponent.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

        for reaction in additionalReactions :
            reaction.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

    def processMultiGroup( self, style, legendreMax, verbosity = 0, indent = '', incrementalIndent = '  ', logFile = None, workDir = None,
                restart = False, additionalReactions = [] ) :

        from LUPY import times as timesModule

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))
        if( not( isinstance( style, stylesModule.heatedMultiGroup ) ) ) : raise( 'Instance is not a heatedMultiGroup style.' )

        t0 = timesModule.times( )

        kwargs = { 'reactionSuite' : self }
        kwargs['failures'] = 0
        if( len( self.reactions ) > 0 ) :
            kwargs['verbosity'] = verbosity
            kwargs['incrementalIndent'] = incrementalIndent
            kwargs['logFile'] = logFile
            kwargs['incidentEnergyUnit'] = self.domainUnit
            kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'
            kwargs['restart'] = restart
            kwargs['legendreMax'] = legendreMax

            projectile = self.PoPs[self.projectile]
            kwargs['projectile'] = projectile
            kwargs['projectileZA'] = chemicalElementMiscPoPsModule.ZA( projectile )
            kwargs['projectileMass'] = projectile.getMass( kwargs['massUnit'] )

            kwargs['targetMass'] = None
            if( not( self.isThermalNeutronScatteringLaw( ) ) ) :
                target = self.PoPs[self.target]
                kwargs['isInfiniteTargetMass'] = isinstance( target, chemicalElementPoPsModule.chemicalElement )
                if( not( isinstance( target, chemicalElementPoPsModule.chemicalElement ) ) ) :
                    if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
                    kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )
                kwargs['target'] = target
                kwargs['targetZA'] = chemicalElementMiscPoPsModule.ZA( target )

            kwargs['masses'] = { 'Projectile'   : kwargs['projectileMass'],
                                 'Product'      : None,
                                 'Residual'     : None,
                                 'Target'       : kwargs['targetMass'] }

            kwargs['zeroPerTNSL'] = self.isThermalNeutronScatteringLaw( )
            projectileEnergyDomain = style.projectileEnergyDomain
            transportable = style.transportables[self.projectile]
            for index, boundary in enumerate( transportable.group.boundaries.values ) :
                if( boundary > projectileEnergyDomain.max ) : break
                kwargs['maximumProjectileGroupIndex'] = index

            if( workDir is None ) : workDir = 'Merced.work'
            kwargs['workDir'] = workDir
            kwargs['workFile'] = []

            style.processMultiGroup( style, kwargs, indent + incrementalIndent )
            kwargs['groupedFlux'] = [ x for x in style.multiGroupFlux.array.constructArray( )[:,0] ]

# BRB FIXME, must have test to determine if reconstructResonances is needed.
#        self.reconstructResonances( styleName = 'reconstructed', accuracy = 1e-3, verbose = False )
            for i1, reaction in enumerate( self.reactions ) :
                kwargs['reactionIndex'] = "%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.orphanProducts ) :
                kwargs['reactionIndex'] = "o%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.productions ) :
                kwargs['reactionIndex'] = "p%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.sums.crossSectionSums ) :
                kwargs['reactionIndex'] = "s%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( self.fissionComponents ) :
                kwargs['reactionIndex'] = "f%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )
            for i1, reaction in enumerate( additionalReactions ) :
                kwargs['reactionIndex'] = "a%.4d" % i1
                reaction.processMultiGroup( style, kwargs, indent + incrementalIndent )

        if logFile is not None: logFile.write( '    ' + str( t0 ) + '\n' )
        if( kwargs['failures'] > 0 ) : raise ValueError( "kwargs['failures'] = %d  > 0" % kwargs['failures'] )

    def processSnElasticUpScatter( self, style, legendreMax, verbosity = 0, indent = '', incrementalIndent = '  ', logFile = None ) :

        if( self.target == IDsPoPsModule.neutron ) : return

        if( not( isinstance( style, stylesModule.SnElasticUpScatter ) ) ) : raise TypeError("style must be an instance of stylesModule.SnElasticUpScatter, not %s" % type(style))

        tempInfo = { 'reactionSuite' : self }
        tempInfo['verbosity'] = verbosity
        tempInfo['incrementalIndent'] = incrementalIndent
        tempInfo['logFile'] = logFile
        incidentEnergyUnit = self.domainUnit
        tempInfo['temperature'] = style.temperature
        tempInfo['legendreOrder'] = legendreMax

        tempInfo['minEval'] = PQUModule.PQU( 1.e-11, 'MeV' ).getValueAs( incidentEnergyUnit )
        tempInfo['maxEval'] = PQUModule.PQU( 20., 'MeV' ).getValueAs( incidentEnergyUnit )

        massUnit = incidentEnergyUnit + '/c**2'
        target = self.PoPs[self.target]
        if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
        tempInfo['targetMassRatio'] = target.getMass( massUnit ) / self.PoPs[self.projectile].getMass( massUnit )
        tempInfo['groupBoundaries'] = style.transportables[self.projectile].group
        tempInfo['zeroPerTNSL'] = False

        elasticReaction = self.getReaction( 'elastic' )
        product = elasticReaction.outputChannel.getProductWithName( IDsPoPsModule.neutron )

        TM_1, TM_E, averageEnergy, maxIncidentGroup = multiGroupUpScatterModule.SnElasticUpScatter( style, tempInfo, comment = '' )
        style.upperCalculatedGroup = maxIncidentGroup

                # multiply constant cross section into the distribution (upscatter code assumes sigma(E) = 1 )
        groupedCrossSec = product.findAttributeInAncestry( 'crossSection' )[style.derivedFromStyle.label].array.values
        TM_1 = multiGroupUpScatterModule.rescaleCrossSection( groupedCrossSec, TM_1 )

        distribution = product.distribution                             # Put matrix into GNDS structures.
        derivedForm = distribution[style.derivedFrom]
        multiGroupSubform = derivedForm.multiGroupSubform
        array = multiGroupSubform.array.constructArray( )
        for incidentGroup in range( style.upperCalculatedGroup + 1, len( array ) ) :
            if( incidentGroup not in TM_1 ) : TM_1[incidentGroup] = {}
            for outgoingGroup in range( len( array[incidentGroup] ) ) :
                if( outgoingGroup not in TM_1[incidentGroup] ) : TM_1[incidentGroup] = {}
                TM_1[incidentGroup][outgoingGroup] = list( array[incidentGroup][outgoingGroup] )

        tempInfo['productName'] = product.id
        tempInfo['groupedFlux'] = [ x for x in style.derivedFromStyle.multiGroupFlux.array.constructArray( )[:,0] ]
        tempInfo['incidentEnergyUnit'] = incidentEnergyUnit
        tempInfo['crossSection']  = style.derivedFromStyle.findFormMatchingDerivedStyle( elasticReaction.crossSection )
        multiGroup = groupModule.TMs2Form( style, tempInfo, TM_1, TM_E )
        distribution.add( multiGroup )

        energyDeposition = product.energyDeposition
        axes = energyDepositionModule.defaultAxes( energyUnit = incidentEnergyUnit )
        averageEnergy = energyDepositionModule.XYs1d( data = averageEnergy, axes = axes, label = style.label )
        averageEnergy = averageEnergy.processMultiGroup( style, tempInfo, indent )
        grouped = averageEnergy.array.constructArray( )
        array = energyDeposition[style.derivedFrom].array.constructArray( )
        for incidentGroup in range( style.upperCalculatedGroup + 1, len( array ) ) : grouped[incidentGroup] = array[incidentGroup]
        energyDeposition.add( groupModule.toMultiGroup1d( energyDepositionModule.gridded1d, style, tempInfo, averageEnergy.axes, grouped, zeroPerTNSL = tempInfo['zeroPerTNSL'] ) )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        self.__styles.removeStyles( styleLabels )
        self.__reactions.removeStyles( styleLabels )
        self.__orphanProducts.removeStyles( styleLabels )
        self.__productions.removeStyles( styleLabels )
        self.__incompleteReactions.removeStyles( styleLabels )
        self.__fissionComponents.removeStyles( styleLabels )
        self.__sums.removeStyles( styleLabels )

    def supportsResonanceReconstruction( self ):
        """
        Checks whether self.resonances contains any resonance parameters that can be reconstructed into
        cross sections and/or product distributions.
        :return: boolean
        """
        if self.resonances is None: return False
        return self.resonances.reconstructCrossSection

    def reconstructResonances( self, style, accuracy = None, thin = True, significantDigits = None, verbose = False ):
        """
        Turn resonance parameters into pointwise cross sections, then merge the results with
        tabulated pointwise cross sections. Resulting pointwise cross sections are stored
        alongside the original 'resonancesWithBackground' data in the reactionSuite.

        @:param style: style - an instance of style crossSectionReconstructed.
        @:param accuracy: float - target accuracy during reconstruction. For example, 0.001.
        @:param thin: boolean - enable/disable thinning after reconstruction (disabling makes summed cross sections consistency checks easier).
        @:param significantDigits: int - energy grid will only include points that can be exactly represented using specified number of significant digits.
        @:param verbose: boolean - turn on/off verbosity
        """

        if( self.resonances is None ) : return
        if not self.resonances.reconstructCrossSection:
            return # nothing to do
        from fudge.processing.resonances import reconstructResonances

        if not isinstance( style, stylesModule.crossSectionReconstructed ):
            raise TypeError("style must be an instance of crossSectionReconstructed, not %s" % type(style))

        xsecs = reconstructResonances.reconstructResonances(self, tolerance = accuracy, verbose = verbose,
                significantDigits = significantDigits, energyUnit = self.domainUnit )
        epsilon = 1e-8  # for joining multiple regions together

        # for each reaction, add tabulated pointwise data (ENDF MF=3) to reconstructed resonances:
# BRB6 hardwired?
        possibleChannels = { 'elastic' : False, 'capture' : True, 'fission' : True, 'total' : False, 'nonelastic' : False }
        derivedFromLabel = ''
        for reaction in self :
            if isinstance( reaction, sumsModule.multiplicitySum ): continue
            evaluatedCrossSection = reaction.crossSection.evaluated
            if not isinstance( evaluatedCrossSection, crossSectionModule.resonancesWithBackground ):
                continue
            # which reconstructed cross section corresponds to this reaction?
            if( derivedFromLabel == '' ) : derivedFromLabel = evaluatedCrossSection.label
            if( derivedFromLabel != evaluatedCrossSection.label ) :
                print ('WARNING derivedFromLabel = "%s" != "%s"' % (derivedFromLabel, evaluatedCrossSection.label))
            RRxsec = None
            if str( reaction ) in xsecs:
                RRxsec = xsecs[ str( reaction ) ]
            else :
                for possibleChannel in possibleChannels :
                    if( possibleChannels[possibleChannel] ) :
                        if( possibleChannel in str( reaction ) ) : RRxsec = xsecs[possibleChannel]
                    if( RRxsec is None ) :
                        if( reaction is self.getReaction( possibleChannel ) ) : RRxsec = xsecs[possibleChannel]
                    if( RRxsec is not None ) : break
            if( RRxsec is None ) :
                if verbose:
                    print( "Warning: couldn't find appropriate reconstructed cross section to add to reaction %s" % reaction )
                continue

            background = evaluatedCrossSection.background
            try:
                background = background.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = epsilon, upperEps = epsilon )
            except Exception as ex:
                print("Encountered an error while linearizing background cross section for reaction '%s'!" % reaction.label)
                raise ex
            RRxsec = RRxsec.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = epsilon, upperEps = epsilon )
            RRxsec.convertUnits( {RRxsec.domainUnit: background.domainUnit,  RRxsec.rangeUnit: background.rangeUnit } )

            background, RRxsec = background.mutualify(0,0,0, RRxsec, -epsilon,epsilon,True)
            RRxsec = background + RRxsec    # result is a crossSection.XYs1d instance
            if( RRxsec.rangeMin < 0 ) :
                # turn any negative xsc to 0
                RRxsec = RRxsec.clip( rangeMin=0 )
                if verbose:
                    print( "Warning: negative cross section encountered for %s; changed to 0 b" % reaction )
            if( RRxsec[0][1] == 0 ) :
                # special handling for threshold reaction:
                #  trim leading zeros + points with cross section < 1e-90 b,
                #  thin aggressively up to xsc = 1e-10 b
                firstNonZero = next((i for i, xy in enumerate(RRxsec) if xy[1]), None)
                RRxsec = RRxsec[firstNonZero-1:]

                points = RRxsec.copyDataToXYs()
                if points[0][0] > evaluatedCrossSection.domainMin:
                    points.insert(0, [evaluatedCrossSection.domainMin, 0])

                indices = []
                priorCrossSection = 0.0
                factor = 1e5
                for idx, xy in enumerate(points):
                    if idx==0: continue
                    if xy[1] > 1e-10: break
                    if xy[1] > 1e-20: factor = 10
                    if xy[1] < 1e-90 or xy[1] < factor * priorCrossSection:
                        indices.append(idx)
                    else:
                        priorCrossSection = xy[1]

                for idx in reversed(indices):
                    del points[idx]

                RRxsec.setData(points)

            if thin:
                RRxsec = RRxsec.thin( accuracy or .01 )
            RRxsec.label = style.label
            reaction.crossSection.add( RRxsec )

        self.styles.add( style )

    def reconstructResonancesAngularDistributions( self, styleName, overwrite=False, accuracy = None, thin = False, verbose = False ):
        """
        Turn resonance parameters into Legendre-expanded angular distributions, then merge the results with
        tabulated angular distributions (overwriting original values in the resonance region). Resulting pointwise
        angular distributions are stored alongside the original angular distributions in the reactionSuite.

        @:param styleName: string - label for reconstructed distribution style
        @:param overwrite: boolean - if style already in use for another distribution, overwrite it with the reconstructed distribution
        @:param accuracy: float - giving target accuracy during reconstruction. For example, 0.001
        @:param thin: boolean - enable/disable thinning after reconstruction
        @:param verbose: boolean - turn on/off verbosity
        """

        if accuracy: raise NotImplementedError("Refining interpolation grid for angular distributions still TBD")
        if thin: raise NotImplementedError("Thinning for angular distributions still TBD")

        if( self.resonances is None ) : return
        from fudge.processing.resonances import reconstructResonances
        from fudge.productData.distributions import angular as angularModule

        newStyle = stylesModule.angularDistributionReconstructed( styleName, 'eval' )       # BRB6, FIXME - 'eval' should not be hardwired.

        distributions = reconstructResonances.reconstructAngularDistributions( self, tolerance=accuracy, verbose=verbose )

        for key in distributions:
            reaction = self.getReaction( key )
            # Kludgy way to discern the product's name so we can look it up
# BRB6 hardwired.
            if key in 'elastic': productName = 'n'
            else:
# BRB6 hardwired.
                pairs = [ ( p.particle.getMass('amu'), p.id ) for p in reaction.outputChannel.products ]
                productName = min(pairs)[1] # pick's the lightest product
            product = reaction.outputChannel.getProductWithName(productName)
            original = product.distribution.evaluated.angularSubform
            reconstructed = distributions[key]

            merged = angularModule.regions2d( axes = original.axes )     # Caleb, FIXME - check that this is correct.
            merged.append( reconstructed )

            if isinstance( original, angularModule.XYs2d ):
                newregion = angularModule.XYs2d( axes = original.axes, interpolation=original.interpolation,
                        interpolationQualifier=original.interpolationQualifier )
                newregion.append( original.getAtEnergy( reconstructed.domainMax ) )
                for val in original:
                    if( val.outerDomainValue <= reconstructed.domainMax ) : continue
                    newregion.append( val )
                merged.append( newregion )
            elif isinstance( original, angularModule.regions2d ):
                for idx,region in enumerate(original):
                    if( region.domainMax > reconstructed.domainMax ) : break
                if( region.domainMin != reconstructed.domainMax ) :
                    newregion = angularModule.XYs2d( axes = region.axes, interpolation=region.interpolation,
                            interpolationQualifier=region.interpolationQualifier )
                    newregion.append( region.getAtEnergy( reconstructed.domainMax ) )
                    for val in region:
                        if( val.value <= reconstructed.domainMax ) : continue
                        newregion.append( val )
                    merged.append( newregion )
                    idx += 1
                for region in original[idx:]:
                    merged.append( region )

            newForm = angularModule.twoBodyForm( label = newStyle.label,
                    productFrame = product.distribution.evaluated.productFrame, angularSubform = merged )
            if overwrite and newStyle.label in product.distribution:
                product.distribution.remove( newStyle.label )
            product.distribution.add( newForm )
        if newStyle.label not in self.styles: self.styles.add( newStyle )

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
                if product.outputChannel is None: continue
                for dproduct in product.outputChannel:
                    removeStyleFromComponent( style, dproduct.energyDeposition )
                    removeStyleFromComponent( style, dproduct.momentumDeposition )
        self.styles.remove( style )

    def saveToOpenedFile( self, fOut, **kwargs ) :

        xmlString = self.toXMLList( **kwargs )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, **kwargs ) :
        """
        Save the reactionSuite in GNDS/xml structure to specified file.
        To suppress extra messages, change the 'verbosity' flag:
        >>>self.saveToFile( "output.xml", flags={'verbosity':0} )
        """

        dirname = os.path.dirname( fileName )
        if( ( len( dirname ) > 0 ) and not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        with open( fileName, "w" ) as fout :
# BRB6 hardwired.
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n' )
            self.saveToOpenedFile( fout, **kwargs )

    def saveToHybrid( self, xmlName, hdfName=None, minLength=4, flatten=False, compress=False, **kwargs ):
        """
        Save the reactionSuite to a hybrid layout, with the data hierarchy in xml
        but most actual data saved in an associated HDF file

        :param xmlName: file name to save.
        :param hdfName: if None, will be the same as xmlName with extension changed to .h5
        :param minLength: minimum number of values before moving a dataset from XML to HDF
            Need to determine best default setting
        :param flatten: if True, GNDS datasets are concatenated into flattened HDF5 datasets
        :param compress: enable gzip + shuffle compression for HDF5 datasets
        :param kwargs:
        :return:
        """
        import numpy
        import h5py
        from . import externalFile as externalFileModule
        from LUPY import checksums

        if hdfName is None:
            if not xmlName.endswith('.xml'):
                hdfName = xmlName + '.h5'
            else:
                hdfName = xmlName.replace('.xml','.h5')

        with h5py.File(hdfName, "w") as h5:
            HDF_opts = {
                'h5file': h5,
                'index': 0,     # assign unique name to each dataset
                'minLength': minLength,
                'flatten': flatten,
                'iData': [],
                'dData': [],
                'compression': {}
            }
            if compress:
                HDF_opts['compression'] = {'compression':'gzip', 'shuffle':'true', 'chunks':True}

            # add filler external file, actual checksum computed below
            relHdfName = hdfName
            if relHdfName == os.path.abspath(relHdfName):
                relHdfName = os.path.basename(relHdfName)
            self.externalFiles.add(
                    externalFileModule.externalFile( "HDF", relHdfName,
                        checksum = "deadbeef", algorithm="sha1" ) )

            xmlString = self.toXMLList( HDF_opts = HDF_opts, **kwargs )

            if len(HDF_opts['iData']) > 0:
                iData = numpy.array( HDF_opts['iData'], dtype=numpy.int32 )
                h5.create_dataset('iData', data=iData, **HDF_opts['compression'])
            if len(HDF_opts['dData']) > 0:
                dData = numpy.array( HDF_opts['dData'], dtype=numpy.float64 )
                h5.create_dataset('dData', data=dData, **HDF_opts['compression'])

        sha1sum = checksums.sha1sum.from_file(hdfName)
        for idx, line in enumerate(xmlString):
            if 'externalFile' in line and 'checksum="deadbeef"' in line:
                xmlString[idx] = line.replace("deadbeef", sha1sum)
                break

        with open( xmlName, "w" ) as fout :
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n')
            fout.write( '\n'.join( xmlString ) )
            fout.write( '\n' )

        self.externalFiles.pop("HDF")


    def thermalNeutronScatteringLawTemperatures( self ) :

        temperatures = {}
        for reaction in self.reactions : reaction.thermalNeutronScatteringLawTemperatures( temperatures )
        return( temperatures )

    def toXML( self, indent = "", **kwargs ) :
        """Returns an GNDS/XML string representation of self."""

        return( '\n'.join( self.toXMLList( **kwargs ) ) )

    def toXMLList( self, indent = "", **kwargs ) :
        """Returns a list of GNDS/XML strings representing self."""

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent
        formatVersion = kwargs.get( 'formatVersion', self.formatVersion )
        kwargs['formatVersion'] = formatVersion
        if( formatVersion not in formatVersionModule.allowed ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )

        xmlString = [ '%s<%s projectile="%s" target="%s" evaluation="%s" format="%s" projectileFrame="%s" interaction="%s">'
            % ( indent, self.moniker, self.projectile, self.target, self.evaluation, formatVersion, self.projectileFrame, self.interaction ) ]

        xmlString += self.externalFiles.toXMLList(indent2, **kwargs)
        xmlString += self.styles.toXMLList( indent2, **kwargs )

        if( formatVersion == formatVersionModule.version_1_10 ) :
            # documentation used to live directly under reactionSuite
            name = None
            _documentation = self.styles.getEvaluatedStyle().documentation
            if( _documentation.body.filled ) :
                name = 'ENDL'
                text = _documentation.body.body
            elif( _documentation.endfCompatible.filled ) :
                name = 'endfDoc'
                text = _documentation.endfCompatible.body
            if( name is not None ) :
                from . import documentation as oldDocumentationModule

                xmlString.append( '%s<documentations>' % indent2 )
                documentation = oldDocumentationModule.documentation( name, text )
                xmlString += documentation.toXMLList( indent3, **kwargs )
                xmlString[-1] += '</documentations>'

        xmlString += self.PoPs.toXMLList( indent2, **kwargs )

        if self.resonances is not None:
            xmlString += self.resonances.toXMLList( indent2, **kwargs )

        xmlString += self.reactions.toXMLList( indent2, **kwargs )
        xmlString += self.orphanProducts.toXMLList( indent2, **kwargs )
        xmlString += self.sums.toXMLList( indent2, **kwargs )
        xmlString += self.fissionComponents.toXMLList( indent2, **kwargs )
        xmlString += self.productions.toXMLList( indent2, **kwargs )
        xmlString += self.incompleteReactions.toXMLList(indent2, **kwargs)
        xmlString += self.applicationData.toXMLList(indent2, **kwargs)

        xmlString.append( '%s</%s>' % ( indent, self.moniker ) )
        return( xmlString )

    def toString( self, indent = '' ) :
        """Returns a string representation of an reactionSuite."""

        s = indent + self.inputParticlesToReactionString( suffix = ' --> ' )
        indent2 = '\n' + len( s ) * ' '
        indent3 = ''
        for reaction in self.reactions :
            s += reaction.toString( indent = indent3 )
            indent3 = indent2
        return( s )

def readXML( gndsFile, numberOfBrokenLinksToPrint = 5 ) :
    """
    Read a GNDS/xml file and create a new reactionSuite instance from the result.

    :param gndsFile: string path to a GNDS file or a file object.
    :return: reactionSuite instance containing all data from the file.
    """

    from xml.etree import cElementTree
    # wrapper around the xml parser:
    from LUPY.xmlNode import xmlNode

    rsElement = cElementTree.parse( gndsFile ).getroot()
    rsElement = xmlNode( rsElement, xmlNode.etree )
    return parseXMLNode( rsElement, sourcePath = gndsFile, numberOfBrokenLinksToPrint = numberOfBrokenLinksToPrint )

def parseXMLNode( rsElement, **kwargs ):
    """Translates a <reactionSuite> xml node into a reactionSuite instance. Users should use the 'readXML' function instead."""

    xPath = ['reactionSuite']   # Keep track of location in the tree, in case errors come up.
    try :
        sourcePath = kwargs.get( 'sourcePath' )
        if( not( isinstance( sourcePath, str ) ) ) : sourcePath = sourcePath.name
        path = os.path.dirname( sourcePath )
        numberOfBrokenLinksToPrint = kwargs.get( 'numberOfBrokenLinksToPrint', 5 )

        formatVersion = rsElement.get( 'format' )
        kwargs['formatVersion'] = formatVersion

        projectile = rsElement.get( 'projectile' )
        target = rsElement.get( 'target' )
        evaluation = rsElement.get( 'evaluation' )
        projectileFrame = rsElement.get( 'projectileFrame' )
        interaction = rsElement.get( 'interaction' )
        guessedInteraction, interaction = guessInteractionModule.guessInteraction(interaction, projectile, target)

        rs = reactionSuite( projectile, target, evaluation, formatVersion = formatVersion, projectileFrame = projectileFrame, interaction = interaction, sourcePath = sourcePath )

        if formatVersion == formatVersionModule.version_1_10:
            rs.guessedInteraction = guessedInteraction
        elif guessedInteraction:
            raise Exception('Function guessInteraction had to guess at the interaction attribute for a non GNDS %s file.' % formatVersion)

        linkData = { 'reactionSuite' : rs, 'unresolvedLinks' : [], 'format' : formatVersion }

        externalFiles = rsElement.find( suitesModule.externalFiles.moniker )
        if externalFiles is not None:
            rs.externalFiles.parseXMLNode( externalFiles, xPath, linkData )

        if "HDF" in rs.externalFiles:   # FIXME hard-coded label
            h5File = os.path.join( path, rs.externalFiles["HDF"].path )
            if not os.path.exists(h5File):
                raise IOError("External HDF file %s not found" % h5File)

            import h5py
            h5File = h5py.File( h5File, 'r', libver='latest' )
            linkData["HDF"] = {
                "h5File": h5File
                }
            for key in ("iData", "dData"):
                if key in h5File:
                    linkData[key] = h5File[key][()]

        styles = rsElement.find( 'styles' )
        if( styles is not None ) : rs.styles.parseXMLNode( styles, xPath, linkData )
        PoPs = rsElement.find( 'PoPs' )
        if( PoPs is not None ) : rs.PoPs.parseXMLNode( PoPs.data, xPath, linkData )

        oldDocumentationElement = None

        for child in rsElement :
            if( child.tag in ( 'styles', 'PoPs', suitesModule.externalFiles.moniker ) ) :
                continue    # already read above
            elif( child.tag == resonancesModule.resonances.resonances.moniker ) :
                rs.resonances = resonancesModule.resonances.resonances.parseXMLNode(child, xPath, linkData)
            elif( child.tag == rs.reactions.moniker ) :
                rs.reactions.parseXMLNode( child, xPath, linkData )
            elif( child.tag == rs.orphanProducts.moniker ) :
                rs.orphanProducts.parseXMLNode( child, xPath, linkData )
            elif( child.tag == rs.sums.moniker ) :
                rs.sums.parseXMLNode( child, xPath, linkData )
            elif( child.tag == rs.productions.moniker ) :
                rs.productions.parseXMLNode( child, xPath, linkData)
            elif( child.tag == rs.incompleteReactions.moniker ) :
                rs.incompleteReactions.parseXMLNode( child, xPath, linkData)
            elif( child.tag == rs.fissionComponents.moniker ) :
                rs.fissionComponents.parseXMLNode( child, xPath, linkData )
            elif( child.tag == rs.applicationData.moniker ) :
                rs.applicationData.parseXMLNode( child, xPath, linkData )
            else :
                if( ( formatVersion == formatVersionModule.version_1_10 ) and ( child.tag == 'documentations' ) ) :
                    oldDocumentationElement = child
                else :
                    print( "Warning: encountered unexpected element '%s' in reactionSuite!" % child.tag )

        if( oldDocumentationElement is not None ) :
            if( len( oldDocumentationElement ) > 1 ) :
                print( 'Too many children of old documentations node: number of children = %s.' % len( oldDocumentationElement ) )
            else :
                oldDocumentationElement = oldDocumentationElement[0]
                text = oldDocumentationElement.text.lstrip( '\n' )
                if( len( text.split( '\n' ) ) < 2 ) :
                    rs.styles.getEvaluatedStyle().documentation.body.body = text
                else :
                    rs.styles.getEvaluatedStyle().documentation.endfCompatible.body = text

    except Exception :
        print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
        raise

    # Fix links.
    unresolvedLinksCounter = 0
    for quant in linkData['unresolvedLinks'] :
        try :
            if quant.path.startswith( '/' + reactionSuite.moniker ) :
                quant.link = quant.follow( rs )
            elif quant.path.startswith( '/covarianceSuite' ) :
                rs._externalLinks.append( quant )
            else:   # relative link
                quant.link = quant.follow( quant )
        except ancestryModule.XPathNotFound :
            if( unresolvedLinksCounter < numberOfBrokenLinksToPrint ) : print( '    Cannot resolve link "%s".' % quant )
            unresolvedLinksCounter += 1
        except :
            raise

    if( unresolvedLinksCounter > 0 ) : print( 'WARNING: %s unresolved links found while parsing file "%s".' % ( unresolvedLinksCounter, rs.sourcePath ) )

    return( rs )
