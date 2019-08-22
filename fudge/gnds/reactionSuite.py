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

from __future__ import print_function

"""
reactionSuite.py contains the 'reactionSuite' class that in turn holds all reactions for a given target/projectile.
reactionSuite is the top-level class for the GNDS structure.
"""

import re
import os

from PoPs import IDs as IDsPoPsModule
from PoPs import database as PoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.groups import chemicalElement as chemicalElementPoPsModule
from PoPs.decays import misc as decaysMiscPoPsModule

from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import standards as standardsModule

from . import suites as suitesModule
from . import sums as sumsModule
from . import styles as stylesModule

from .reactionData import crossSection as crossSectionModule
from .reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule

from .version import GNDS_VERSION

import fudge

__metaclass__ = type

class nullDevice :
    """For internal use. Used in methods to set logFile when one is not entered."""

    def write( self, **kwargs ) : pass

class reactionSuite( ancestryModule.ancestry ) :
    """
    This is the main class for a gnds projectile/target object. It contains
        * documentation
        * a list of all particles encountered in the file
        * resonance data
        * a list of reactions
    """

    moniker = 'reactionSuite'
    ancestryMembers = ( 'styles', 'reactions', 'orphanProducts', 'sums', 'productions', 'fissionComponents' )

    def __init__( self, projectile, target, evaluation, GNDS_version = GNDS_VERSION, style = None,
            documentation = None, projectileFrame = standardsModule.frames.labToken, MAT = None, PoPs = None ) :
        """
        Creates a new reactionSuite object, containing all reaction data for a projectile/target combination.

        :param projectile: projectile ID (string)
        :param target: target ID (string)
        :param evaluation: the evaluation for this protare (string)
        :param GNDS_version: GNDS format version (string)
        :param style: Indicates style of data in this reactionSuite (evaluated, linearized, processed, etc.).  Instance of styles.style
        :param documentation: Top-level documentation for the entire reactionSuite (instance of documentation.documentation)
        :param projectileFrame: frame of the projectile. Allowed values are 'lab' or 'centerOfMass'
        :param MAT: integer ENDF MAT number (only needed for writing back to ENDF-6)
        :param PoPs: a PoPs database that is copied to self's PoPs database
        :return:
        """

        ancestryModule.ancestry.__init__( self )

        if( GNDS_version not in ( GNDS_VERSION, ) ) : raise Exception( "Unsupported GNDS structure '%s'!" % str(GNDS_version) )

        if( not( isinstance( projectile, str ) ) ) : raise TypeError( 'projectile ID not a string' )
        self.projectile = projectile
        self.projectileFrame = projectileFrame
        self.MAT = MAT

        if( not( isinstance( target, str ) ) ) : raise TypeError( 'target ID not a string' )
        self.target = target

        if( not( isinstance( evaluation, str ) ) ) : raise TypeError( 'evaluation not a string' )
        self.__evaluation = evaluation

        self.GNDS_version = GNDS_version

        if( PoPs is not None ) :
            self.__PoPs = PoPs.copy( )
        else :
            self.__PoPs = PoPsModule.database( 'protare_internal', '1.0' )
        self.__PoPs.setAncestor( self )

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )
        if( style is not None ) : self.styles.add( style )

        self.__externalFiles = suitesModule.externalFiles()
        self.__externalFiles.setAncestor( self )

        self.resonances = None

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

        self.documentation = {}
        if( documentation is not None ) : self.addDocumentation( documentation )

        self._externalLinks = []    # keep track of links to external files

    def __iter__( self ) :

        for reaction in self.reactions : yield reaction
        for orphanProduct in self.orphanProducts : yield orphanProduct
        for _sum in self.sums.crossSections : yield _sum
        for reaction in self.fissionComponents : yield reaction
        for reaction in self.productions : yield reaction
        for reaction in self.incompleteReactions : yield reaction

    def getProductionReactions( self ) :

        return( self.productions )

    def __str__( self ) :
        """See method toString."""

        return( self.toString( ) )

    @property
    def compound_nucleus(self):
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
    def evaluation( self ) :
        """Returns self's evaluation."""

        return( self.__evaluation )

    @property
    def styles(self):
        return self.__styles

    @property
    def externalFiles(self):
        return self.__externalFiles

    @property
    def reactions(self):
        return self.__reactions

    @property
    def orphanProducts(self):
        return self.__orphanProducts

    @property
    def sums(self):
        return self.__sums

    @property
    def productions(self):
        return self.__productions

    @property
    def incompleteReactions(self):
        return self.__incompleteReactions

    @property
    def fissionComponents(self):
        return self.__fissionComponents

    @property
    def applicationData(self):
        return self.__applicationData

    @property
    def PoPs( self ) :

        return( self.__PoPs )

    def addDocumentation( self, documentation ) :

        self.documentation[documentation.name] = documentation

    def addResonances( self, resonances ):
        """Add resonance parameter data (gnds.resonances.resonances instance) to the reactionSuite. """
        self.resonances = resonances
        self.resonances.setAncestor( self )

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

        from fudge.gnds import warning

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
                        'momentumDepositionUnit': incidentEnergyUnit + '/c',
                        'projectileMass': projectileMass, 'targetMass': targetMass,
                        'reactionSuite':self, 'photonBranchingData': self.photonBranchingData() }

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
                    if isinstance( elastic_distribution, fudge.gnds.productData.distributions.angular.XYs2d):
                        linearized = elastic_distribution.toPointwise_withLinearXYs()
                        forward_scattering = crossSectionModule.XYs1d( data=[], axes=crossSectionModule.XYs1d.defaultAxes() )
                        for energy_in in linearized.getEnergyArray('eV'):
                            forward_scattering.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                    elif isinstance( elastic_distribution, fudge.gnds.productData.distributions.angular.regions2d):
                        forward_scattering = crossSectionModule.regions1d( axes=crossSectionModule.XYs1d.defaultAxes() )
                        for region in elastic_distribution:
                            ptw = crossSectionModule.XYs1d( data=[] )
                            linearized = region.toPointwise_withLinearXYs()
                            for energy_in in linearized.getEnergyArray('eV'):
                                ptw.setValue( energy_in, linearized.interpolateAtValue( energy_in ).evaluate(1.0) )
                            forward_scattering.append( ptw )
                        forward_scattering = forward_scattering.toPointwise_withLinearXYs( accuracy = 1e-5, 
                                lowerEps = 1e-8, upperEps = 1e-8 )

                    mutualDomain = zip( *[dat.domain() for dat in (elastic_xsc, total_xsc, forward_scattering)])
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

        for reaction in self :
            reactionWarnings = reaction.check( info )
            if reactionWarnings: warnings.append( warning.context('%s label %s'
                % (reaction.moniker, reaction.label), reactionWarnings ) )

        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info

        if options['cleanUpTempFiles']:
            if info['reconstructedStyleName'] is not None: self.removeStyle( info['reconstructedStyleName'] )

        return result

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

    def hasParticle( self, name ) :

        return( name in self.PoPs )

    def heatCrossSections( self, style, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.001, heatAllPoints = False,
            doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True, setThresholdToZero = False ) :

        EMin = self.reactions[0].crossSection.domainMin
        for reaction in self.reactions :
            reaction.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit, 
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero )

        for orphanProduct in self.orphanProducts :
            orphanProduct.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit, 
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero )

        for production in self.productions :
            production.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit, 
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero )

        for _sum in self.sums.crossSections :
            _sum.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit,
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero )

        for fissionComponent in self.fissionComponents :
            fissionComponent.heatCrossSection( style, EMin, lowerlimit = lowerlimit, upperlimit = upperlimit, 
                    interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                    heatBelowThreshold = heatBelowThreshold, heatAllEDomain = heatAllEDomain, setThresholdToZero = setThresholdToZero )

    def getProjectileFrame( self ) :

        return( self.projectileFrame )

    def getDocumentation( self, name ) :

        return( self.documentation[name] )

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

        # If channel is an int, it might be an ENDF MT, so check those first
        if type(channel)==int: # must be an ENDF MT
            for reactionList in self.reactions, self.sums.crossSections, self.fissionComponents, self.productions:
                for reaction in reactionList:
                    if reaction.ENDF_MT == channel: return reaction
            return None

        # translate special channel names:
        if channel=='elastic': channel = channel_tr = '%s + %s' % (self.projectile,self.target)
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
            retVal  = [ r1 for r1 in self.reactions         if r1.outputChannel.fissionGenre == genre ]
            retVal += [ r1 for r1 in self.fissionComponents if r1.outputChannel.fissionGenre == genre ]
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

    def getTemperature( self, style ) :

        return( self.styles[style].temperature )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """
        Calculate average energy and momentum data for all products of all reactions.
        Resulting data are stored within each product. Example usage is:

        from fudge.gnds import reactionSuite as reactionSuiteModule
        from fudge.gnds import styles as stylesSuiteModule
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

        kwargs['incidentEnergyUnit'] = self.reactions[0].crossSection.domainUnit
        kwargs['momentumDepositionUnit'] = kwargs['incidentEnergyUnit'] + '/c'
        kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'
        kwargs['projectileMass'] = self.PoPs[self.projectile].getMass( kwargs['massUnit'] )
        target = self.PoPs[self.target]
        if( isinstance( target, chemicalElementPoPsModule.chemicalElement ) ) :
            kwargs['targetMass'] = None
        else :
            if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
            kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        kwargs['reactionSuite'] = self

        for reaction in self.reactions :
            reaction.calculateAverageProductData( style, indent = indent2, **kwargs )
        for reaction in self.orphanProducts :
            reaction.calculateAverageProductData( style, indent = indent2, **kwargs )

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

        CoulombPlusNuclearElastic = None
        for reaction in self.reactions :
            for doubleDifferentialCrossSection in reaction.doubleDifferentialCrossSection :
                if( isinstance( doubleDifferentialCrossSection, CoulombPlusNuclearElasticModule.form ) ) :
                    CoulombPlusNuclearElastic = doubleDifferentialCrossSection
                    break
            if( CoulombPlusNuclearElastic is not None ) : break

        if( CoulombPlusNuclearElastic is None ) : return( None )

        muCutoffs = []
        if( CoulombPlusNuclearElastic.nuclearPlusInterference is not None ) :
            muCutoffs.append( CoulombPlusNuclearElastic.nuclearPlusInterference.muCutoff )

        for style in self.styles :
            if( isinstance( style, stylesModule.CoulombPlusNuclearElasticMuCutoff ) ) : muCutoffs.append( style.muCutoff )

        return( list( set( muCutoffs ) ) )

    def processCoulombPlusNuclearMuCutoff( self, style ) :

        for reaction in self.reactions : reaction.processCoulombPlusNuclearMuCutoff( style )

    def processMC_cdf( self, style, verbosity = 0, indent = '', incrementalIndent = '  ' ) :

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))

        tempInfo = { 'reactionSuite' : self }
        tempInfo['verbosity'] = verbosity
        tempInfo['incrementalIndent'] = incrementalIndent
        tempInfo['brokenLinks'] = []

        for reaction in self.reactions : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )
        for reaction in self.orphanProducts : reaction.processMC_cdf( style, tempInfo, indent + incrementalIndent )

        for original, newReference in tempInfo['brokenLinks']:
            newReference.link = original.referenceInstance.ancestor[ newReference.label ]

    def processGriddedCrossSections( self, style, verbosity = 0, indent = '', incrementalIndent = '  ' ) :
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

        domainMin = min( [ reaction.domainMin for reaction in self.reactions ] )
        domainMax = max( [ reaction.domainMax for reaction in self.reactions ] )

        axes = style.findFormMatchingDerivedStyle( reaction.crossSection ).axes.copy( )

        target = self.PoPs[self.target]
        isChemicalElement = isinstance( target, chemicalElementPoPsModule.chemicalElement )
        isPhotoAtomic = ( self.projectile == IDsPoPsModule.photon ) and isChemicalElement

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
                
            for sum in self.sums.crossSections :
                sum.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

            for fissionComponent in self.fissionComponents :
                fissionComponent.processGriddedCrossSections( style, verbosity = verbosity, indent = indent, incrementalIndent = incrementalIndent )

    def processMultiGroup( self, style, verbosity = 0, indent = '', incrementalIndent = '  ', logFile = None, workDir = None ) :

        from fudge.core.utilities import times as timesModule

        if( verbosity > 0 ) : print ('%s%s' % (indent, self.inputParticlesToReactionString(suffix=" -->")))
        if( not( isinstance( style, stylesModule.heatedMultiGroup ) ) ) : raise( 'Instance is not a heatedMultiGroup style.' )

        t0 = timesModule.times( )

        kwargs = { 'reactionSuite' : self }
        kwargs['failures'] = 0
        if( len( self.reactions ) > 0 ) :
            kwargs['verbosity'] = verbosity
            kwargs['incrementalIndent'] = incrementalIndent
            kwargs['logFile'] = logFile
            kwargs['incidentEnergyUnit'] = self.reactions[0].crossSection.domainUnit
            kwargs['massUnit'] = kwargs['incidentEnergyUnit'] + '/c**2'

            projectile = self.PoPs[self.projectile]
            kwargs['projectile'] = projectile
            kwargs['projectileZA'] = chemicalElementMiscPoPsModule.ZA( projectile )
            kwargs['projectileMass'] = projectile.getMass( kwargs['massUnit'] )

            target = self.PoPs[self.target]
            if( isinstance( target, chemicalElementPoPsModule.chemicalElement ) ) :
                kwargs['targetMass'] = None
            else :
                if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
                kwargs['targetMass'] = target.getMass( kwargs['massUnit'] )
            kwargs['target'] = target
            kwargs['targetZA'] = chemicalElementMiscPoPsModule.ZA( target )

            kwargs['masses'] = {
                'Projectile' : kwargs['projectileMass'],
                'Product' : None,
                'Residual' : None }
            kwargs['masses']['Target'] = kwargs['targetMass']
            if( workDir is None ) : 'Merced.work'
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

        logFile.write( str( t0 ) + '\n' )
        if( kwargs['failures'] > 0 ) : raise ValueError( "kwargs['failures'] = %d  > 0" % kwargs['failures'] )

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

        @:param styleName: string - label for reconstructed cross section style
        @:param accuracy: float - target accuracy during reconstruction. For example, 0.001
        @:param thin: boolean - enable/disable thinning after reconstruction
                Disabling thinning makes it easier to check for consistency of summed cross sections.
        @:param significantDigits: int - energy grid will only include points that can be exactly represented using
                specified number of significant digits
        @:param verbose: boolean - turn on/off verbosity
        """

        if( self.resonances is None ) : return
        if not self.resonances.reconstructCrossSection:
            return # nothing to do
        from fudge.processing.resonances import reconstructResonances

        if not isinstance( style, stylesModule.crossSectionReconstructed ):
            raise TypeError("style must be an instance of crossSectionReconstructed, not %s" % type(style))

        xsecs = reconstructResonances.reconstructResonances(self, tolerance = accuracy, verbose = verbose,
                significantDigits = significantDigits)
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
            background = background.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = epsilon, upperEps = epsilon )
            RRxsec = RRxsec.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = epsilon, upperEps = epsilon )
            RRxsec.convertUnits( {RRxsec.domainUnit: background.domainUnit,  RRxsec.rangeUnit: background.rangeUnit } )

            background, RRxsec = background.mutualify(0,0,0, RRxsec, -epsilon,epsilon,True)
            RRxsec = background + RRxsec    # result is a crossSection.XYs1d instance
            if( RRxsec.rangeMin < 0 ) :
                # turn any negative xsc to 0
                RRxsec = RRxsec.clip( rangeMin=0 )
                if verbose:
                    print( "Warning: negative cross section encountered for %s; changed to 0 b" % reaction )
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
        from fudge.gnds.productData.distributions import angular as angularModule

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
                    if( val.value <= reconstructed.domainMax ) : continue
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
        """Save the reactionSuite in GNDS/xml structure to specified file.
        To suppress extra messages, change the 'verbosity' flag:
            >>>self.saveToFile( "output.xml", flags={'verbosity':0} )
        """

        dirname = os.path.dirname( fileName )
        if( ( len( dirname ) > 0 ) and not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        with open( fileName, "w" ) as fout :
# BRB6 hardwired.
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n' )
            self.saveToOpenedFile( fout, **kwargs )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( **kwargs ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xmlString = [ '%s<%s projectile="%s" target="%s" evaluation="%s" format="%s" projectileFrame="%s">'
            % ( indent, self.moniker, self.projectile, self.target, self.evaluation, self.GNDS_version, self.projectileFrame ) ]

        xmlString += self.styles.toXMLList( indent2, **kwargs )
        xmlString += self.externalFiles.toXMLList(indent2, **kwargs)

# BRB6 documentation should have its own toXMLList.
        xmlString.append( '%s<documentations>' % indent2 )
        for doc in self.documentation : xmlString += self.documentation[doc].toXMLList( indent3, **kwargs )
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
        indent2 = len( s ) * ' '
        indent3 = ''
        for channel in self.reactions :
            s += channel.toString( indent = indent3 )
            indent3 = indent2
        return( s )

def readXML( gndsFile ):
    """
    Read a GNDS/xml file and create a new reactionSuite instance from the result.

    :param gndsFile: path to a GNDS file, as a string.
    :return: reactionSuite instance containing all data from the file.
    """

    from xml.etree import cElementTree
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode

    rsElement = cElementTree.parse( gndsFile ).getroot()
    rsElement = xmlNode( rsElement, xmlNode.etree )
    return parseXMLNode( rsElement )

def parseXMLNode( rsElement ):
    """Translates a <reactionSuite> xml node into a reactionSuite instance. Users should use the 'readXML' function instead."""

    xPath = ['reactionSuite']   # Keep track of location in the tree, in case errors come up.
    try :
        GNDS_version = rsElement.get( 'format' )

        projectile = rsElement.get( 'projectile' )
        target = rsElement.get( 'target' )
        evaluation = rsElement.get( 'evaluation' )
        projectileFrame = rsElement.get( 'projectileFrame' )

        rs = reactionSuite( projectile, target, evaluation, GNDS_version = GNDS_version, projectileFrame = projectileFrame )

        linkData = { 'reactionSuite' : rs, 'unresolvedLinks' : [], 'format' : GNDS_version }
        rs.styles.parseXMLNode( rsElement.find('styles'), xPath, linkData )
        rs.PoPs.parseXMLNode( rsElement.find('PoPs').data, xPath, linkData )
        documents = rsElement.find( 'documentations' )
        if( documents is not None ) :
            for doc in documents :
                rs.addDocumentation( fudge.gnds.documentation.documentation.parseXMLNode( doc, xPath, linkData))

        for child in rsElement :
            if( child.tag in ( 'styles', 'documentations', 'PoPs' ) ) :
                continue    # already read above
            elif(child.tag == suitesModule.externalFiles.moniker) :
                rs.externalFiles.parseXMLNode( child, xPath, linkData )
            elif(child.tag == fudge.gnds.resonances.resonances.moniker) :
                rs.addResonances( fudge.gnds.resonances.resonances.parseXMLNode( child, xPath, linkData))
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
                print( "Warning: encountered unexpected element '%s' in reactionSuite!" % child.tag )
    except Exception :
        print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
        raise

    # Fix links.
    for quant in linkData['unresolvedLinks'] :
        if quant.path.startswith( '/' + reactionSuite.moniker ) :
            quant.link = quant.follow( rs )
        elif quant.path.startswith( '/covarianceSuite' ) :
            rs._externalLinks.append( quant )
        else:   # relative link
            quant.link = quant.follow( quant )

    return( rs )
