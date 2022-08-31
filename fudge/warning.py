# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge import styles as stylesModule

"""
Store and report warnings and errors in GNDS data, in order to discover problems in the library.
The reactionSuite.check() function returns a nested list of Warning objects:

    >>> warnings = reactionSuite.check()
    >>> print(warnings)

May include or exclude specific classes of Warning using the filter command.
filter() returns a new context instance:

    >>> warnings2 = warnings.filter( exclude=[Warning.EnergyImbalance] )

Or, for easier searching you may wish to flatten the list (to get warnings alone without context messages):

    >>> flat = warnings.flatten()
"""

class Context:
    """
    Store warnings in Context. This class contains location information (reactionSuite, reaction, etc)
    plus a nested list of warnings or other Context instances
    """
    def __init__(self, message='', warningList=None):
        self.message = message
        self.warningList = warningList or []

    def __len__(self):
        return len( self.warningList )

    def __getitem__(self, idx):
        return self.warningList[idx]

    def __str__(self):
        if len(self.warningList)==0:
            return self.message + ": no problems encountered"
        return '\n'.join( self.toStringList() )

    def __eq__(self, other):
        return self.message == other.message and self.warningList == other.warningList

    def filter(self, include=None, exclude=None):
        """Filter warning list to only include (or exclude) specific classes of Warning. For example:

            >>> newWarnings = warnings.filter( exclude=[Warning.EnergyImbalance, Warning.Q_mismatch] )

        Note that if both 'include' and 'exclude' lists are provided, exclude is ignored."""

        if include is None and exclude is None: return self
        newWarningList = []
        for warning in self.warningList:
            if isinstance( warning, Context ):
                newContext = warning.filter( include, exclude )
                if newContext: newWarningList.append( newContext )
            elif include is not None:
                if warning.__class__ in include:
                    newWarningList.append( warning )
            else: # exclude is not None:
                if warning.__class__ not in exclude:
                    newWarningList.append( warning )
        return Context( self.message, newWarningList )

    def flatten(self):
        """From a nested hierarchy of warnings, get back a flat list for easier searching:

            >>> w = reactionSuite.check()
            >>> warningList = w.flatten()"""

        List = []
        for val in self.warningList:
            if isinstance(val, Warning):
                List.append(val)
            else:
                List += val.flatten()
        return List

    def toStringList( self, indent='', dIndent='    ' ):
        """ Format warnings for printing. Returns a list of warning strings with indentation. """
        s = ['%s%s' % (indent, self.message)]
        for warning in self.warningList:
            s += warning.toStringList( indent+dIndent )
        return s

class DiffResult :

    def __init__( self, code, message, link1, link2 ) :

        self.code = code
        self.message = message
        self.link1 = link1
        self.link2 = link2

    def __repr__( self ) :

        return( "%s: %s\n    %s\n    %s" % ( self.code, self.message, self.link1, self.link2 ) )

class DiffResults :

    def __init__( self, protare1, protare2 ) :

        self.protare1 = protare1
        self.protare2 = protare2
        self.diffs = []

    def __len__( self ) :
        """Returns the number of differences found."""

        return( len( self.diffs ) )

    def __iter__( self ) :

        n1 = len( self.diffs )
        for i1 in range( n1 ) : yield self.diffs[i1]

    def __repr__( self ) :
        """Returns a string representation of the differences found."""

        return( '\n'.join( [ str( diff ) for diff in self.diffs ] ) )

    def __copy__( self ) :
        """Make a copy of self."""

        diffResults1 = DiffResults( self.protare1, self.protare2 )
        for item in self : diffResults1.append( DiffResult( item.code, item.message, item.link1, item.link2 ) )

        return( diffResults1 )

    def append( self, code, message, link1, link2 ) :
        """Adds a DiffResult instance to self."""

        self.diffs.append( DiffResult( code, message, link1, link2 ) )

    def cull( self, code ) :
        """Removes all DiffResult instance that contain the specified code."""

        cullList = []
        for i1, item in enumerate( self.diffs ) :
            if( item.code == code ) : cullList.append( i1 )
        cullList.reverse( )
        for i1 in cullList : del self.diffs[i1]

    def patch( self, label, version, library = None ) :

        chains = self.protare2.styles.chains( True )
        if( len( chains ) != 1 ) : raise Exception( "Only one styles chain is supported." )

        chain = chains[0]
        for style in chain :
            if( isinstance( style, stylesModule.Evaluated ) ) : break
        if( not( isinstance( style, stylesModule.Evaluated ) ) ) : raise Exception( "At least one evaluated style must be present in file to patch." )

        if( len( self.diffs ) > 0 ) :
            if( label in self.protare2.styles ) : raise ValueError( "Label already in file to patch." )
            if( library is None ) : library = style.library
            evaluateStyle = stylesModule.Evaluated( label, style.label, None, None, library, version )
            if( style.library == evaluateStyle.library ) :
                if( style.version == evaluateStyle.version ) : raise ValueError( "Library and version must be different than one in patch file." )
            self.protare2.styles.add( evaluateStyle )

            for diff in self :
                if( diff.code == "Distribution unspecified - 1" ) :
                    pass
                elif( diff.code == "Distribution unspecified - 2" ) :
                    fromComponent = self.protare1.followXPath( diff.link1 )
                    toComponent = self.protare2.followXPath( diff.link2 )
                    for key in toComponent.keys( ) : toComponent.remove( key )
                    toComponent.add( fromComponent.evaluated )
                elif( diff.code == "reactions missing" ) :
                    if( diff.link1 != "" ) :
                        fromReaction = self.protare1.followXPath( diff.link1 )
                        fromReaction.amendForPatch( style.label, evaluateStyle.label )
                        self.protare2.reactions.add( fromReaction )
                else :
                    print( """Unsupported diff code "%s".""" % diff.code )

        crossSectionReconstructedStyles = self.protare2.styles.getStylesOfClass( stylesModule.CrossSectionReconstructed )
        if( len( crossSectionReconstructedStyles ) > 0 ) :
            if( len( crossSectionReconstructedStyles ) > 1 ) : raise Exception( "Multiple cross section reconstructed styles present which is not supported." )
            crossSectionReconstructedStyle = crossSectionReconstructedStyles[0]
            self.protare2.removeStyle( crossSectionReconstructedStyle.label )
            crossSectionReconstructedStyle = stylesModule.CrossSectionReconstructed( crossSectionReconstructedStyle.label, label )
            self.protare2.reconstructResonances(crossSectionReconstructedStyle, 1e-3)

class Warning(BaseException):
    """
    General Warning class. Contains link to problem object,
    xpath in case the object leaves memory,
    and information about the warning or error.
    """
    def __init__(self, obj=None):
        self.obj = obj
        self.xpath = ''
        if hasattr( obj, 'toXLink' ):
            self.xpath = obj.toXLink()

    def __str__(self):
        return "Generic warning for %s" % self.xpath

    def __eq__(self, other):
        return self.xpath == other.xpath

    def toStringList( self, indent='' ):
        return ['%sWARNING: %s' % (indent, self)]

#
# specific Warning classes:
#

class NotImplemented( Warning ):
    def __init__(self, form, obj=None):
        Warning.__init__(self, obj)
        self.form = form

    def __str__(self):
        return "Checking not yet implemented for %s type data" % self.form

    def __eq__(self, other):
        return (self.form == other.form and self.xpath == other.xpath)

class UnorthodoxParticleNotImplemented(Warning):
    def __init__(self, particleId, obj=None):
        Warning.__init__(self, obj)
        self.particleId = particleId

    def __str__(self):
        return "Checking not yet implemented for unorthodox particle '%s'" % self.particleId

    def __eq__(self, other):
        return (self.particleId == other.particleId)


class ElementalTarget(Warning):
    def __init__(self, particleId, obj=None):
        Warning.__init__(self, obj)
        self.particleId = particleId

    def __str__(self):
        return "Elemental target '%s'!  Some checks (e.g. energy balance and Q-values) will be skipped" % self.particleId

    def __eq__(self, other):
        return (self.particleId == other.particleId)


# external files:

class MissingExternalFile(Warning):
    def __init__(self, path):

        Warning.__init__(self)
        self.path = path

    def __str__(self):
        return "External file '%s' is missing" % self.path

class WrongExternalFileChecksum(Warning):

    def __init__(self, path, expected, computed):

        Warning.__init__(self)
        self.path = path
        self.expected = expected
        self.computed = computed

    def __str__(self):
        s = "Computed checksum for external file '%s' does not match: %s vs. %s" % (
            self.path, self.computed, self.expected
        )
        return s

# resonance region:

class BadScatteringRadius(Warning):

    def __init__(self, factor=3.0, gotAP=None, expectedAP=None, L=None, E=None ):

        Warning.__init__(self)
        self.factor=factor
        self.gotAP=gotAP
        self.expectedAP=expectedAP
        self.L=L
        self.E=E

    def __str__(self):
        direction = 'high'
        if self.gotAP<self.expectedAP: direction = 'low'
        s = "Scattering radius (%s) at least a factor of %s too %s compared to expectations (%s)" % ( list( map( str,( self.gotAP,self.factor,direction,self.expectedAP ) ) ) )
        if self.L is not None: s+=" for L=%i"%self.L
        if self.E is not None: s+=" for E=%s"%str(self.E)
        return s

class BadSpinStatisticalWeights(Warning):

    def __init__(self, L, gJ, expectedgJ, reaction=None):

        Warning.__init__(self)
        self.L=L
        self.gJ=gJ
        self.expectedgJ=expectedgJ
        self.reaction=reaction

    def __str__(self):
        direction = 'many'
        if self.gJ<self.expectedgJ: direction = 'few'
        s="The spin statical weights for L=%i sums to %s, but should sum to %s.  You have too %s channels " % ( self.L, str(self.gJ), str(self.expectedgJ), direction )
        if self.reaction is not None:
            s+='for reaction '+self.reaction
        return s

class InvalidAngularMomentaCombination(Warning):

    def __init__(self,L,S,J,name):

        Warning.__init__(self)
        self.L=L
        self.S=S
        self.J=J
        self.name=name

    def __str__(self):
        return 'Invalid angular momenta combination: cannot couple L = %s and S = %s up to J = %s for channel "%s"' % (str(self.L),str(self.S),str(self.J),self.name)

class InvalidParity(Warning):

    def __init__(self,L,S,J,Pi,name):

        Warning.__init__(self)
        self.L=L
        self.S=S
        self.J=J
        self.Pi=Pi
        self.name=name

    def __str__(self):
        return 'Invalid parity %s for channel "%s" with L=%s, S=%s, J=%s' % (str(self.Pi),self.name,str(self.L),str(self.S),str(self.J))

class UnknownSpinParity(Warning):

    def __init__(self,reaction):

        Warning.__init__(self)
        self.reaction = reaction

    def __str__(self):
        return "Could not determine product spin/parities for reaction '%s'" % self.reaction

class MissingResonanceChannel(Warning):

    def __init__(self,L,S,J,name):

        Warning.__init__(self)
        self.L=L
        self.S=S
        self.J=J
        self.name=name

    def __str__(self):
        return 'Missing a channel with angular momenta combination L = %s, S = %s and J = %s for "%s"' % (
            self.L, self.S, self.J, self.name)

class InvalidSpinCombination(Warning):

    def __init__(self,Ia,Ib,S,name):

        Warning.__init__(self)
        self.Ia=Ia
        self.Ib=Ib
        self.S=S
        self.name=name

    def __str__(self):
        return 'Invalid spin combination: cannot couple Ia = %s and Ib = %s up to S = %s for channel "%s"' % (str(self.Ia),str(self.Ib),str(self.S),self.name)

class PotentialScatteringNotConverged(Warning):

    def __init__(self, L, E, fom, fomTarget):

        Warning.__init__(self)
        self.L=L
        self.E=E
        self.fom=fom
        self.fomTarget=fomTarget

    def __str__(self):
        s="Potential scattering hasn't converged by L=%i at E=%s eV, xs[%i]/xs[0]=%s > %s"%\
          (self.L,str(self.E),self.L,str(100.*self.fom)+'%',str(100.*self.fomTarget)+'%')
        return s


class RRmultipleRegions(Warning):

    def __init__(self, obj=None):
        Warning.__init__(self, obj)

    def __str__(self):
        return "Use of more than one resolved resonance regions is deprecated"

class URRdomainMismatch(Warning):

    def __init__(self, Lval, Jval, obj=None):

        Warning.__init__(self, obj)
        self.Lval = Lval
        self.Jval = Jval

    def __str__(self):
        return "Unresolved L=%i / J=%.1f widths don't span URR energy limits" % (self.Lval, self.Jval)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval)

class URRinsufficientEnergyGrid(Warning):

    def __init__(self, Lval, Jval, eLow, eHigh, obj=None):

        Warning.__init__(self, obj)
        self.Lval = Lval
        self.Jval = Jval
        self.eLow = eLow
        self.eHigh = eHigh

    def __str__(self):
        return "More points needed in L=%i J=%.1f unresolved widths between %s and %s" % (self.Lval, self.Jval, self.eLow, self.eHigh)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval and self.eLow == other.eLow and self.eHigh == other.eHigh)

class UnresolvedLink(Warning):

    def __init__(self,link):

        Warning.__init__(self, obj)
        self.link=link

    def __str__(self):
        return "Unresolved link to %s" % (str(self.link))

# reaction objects:

class WicksLimitError(Warning):

    def __init__(self, percentErr, energy_in, obj=None):

        Warning.__init__(self, obj)
        self.percentErr = percentErr
        self.energy_in = energy_in

    def __str__(self):
        return "Wick's limit too low by %.3f%% at %s" % (100*self.percentErr, self.energy_in)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.percentErr == other.percentErr and self.energy_in == other.energy_in)

class ZAbalanceWarning(Warning):

    def __str__(self):
        return "ZA doesn't balance for this reaction!"

class Q_mismatch(Warning):

    def __init__(self, Qcalc, Qactual, obj=None):

        Warning.__init__(self, obj)
        self.Qcalc = Qcalc
        self.Qactual = Qactual

    def __str__(self):
        return "Calculated and tabulated Q-values disagree: %s vs %s!" % (self.Qcalc, self.Qactual)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.Qcalc == other.Qcalc and self.Qactual == other.Qactual)

class BadFissionEnergyRelease(Warning):

    def __init__(self, worstCase, obj=None):

        Warning.__init__(self, obj)
        self.worstCase = worstCase

    def __str__(self):
        return "Fission energy release seems unphysical. Worst case: %s" % (self.worstCase)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.worstCase == other.worstCase)

class NegativeDelayedRate(Warning):

    def __init__(self, rate, obj=None):

        Warning.__init__(self, obj)
        self.rate = rate

    def __str__(self):
        return "Negative delayed neutron decay constant: %s" % (self.rate)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.rate == other.rate)

class Threshold_mismatch(Warning):

    def __init__(self, threshold, thresholdCalc, obj=None):

        Warning.__init__(self, obj)
        self.threshold = threshold
        self.thresholdCalc = thresholdCalc

    def __str__(self):
        return "Calculated and tabulated thresholds disagree: %s vs %s!" % (self.thresholdCalc, self.threshold)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.threshold == other.threshold
                and self.thresholdCalc == other.thresholdCalc)

class Coulomb_threshold_mismatch(Threshold_mismatch):

    def __str__(self):
        return "Tabulated threshold %s below calculated threshold %s!" % (self.threshold, self.thresholdCalc)

class NonZero_crossSection_at_threshold(Warning):

    def __init__(self, value, obj=None):

        Warning.__init__(self, obj)
        self.value = value

    def __str__(self):
        return "First cross section point for threshold reaction should be 0, not %s" % self.value

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.value == other.value)

class GapInCrossSection(Warning):

    def __init__(self, minGapEnergy, maxGapEnergy, obj=None):

        Warning.__init__(self, obj)
        self.minGapEnergy = minGapEnergy
        self.maxGapEnergy = maxGapEnergy

    def __str__(self):
        return "Gap in cross section data from %s to %s" % (self.minGapEnergy, self.maxGapEnergy)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.minGapEnergy == other.minGapEnergy
                and self.maxGapEnergy == other.maxGapEnergy)

class NegativeCrossSection(Warning):

    def __init__(self, energy_in, index, obj=None):

        Warning.__init__(self, obj)
        self.energy_in = energy_in
        self.index = index

    def __str__(self):
        return "Negative cross section encountered at %s (index %i)!" % (self.energy_in, self.index)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index)

class BadCrossSectionReference(Warning):

    def __str__(self):
        return "Cross section reference outside of current reactionSuite not allowed!"

class SummedCrossSectionMismatch(Warning):

    def __init__(self, maxPercentDiff, obj=None):

        Warning.__init__(self, obj)
        self.maxPercentDiff = maxPercentDiff

    def __str__(self):
        msg = ("Cross section does not match sum of linked reaction cross sections! Max diff: %.2f%%" %
                self.maxPercentDiff)
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.maxPercentDiff == other.maxPercentDiff)

class SummedCrossSectionDomainMismatch(Warning):

    def __str__(self):
        return "Cross section domain does not match domain of linked reaction cross section sum"

class SummedCrossSectionZeroDivision(Warning):

    def __str__(self):
        return "Divide by 0 error when comparing summed cross section to components!"

class SummedMultiplicityMismatch(Warning):

    def __init__(self, maxPercentDiff, obj=None):

        Warning.__init__(self, obj)
        self.maxPercentDiff = maxPercentDiff

    def __str__(self):
        msg = ("Multiplicity does not match sum of linked product multiplicities! Max diff: %.2f%%" %
                self.maxPercentDiff)
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.maxPercentDiff == other.maxPercentDiff)

class SummedMultiplicityDomainMismatch(Warning):
    def __str__(self):
        return "Multiplicity domain does not match domain of linked product multiplicity sum"

class SummedMultiplicityZeroDivision(Warning):
    def __str__(self):
        return "Divide by 0 error when comparing summed multiplicity to components!"

class NegativeMultiplicity(Warning):
    def __init__(self, value, obj=None):
        Warning.__init__(self, obj)
        self.value = value

    def __str__(self):
        return "Encountered negative multiplicity (%s)!" % self.value

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.value == other.value)

class Domain_mismatch(Warning):
    def __init__(self, lowBound, highBound, xscLowBound, xscHighBound, obj=None):
        Warning.__init__(self, obj)
        self.lowBound, self.highBound = lowBound, highBound
        self.xscLowBound, self.xscHighBound = xscLowBound, xscHighBound

    def __str__(self):
        return "Domain doesn't match the cross section domain: (%s -> %s) vs (%s -> %s)" % (self.lowBound,
                self.highBound, self.xscLowBound, self.xscHighBound)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.lowBound == other.lowBound
                and self.highBound == other.highBound)

class MissingDistribution(Warning):
    def __init__(self, productName, obj=None):
        Warning.__init__(self, obj)
        self.productName = productName

    def __str__(self):
        return "Missing distribution (required for all '%s' products)!" % self.productName

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.productName == other.productName)

""" eliminate this one?
class NoDistributions(Warning):
    # more serious than missingDistributions:
    def __str__(self):
        return "No distribution data available for *any* of this reaction's products!"
"""

class MissingRecoilDistribution(Warning):
    def __str__(self):
        return "Recoil distribution type specified, but recoil partner has unsupported distribution type!"

class WrongDistributionComponent(Warning):
    def __init__(self, component, reactionType='2-body', obj=None):
        Warning.__init__(self, obj)
        self.component = component
        self.reactionType = reactionType

    def __str__(self):
        return "%s is not a valid distribution component for %s reaction" % (self.component, self.reactionType)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.component == other.component
                and self.reactionType == other.reactionType)

class Wrong2BodyFrame(Warning):
    def __str__(self):
        return ("2-body reaction not in center-of-mass frame!")

class UncorrelatedFramesMismatch(Warning):
    def __init__(self, angleFrame, energyFrame, obj=None):
        Warning.__init__(self, obj)
        self.angleFrame = angleFrame
        self.energyFrame = energyFrame

    def __str__(self):
        return ("For uncorrelated energy/angle distributions, frame must be identical! Currently, angle given in %s, but energy in %s"
                % (self.angleFrame, self.energyFrame))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.angleFrame == other.angleFrame
                and self.energyFrame == other.energyFrame)

class FlatIncidentEnergyInterpolation(Warning):
    def __str__(self):
        return ("For distributions, flat interpolation along incident energy is unphysical!")

class EnergyDistributionBadU(Warning):
    def __str__(self):
        return ("For energy distribution functional form, parameter 'U' results in negative outgoing energies!")

class NegativeDiscreteGammaEnergy(Warning):
    def __str__(self):
        return ("Discrete gamma energy <= 0")

class PrimaryGammaEnergyTooLarge(Warning):
    """Primary gamma energy at threshold should be <= available energy (depending on which discrete level it ends up in)"""
    def __init__(self, energy, fraction, obj=None):
        Warning.__init__(self, obj)
        self.energy = energy
        self.fraction = fraction

    def __str__(self):
        return ( "Primary gamma energy %s exceeds available energy by %s%%"
                 % (self.energy, self.fraction) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy == other.energy
                and self.fraction == other.fraction)

class MadlandNixBadParameters(Warning):
    def __init__(self, EFL, EFH, minTm, obj=None):
        Warning.__init__(self, obj)
        self.EFL = EFL
        self.EFH = EFH
        self.minTm = minTm

    def __str__(self):
        return ( "Madland-Nix fission spectrum contains negative parameters: EFL=%s, EFH=%s, min(Tm)=%s"
                % (self.EFL, self.EFH, self.minTm) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.EFL == other.EFL
                and self.EFH == other.EFH and self.minTm == other.minTm)

class WeightsDontSumToOne(Warning):
    def __str__(self):
        return ("Weights don't sum to 1.0!")

class UnnormalizedDistribution(Warning):
    def __init__(self, energy_in, index, integral, obj=None):
        Warning.__init__(self, obj)
        self.energy_in = energy_in
        self.index = index
        self.integral = integral

    def __str__(self):
        return "Unnormalized distribution! At energy_in = %s (index %i), integral = %s" % (self.energy_in,
                self.index, self.integral)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index and self.integral == other.integral)

class UnnormalizedDistributionAtMu(Warning):
    """
    Only appears when checking legacy ENDL data (I=1 and 3)
    """
    def __init__(self, mu, integral, obj=None):
        Warning.__init__(self, obj)
        self.mu = mu
        self.integral = integral

    def __str__(self):
        return "Unnormalized distribution! At mu = %s, integral = %s" % (self.mu, self.integral)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.mu == other.mu
                and self.integral == other.integral)

class UnnormalizedKMDistribution(UnnormalizedDistribution):
    pass    # identical to UnnormalizedDistribution, except it needs a special fix() function

class KalbachMannDomainMismatch(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, obj)

    def __str__(self):
        return "Kalbach-Mann terms do not span the same incident energy domain!"

class IncompleteDistribution(Warning):
    def __init__(self, energy_in, lowerMu, upperMu, obj=None):
        Warning.__init__(self, obj)
        self.energy_in = energy_in
        self.lowerMu = lowerMu
        self.upperMu = upperMu

    def __str__(self):
        return ("Incomplete angular coverage (mu = %.2f to %.2f) for distribution at incident energy %s"
                % (self.lowerMu, self.upperMu, self.energy_in))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.lowerMu == other.lowerMu and self.upperMu == other.upperMu)

class NegativeProbability(Warning):

    def __init__(self, energy_in, energy_out=None, mu=None, value=None, obj=None):
        Warning.__init__(self, obj)
        self.energy_in = energy_in
        self.energy_out = energy_out
        self.mu = mu
        self.value = value

    def __str__(self):
        msg = "Negative probabilities encountered. Incident energy: %s" % self.energy_in
        if self.mu is not None: msg += ", mu: %s" % self.mu
        if self.energy_out is not None: msg += ", outgoing energy: %s" % self.energy_out
        if self.value is not None: msg += ", worst case: %s" % self.value
        if self.obj is not None: msg += ', %s' % self.obj.moniker
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.energy_out == other.energy_out and self.mu == other.mu)

class ExtraOutgoingEnergy(Warning):
    """If an outgoing energy distribution ends with more than one energy with probability=0,
    proper unitbase treatment is unclear. Distribution should end with exactly one P=0 point."""

    def __init__(self, energy_in, obj=None):
        Warning.__init__(self, obj)
        self.energy_in = energy_in

    def __str__(self):
        msg = "Extra zero-probability outgoing energies found at incident energy %s" % self.energy_in
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in)

class MissingCoulombIdenticalParticlesFlag(Warning):
    def __init__(self, obj=None):
        Warning.__init__(self, obj)

    def __str__(self):
        return "Need 'identicalParticles=\"true\"' when target==projectile"

class IncorrectCoulombIdenticalParticlesFlag(Warning):
    def __init__(self, projectile, target, obj=None):
        Warning.__init__(self, obj)
        self.projectile = projectile
        self.target = target

    def __str__(self):
        return "doubleDifferentialCrossSection claims that target (%s) & projectile (%s) are identical!"\
               % (self.target, self.projectile)

    def __eq__(self, other):
        return (self.target == other.target and self.projectile == other.projectile)

class EnergyImbalance(Warning):
    def __init__(self, energy_in, index, availableEnergy, deposition_per_product, obj=None):
        Warning.__init__(self, obj)
        self.energy_in = energy_in
        self.index = index
        self.availableEnergy = availableEnergy
        self.deposition_per_product = deposition_per_product
        self.total_deposited = sum( [val[1] for val in deposition_per_product] )
        if self.availableEnergy == 0: self.total_deposited = float('inf')

    def __str__(self):
        per_product = ', '.join(["%s = %.4g%%" % (key,val) for key,val in self.deposition_per_product[:5]])
        if len( self.deposition_per_product ) > 5: per_product += ', ...'
        return ("Energy imbalance at incident energy %s (index %i). Total deposited = %.4g%% (%s)" %
                (self.energy_in, self.index, self.total_deposited, per_product))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index and self.availableEnergy == other.availableEnergy
                and self.deposition_per_product == other.deposition_per_product
                and self.total_deposited == other.total_deposited)

class FissionEnergyImbalance(EnergyImbalance):
    def __str__(self):
        per_product = ', '.join(["%s = %.4g%%" % (key,val) for key,val in self.deposition_per_product[:5]])
        if len( self.deposition_per_product ) > 5: per_product += ', ...'
        return ("Fission energy imbalance at incident energy %s (index %i). Total deposited = %.4g%% (%s), leaving insufficient energy for fission products!" %
                (self.energy_in, self.index, self.total_deposited, per_product))

class ValueOutOfRange(Warning):
    def __init__(self, contextMessage, value, lowerBound, upperBound, obj=None):
        Warning.__init__(self, obj)
        self.contextMessage = contextMessage
        self.value = value
        self.lowerBound = lowerBound
        self.upperBound = upperBound

    def __str__(self):
        return ("%s. Value=%s, should be in range %s -> %s" % (self.contextMessage, self.value,
            self.lowerBound, self.upperBound) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.contextMessage == other.contextMessage
                and self.value == other.value and self.lowerBound == other.lowerBound
                and self.upperBound == other.upperBound)

class TestSkipped(Warning):
    """ indicate if test was skipped due to missing information """
    def __init__(self, testName, reason, obj=None):
        Warning.__init__(self, obj)
        self.testName = testName
        self.reason = reason

    def __str__(self):
        return ("Skipped test %s: %s" % (self.testName, self.reason))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.testName == other.testName and self.reason == other.reason)

class ExceptionRaised(Warning):

    """If we run into an exception when running check(), try to exit gracefully and return this warning."""

    def __init__(self, Exception_String, obj=None):
        Warning.__init__(self, obj)
        self.Exception_String = Exception_String

    def __str__(self):
        return ("Encountered Exception: %s" % self.Exception_String)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.Exception_String == other.Exception_String)

class EnergyDepositionExceptionRaised( ExceptionRaised ):
    def __str__(self):
        return ("Exception raised when calculating energy deposition: %s" % self.Exception_String)

class SkippedCoulombElasticEnergyDeposition(Warning):
    def __str__(self):
        return ("Energy/momentum deposition cannot be computed for charged-particle elastic")

### covarianceSuite warnings: ###

class CyclicDependency(Warning):
    def __init__(self, cycle, obj=None):
        Warning.__init__(self, obj)
        self.cycle = tuple(cycle)

    def __str__(self):
        if len(self.cycle) == 2:
            return ("Cyclic dependency in summed covariances for sections %s and %s" % self.cycle)
        else:
            return ("Cyclic dependency in summed covariances for section %s" % self.cycle)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.cycle == other.cycle)

class VarianceTooSmall(Warning):
    def __init__(self, index, variance, obj=None):
        Warning.__init__(self, obj)
        self.index = index
        self.variance = variance

    def __str__(self):
        return ("Encountered very small variance (%e%%) at index %i." % (self.variance, self.index))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.variance == other.variance)

class VarianceTooLarge(Warning):
    def __init__(self, index, variance, obj=None):
        Warning.__init__(self, obj)
        self.index = index
        self.variance = variance

    def __str__(self):
        return ("Encountered very large variance (%e%%) at index %i." % (self.variance, self.index))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.variance == other.variance)

class NegativeEigenvalues(Warning):
    def __init__(self, negativeCount, worstCase, obj=None):
        Warning.__init__(self, obj)
        self.negativeCount = negativeCount
        self.worstCase = worstCase

    def __str__(self):
        return ("%i negative eigenvalues! Worst case = %e" % (self.negativeCount, self.worstCase) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.negativeCount == other.negativeCount
                and self.worstCase == other.worstCase)

class BadEigenvalueRatio(Warning):
    def __init__(self, ratio, obj=None):
        Warning.__init__(self, obj)
        self.ratio = ratio

    def __str__(self):
        return ("Ratio of smallest/largest eigenvalue (%e) is too small" % (self.ratio) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.ratio == other.ratio)

class InvalidShortRangeVarianceData(Warning):
    def __init__(self, matrixType, obj=None):
        Warning.__init__(self, obj)
        self.matrixType = matrixType

    def __str__(self):
        return ("shortRangeSelfScalingVariance should contain a diagonal array, contains %s" % self.matrixType)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.matrixType == other.matrixType)

class ParameterCovarianceMismatch(Warning):
    def __init__(self, nParams, matrixShape, obj=None):
        Warning.__init__(self, obj)
        self.nParams = nParams
        self.matrixShape = matrixShape

    def __str__(self):
        return ("Number of linked parameters (%d) does not match covariance matrix dimension (%s)" % (self.nParams, self.matrixShape))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.nParams == other.nParams and self.matrixShape == other.matrixShape)
