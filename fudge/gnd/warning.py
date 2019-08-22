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
"""
Store and report warnings and errors in GND data, in order to discover problems in the library.
The reactionSuite.check() function returns a nested list of warning objects:

    >>> warnings = reactionSuite.check()
    >>> print warnings

May include or exclude specific classes of warning using the filter command.
filter() returns a new context instance:

    >>> warnings2 = warnings.filter( exclude=[warning.energyImbalance] )

Or, for easier searching you may wish to flatten the list (to get warnings alone without context messages):

    >>> flat = warnings.flatten()

cmattoon, Jan 2012
"""

__metaclass__ = type

class context:
    """ 
    Store warnings in context. This class contains location information (reactionSuite, reaction, etc)
    plus a nested list of warnings or other context instances
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
        """Filter warning list to only include (or exclude) specific classes of warning. For example:
        
            >>> newWarnings = warnings.filter( exclude=[warning.energyImbalance, warning.Q_mismatch] )

        Note that if both 'include' and 'exclude' lists are provided, exclude is ignored."""

        if include is None and exclude is None: return self
        newWarningList = []
        for warning in self.warningList:
            if isinstance( warning, context ):
                newContext = warning.filter( include, exclude )
                if newContext: newWarningList.append( newContext )
            elif include is not None: 
                if warning.__class__ in include:
                    newWarningList.append( warning )
            else: # exclude is not None:
                if warning.__class__ not in exclude:
                    newWarningList.append( warning )
        return context( self.message, newWarningList )

    def flatten(self):
        """From a nested hierarchy of warnings, get back a flat list for easier searching:
        
            >>> w = reactionSuite.check()
            >>> warningList = w.flatten()"""

        List = []
        for val in self.warningList:
            if isinstance(val, warning):
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

class warning:
    """
    General warning class. Contains link to problem object,
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
# specific warning classes:
#

class NotImplemented( warning ):
    def __init__(self, form, obj=None):
        warning.__init__(self, obj)
        self.form = form

    def __str__(self):
        return "Checking not yet implemented for %s type data" % self.form

    def __eq__(self, other):
        return (self.form == other.form and self.xpath == other.xpath)

# resonance region:

class RRmultipleRegions( warning ):
    def __init__(self, obj=None):
        warning.__init__(self, obj)

    def __str__(self):
        return "Use of more than one resolved resonance regions is deprecated"

class URRmissingEnergyList( warning ):
    def __init__(self, Lval, Jval, obj=None):
        warning.__init__(self, obj)
        self.Lval = Lval
        self.Jval = Jval

    def __str__(self):
        return "Missing energy list for unresolved region L=%i / J=%.1f" % (self.Lval, self.Jval)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval)

class URRdomainMismatch( warning ):
    def __init__(self, Lval, Jval, obj=None):
        warning.__init__(self, obj)
        self.Lval = Lval
        self.Jval = Jval

    def __str__(self):
        return "Unresolved L=%i / J=%.1f widths don't span URR energy limits" % (self.Lval, self.Jval)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval)

class URRinsufficientEnergyGrid( warning ):
    def __init__(self, Lval, Jval, eLow, eHigh, obj=None):
        warning.__init__(self, obj)
        self.Lval = Lval
        self.Jval = Jval
        self.eLow = eLow
        self.eHigh = eHigh

    def __str__(self):
        return "More points needed in L=%i J=%.1f unresolved widths between %s and %s" % (self.Lval, self.Jval, self.eLow, self.eHigh)

    def __eq__(self, other):
        return (self.Lval == other.Lval and self.Jval == other.Jval and self.eLow == other.eLow and self.eHigh == other.eHigh)

# particle list:

class discreteLevelsOutOfOrder( warning ):
    def __init__(self, lidx, obj=None):
        warning.__init__(self, obj)
        self.lidx = lidx

    def __str__(self):
        return "Discrete level %i is out of order" % self.lidx

    def __eq__(self, other):
        return (self.lidx == other.lidx)

class unnormalizedGammas( warning ):
    def __init__(self, branchingSum, obj=None):
        warning.__init__(self, obj)
        self.branchingSum = branchingSum

    def __str__(self):
        return "Gamma branching = %s, should be 1.0!" % (self.branchingSum)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.branchingSum == other.branchingSum)

# reaction objects:

class WicksLimitError( warning ):
    def __init__(self, percentErr, energy_in, obj=None):
        warning.__init__(self, obj)
        self.percentErr = percentErr
        self.energy_in = energy_in

    def __str__(self):
        return "Wick's limit too low by %.3f%% at %s" % (100*self.percentErr, self.energy_in)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.percentErr == other.percentErr and self.energy_in == other.energy_in)

class ZAbalanceWarning( warning ):
    def __str__(self):
        return "ZA doesn't balance for this reaction!"

class Q_mismatch( warning ):
    def __init__(self, Qcalc, Qactual, obj=None):
        warning.__init__(self, obj)
        self.Qcalc = Qcalc
        self.Qactual = Qactual

    def __str__(self):
        return "Calculated and tabulated Q-values disagree: %s vs %s!" % (self.Qcalc, self.Qactual)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.Qcalc == other.Qcalc and self.Qactual == other.Qactual)

class threshold_mismatch( warning ):
    def __init__(self, threshold, thresholdCalc, obj=None):
        warning.__init__(self, obj)
        self.threshold = threshold
        self.thresholdCalc = thresholdCalc

    def __str__(self):
        return "Calculated and tabulated thresholds disagree: %s vs %s!" % (self.thresholdCalc, self.threshold)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.threshold == other.threshold
                and self.thresholdCalc == other.thresholdCalc)

class nonZero_crossSection_at_threshold( warning ):
    def __init__(self, value, obj=None):
        warning.__init__(self, obj)
        self.value = value

    def __str__(self):
        return "First cross section point for threshold reaction should be 0, not %s" % self.value

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.value == other.value)

class gapInCrossSection( warning ):
    def __init__(self, minGapEnergy, maxGapEnergy, obj=None):
        warning.__init__(self, obj)
        self.minGapEnergy = minGapEnergy
        self.maxGapEnergy = maxGapEnergy

    def __str__(self):
        return "Gap in cross section data from %s to %s" % (self.minGapEnergy, self.maxGapEnergy)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.minGapEnergy == other.minGapEnergy
                and self.maxGapEnergy == other.maxGapEnergy)

class negativeCrossSection( warning ):
    def __init__(self, energy_in, index, obj=None):
        warning.__init__(self, obj)
        self.energy_in = energy_in
        self.index = index

    def __str__(self):
        return "Negative cross section encountered at %s (index %i)!" % (self.energy_in, self.index)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index)

class badCrossSectionReference( warning ):
    def __str__(self):
        return "Cross section reference outside of current reactionSuite not allowed!"

class summedCrossSectionMismatch( warning ):
    def __init__(self, maxPercentDiff, obj=None):
        warning.__init__(self, obj)
        self.maxPercentDiff = maxPercentDiff

    def __str__(self):
        msg = ("Cross section does not match sum of linked reaction cross sections! Max diff: %.2f%%" %
                self.maxPercentDiff)
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.maxDiff == other.maxDiff)

class summedCrossSectionDomainMismatch( warning ):
    def __str__(self):
        return "Cross section domain does not match domain of linked reaction cross section sum"

class negativeMultiplicity( warning ):
    def __init__(self, value, obj=None):
        warning.__init__(self, obj)
        self.value = value

    def __str__(self):
        return "Encountered negative multiplicity (%s)!" % self.value

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.value == other.value)

class domain_mismatch( warning ):
    def __init__(self, lowBound, highBound, xscLowBound, xscHighBound, obj=None):
        warning.__init__(self, obj)
        self.lowBound, self.highBound = lowBound, highBound
        self.xscLowBound, self.xscHighBound = xscLowBound, xscHighBound

    def __str__(self):
        return "Domain doesn't match the cross section domain: (%s -> %s) vs (%s -> %s)" % (self.lowBound,
                self.highBound, self.xscLowBound, self.xscHighBound)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.lowBound == other.lowBound
                and self.highBound == other.highBound)

class badFissionEnergyRelease( warning ):
    def __init__(self, key, energy, val, obj=None):
        warning.__init__(self, obj)
        self.key = key
        self.energy = energy
        self.val = val

    def __str__(self):
        return "Unphysical %s fission energy release: %s at incident energy %s" % (self.key, self.val, self.energy)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.key == other.key and self.val == other.val
                and self.energy == other.energy)

class missingDistribution( warning ):
    def __init__(self, productName, obj=None):
        warning.__init__(self, obj)
        self.productName = productName

    def __str__(self):
        return "Missing distribution (required for all '%s' products)!" % self.productName

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.productName == other.productName)

""" eliminate this one?
class noDistributions( warning ):
    # more serious than missingDistributions:
    def __str__(self):
        return "No distribution data available for *any* of this reaction's products!"
"""

class missingRecoilDistribution( warning ):
    def __str__(self):
        return "Recoil distribution type specified, but recoil partner has no distribution data!"

class wrongDistributionComponent( warning ):
    def __init__(self, component, reactionType='2-body', obj=None):
        warning.__init__(self, obj)
        self.component = component
        self.reactionType = reactionType

    def __str__(self):
        return "%s is not a valid distribution component for %s reaction" % (self.component, self.reactionType)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.component == other.component
                and self.reactionType == other.reactionType)

class wrong2BodyFrame( warning ):
    def __str__(self):
        return ("2-body reaction not in center-of-mass frame!")

class uncorrelatedFramesMismatch( warning ):
    def __init__(self, angleFrame, energyFrame, obj=None):
        warning.__init__(self, obj)
        self.angleFrame = angleFrame
        self.energyFrame = energyFrame

    def __str__(self):
        return ("For uncorrelated energy/angle distributions, frame must be identical! Currently, angle given in %s, but energy in %s"
                % (self.angleFrame, self.energyFrame))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.angleFrame == other.angleFrame
                and self.energyFrame == other.energyFrame)

class energyDistributionBadU( warning ):
    def __str__(self):
        return ("For energy distribution functional form, parameter 'U' results in negative outgoing energies!")

class MadlandNixBadParameters( warning ):
    def __init__(self, EFL, EFH, minTm, obj=None):
        warning.__init__(self, obj)
        self.EFL = EFL
        self.EFH = EFH
        self.minTm = minTm

    def __str__(self):
        return ( "Madland-Nix fission spectrum contains negative parameters: EFL=%s, EFH=%s, min(Tm)=%s"
                % (self.EFL, self.EFH, self.minTm) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.EFL == other.EFL
                and self.EFH == other.EFH and self.minTm == other.minTm)

class weightsDontSumToOne( warning ):
    def __str__(self):
        return ("Weights don't sum to 1.0!")

class unnormalizedDistribution( warning ):
    def __init__(self, energy_in, index, integral, obj=None):
        warning.__init__(self, obj)
        self.energy_in = energy_in
        self.index = index
        self.integral = integral

    def __str__(self):
        return "Unnormalized distribution! At energy_in = %s (index %i), integral = %s" % (self.energy_in,
                self.index, self.integral)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.index == other.index and self.integral == other.integral)

class unnormalizedKMDistribution( unnormalizedDistribution ):
    pass    # identical to unnormalizedDistribution, except it needs a special fix() function

class incompleteDistribution( warning ):
    def __init__(self, energy_in, lowerMu, upperMu, obj=None):
        warning.__init__(self, obj)
        self.energy_in = energy_in
        self.lowerMu = lowerMu
        self.upperMu = upperMu

    def __str__(self):
        return ("Incomplete angular coverage (mu = %.2f to %.2f) for distribution at incident energy %s"
                % (self.lowerMu, self.upperMu, self.energy_in))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.lowerMu == other.lowerMu and self.upperMu == other.upperMu)

class negativeProbability( warning ):
    def __init__(self, energy_in, energy_out=None, mu=None, value=None, obj=None):
        warning.__init__(self, obj)
        self.energy_in = energy_in
        self.energy_out = energy_out
        self.mu = mu
        self.value = value

    def __str__(self):
        msg = "Negative probabilities encountered in distribution. Incident energy: %s" % self.energy_in
        if self.mu is not None: msg += ", mu: %s" % self.mu
        if self.energy_out is not None: msg += ", outgoing energy: %s" % self.energy_out
        if self.value is not None: msg += ", worst case: %s" % self.value
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in
                and self.energy_out == other.energy_out and self.mu == other.mu)

class extraOutgoingEnergy( warning ):
    """If an outgoing energy distribution ends with more than one energy with probability=0,
    proper unitbase treatment is unclear. Distribution should end with exactly one P=0 point."""

    def __init__(self, energy_in, obj=None):
        warning.__init__(self, obj)
        self.energy_in = energy_in

    def __str__(self):
        msg = "Extra zero-probability outgoing energies found at incident energy %s" % self.energy_in
        return msg

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.energy_in == other.energy_in)

class energyImbalance( warning ):
    def __init__(self, energy_in, index, availableEnergy, deposition_per_product, obj=None):
        warning.__init__(self, obj)
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

class fissionEnergyImbalance( energyImbalance ):
    def __str__(self):
        per_product = ', '.join(["%s = %.4g%%" % (key,val) for key,val in self.deposition_per_product[:5]])
        if len( self.deposition_per_product ) > 5: per_product += ', ...'
        return ("Fission energy imbalance at incident energy %s (index %i). Total deposited = %.4g%% (%s), leaving insufficient energy for fission products!" %
                (self.energy_in, self.index, self.total_deposited, per_product))

class valueOutOfRange( warning ):
    def __init__(self, contextMessage, value, lowerBound, upperBound, obj=None):
        warning.__init__(self, obj)
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

class ExceptionRaised( warning ):
    """ if we run into an exception when running check(), try to exit gracefully and return this warning """
    def __init__(self, Exception_String, obj=None):
        warning.__init__(self, obj)
        self.Exception_String = Exception_String

    def __str__(self):
        return ("Encountered Exception: %s" % self.Exception_String)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.Exception_String == other.Exception_String)

class EnergyDepositionExceptionRaised( ExceptionRaised ):
    def __str__(self):
        return ("Exception raised when calculating energy deposition: %s" % self.Exception_String)

### covarianceSuite warnings: ###

class cyclicDependency( warning ):
    def __init__(self, cycle, obj=None):
        warning.__init__(self, obj)
        self.cycle = cycle

    def __str__(self):
        return ("Cyclic dependency in summed covariances for sections %s and %s" % tuple(self.cycle))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.cycle == other.cycle)

class varianceTooSmall( warning ):
    def __init__(self, index, variance, obj=None):
        warning.__init__(self, obj)
        self.index = index
        self.variance = variance

    def __str__(self):
        return ("Encountered very small variance (%e%%) at index %i." % (self.variance, self.index))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.variance == other.variance)

class varianceTooLarge( warning ):
    def __init__(self, index, variance, obj=None):
        warning.__init__(self, obj)
        self.index = index
        self.variance = variance

    def __str__(self):
        return ("Encountered very large variance (%e%%) at index %i." % (self.variance, self.index))

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.index == other.index
                and self.variance == other.variance)

class negativeEigenvalues( warning ):
    def __init__(self, negativeCount, worstCase, obj=None):
        warning.__init__(self, obj)
        self.negativeCount = negativeCount
        self.worstCase = worstCase

    def __str__(self):
        return ("%i negative eigenvalues! Worst case = %e" % (self.negativeCount, self.worstCase) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.negativeCount == other.negativeCount
                and self.worstCase == other.worstCase)

class badEigenvalueRatio( warning ):
    def __init__(self, ratio, obj=None):
        warning.__init__(self, obj)
        self.ratio = ratio

    def __str__(self):
        return ("Ratio of smallest/largest eigenvalue (%e) is too small" % (self.ratio) )

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.ratio == other.ratio)
