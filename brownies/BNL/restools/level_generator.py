import math
import numpy

from fudge.core.math.pdf import UnivariatePDF, WignerDistribution, PoissonDistribution, \
    BrodyDistribution, GOEDistribution
from xData import XYs1d as XYs1d

"""
Collection of fake level sequence generators, including the One True Generator: getGOEFakeLevelSequence
"""

fakeLevelStyles = ['wigner', 'picket fence', 'poisson', 'brody', 'goe']


def getFakeLevelSequence(E0=0.0, aveD=None, numLevels=None, style='goe', BrodyW=0.5, levelDensity=None):
    """
    wrapper function

    :param E0: see documentation for individual styles; note: for GOE is best to use E0=0.0
    :param aveD: see documentation for individual styles
    :param numLevels: see documentation for individual styles
    :param style: one of ['wigner','picket fence', 'poisson', 'brody', 'goe']
    :param BrodyW: see documentation for individual styles
    :param levelDensity: see documentation for individual styles
    :return:
    """

    # To generate non-GOE levels, we need to know the average level spacing, the starting energy and the number
    # of levels to build.  We can get this information a number of different ways.  These are the options:
    #   * Easiest way: aveD, E0 and numLevels
    #   * aveD, E0 and upperBound, compute numLevels = 1+(upperBound - E0)/D
    #   * levelDensity and E0.  Internally getFakeLevelSequence converts this to aveD and numLevels.
    if style != 'goe':
        if aveD is None:
            if levelDensity is not None:
                aveD = 1 / levelDensity.evaluate(0.5 * (levelDensity.domainMin + levelDensity.domainMax))
            else:
                raise ValueError("Not enough information to determine aveD")
        if numLevels is None:
            if levelDensity is not None:
                numLevels = 1+(levelDensity.domainMax - E0)/aveD
            else:
                raise ValueError("Not enough information to determine numLevels")

    # To generate GOE levels, we need basically the same information, but the GOE routine uses the levelDensity
    # instead of the aveD for a more finely tuned reproduction of the fake level scheme.
    if style == 'goe':
        if levelDensity is None:
            raise ValueError("For GOE style, need a level density")
        if numLevels is None:
            numLevels = numpy.random.poisson(levelDensity.integrate())
            print("setting numLevels to", numLevels)

    if style == 'wigner':
        return getWignerFakeLevelSequence(E0, 1/levelDensity)
    elif style == 'picket fence':
        return getPicketFenceFakeLevelSequence(E0, aveD, numLevels)
    elif style == 'poisson':
        return getPoissonFakeLevelSequence(E0, aveD, numLevels)
    elif style == 'brody':
        return getBrodyFakeLevelSequence(E0, aveD, numLevels, BrodyW)
    elif style == 'goe':
        return getGOEFakeLevelSequence(E0, numLevels, levelDensity)
    else:
        raise ValueError("style must be one of " + str(fakeLevelStyles))


def sample_goe_matrix(ndim, scale=1.0):
    """
    Create a ndim x ndim GOE "Hamiltonian" matrix following Declan Mulhall's algorithm.

    :param ndim: an int, the dimension of the matrix
    :param scale: an energy scale factor.  This sets the variance of the Gaussian used to generate the GOE matrix
    :return: a numpy ndim x ndim GOE
    """
    goe = numpy.random.normal(loc=0.0, scale=scale, size=(ndim, ndim))
    goe = (goe + goe.T) / math.sqrt(2.0)
    for i in range(ndim):
        goe[i, i] *= math.sqrt(2.0)
    return goe


def sample_goe_eigenvalues(ndim, normalize=True):
    """
    Generate a GOE spectrum by making a GOE "Hamiltonian" and then diagonalizing it

    Note: on Dave's MacBook Pro, we can't handle too many levels (<500 safer).  Don't know why

    :param ndim: number of eigenvalues (equivalently the size of the GOE "Hamiltonian")
    :param normalize: Ensure the eigenmodes are on the interval [-1,1] instead of [-ndim/2,ndim/2]
    :return:list of eigenvalues
    """
    if normalize:
        scale = 0.5 / math.sqrt(ndim)
    else:
        scale = 1.0
    sample = numpy.linalg.eigvals(sample_goe_matrix(ndim, scale=scale))
    sample.sort()
    return sample


def getWignerFakeLevelSequence(E0, levelSpacing):
    """
    Random Matrix Theory (RMT) predicts that the Nearest Neighbor Spacing Distribution (NNSD) will have the
    shape of a Wigner distribution.  If you make levels by drawing spacings from a Wigner distribution,
    by construction you have the correct NNSD.

    :param E0: first level of the sequence
    :param levelSpacing: energy-dependent level spacing, assumed to be in same units as E0
    :return: the list of level energies
    """
    result = [E0]
    WD = WignerDistribution()
    domainMax = levelSpacing.domainMax
    while True:
        s = WD.drawSample()
        result.append(result[-1] + s * levelSpacing.evaluate(result[-1]))
        if result[-1] > domainMax: break
    result.pop()    # last point is beyond the levelSpacing domain
    return result


def getModifiedWignerFakeLevelSequence(E0, levelSpacing):
    """
    Because creating levels using the getWignerFakeLevelSequence above gets the right NNSD,
    but fails to create the longer range spectral correlations (e.g. spectral stiffness),
    I kludged together a scheme that builds some stiffness into the generated sequence.

    :param E0: first level of the sequence
    :param levelSpacing: energy-dependent level spacing, assumed to be in same units as E0
    :return: the list of level energies
    """
    result = [E0]
    WD = WignerDistribution()
    domainMax = levelSpacing.domainMax
    while True:
        s = WD.drawSample()
        result.append(result[-1] + s * levelSpacing.evaluate(result[-1]))
        if result[-1] > domainMax: break
        result.append(result[-1] + (1. - s) * levelSpacing.evaluate(result[-1]))
        if result[-1] > domainMax: break
    result.pop()    # last point is beyond the levelSpacing domain
    return result


def getPicketFenceFakeLevelSequence(E0, aveD, numLevels):
    """
    An evenly spaced set of fake resonances, separated by energy aveD.  This gets the
    level repulsion right, but otherwise it is so so wrong.

    :param E0: first level of the sequence
    :param aveD: average level spacing, assumed to be in same units as E0
    :param numLevels: number of levels to manufacture
    :return: the list of level energies
    """
    return [E0 + s * aveD for s in range(numLevels)]


def getPoissonFakeLevelSequence(E0, aveD, numLevels):
    """
    A Poisson distribution keeps the energies positive, but ignores the level repulsion built into sampling
    from a Wigner distribution.

    :param E0: first level of the sequence
    :param aveD: average level spacing, assumed to be in same units as E0
    :param numLevels: number of levels to manufacture
    :return: the list of level energies
    """
    result = [E0]
    for s in PoissonDistribution().drawSample(numLevels - 1):
        result.append(result[-1] + s * aveD)
    return result


def getBrodyFakeLevelSequence(E0, aveD, numLevels, BrodyW=0.5):
    """
    Brody's scheme to interpolate between Wigner and Poisson distributions (need reference, I lost it)

    :param E0: first level of the sequence
    :param aveD: average level spacing, assumed to be in same units as E0
    :param BrodyW: trade-off parameter
    :param numLevels: number of levels to manufacture
    :return: the list of level energies
    """
    result = [E0]
    for s in BrodyDistribution(w=BrodyW).drawSample(numLevels - 1):
        result.append(result[-1] + s * aveD)
    return result


def getCLDInverse(levelDensity):
    if not isinstance(levelDensity, XYs1d.XYs1d):
        raise TypeError("For GOE, levelDensity must be a XYs instance")

    DOPLOTS = False

    # Normalized level density as a PDF
    totalNumLevels = int(levelDensity.integrate()) #.getValue())
    fakePDF = levelDensity / totalNumLevels
    fakePDF.axes[0].label = "PDF(E)"
    if DOPLOTS:
        fakePDF.plot(title='fake PDF')

    # Convert it to a CDF
    fakeCDF = fakePDF.indefiniteIntegral()
    fakeCDF.axes[0].unit = ""
    fakeCDF.axes[0].label = "CDF(E)"
    if DOPLOTS:
        fakeCDF.plot(title="fake CDF")

    # Invert it to make a probability -> energy converter
    return fakeCDF.inverse()


def getGOEFakeLevelSequence(E0, totalNumLevels, levelDensity, paddingNumLevels=100, keepLevelsAboveDomainMax=False,
                            DOPLOTS=False):
    """
    This generates a sequence of energies that is both consistent with the GOE of RMT and has the correct
    secular variation contained in the levelDensity

    :param E0: Levels are generated with E0 as an energy offset. E0=0 is recommended for GOE realizations
    :param totalNumLevels: Number of fake levels to generate
    :param levelDensity: Level density we're trying to emulate with a fake GOE inspired level scheme.
                         Should be an XYs1d
    :param paddingNumLevels: We want to grab eigenvalues from the center of the GOE spectrum to avoid edge effects.
                             This is the number of extra eigenvalues to generate to pad the ends of the GOE spectrum.
    :param keepLevelsAboveDomainMax: If True, keep any levels above the domainMax of the levelDensity.
                                     If False, the resulting level scheme may (or may not) have the same number of
                                     levels as totalNumLevels.
    :param DOPLOTS: If True, make some plots
    :return: the list of level energies
    """
    # Get the GOE single eigenvalue distribution as a PDF, this is the secular variation of
    # the eigenvalues of the full GOE
    if E0 != 0:
        print("WARNING: non-zero offset is discouraged when using GOE realizations")

    a = 1.0

    # Sample a GOE set of "energies"
    dimension = totalNumLevels + 2*paddingNumLevels
    goeSample = numpy.linalg.eigvalsh( sample_goe_matrix( dimension, scale=a/2.0/math.sqrt(dimension) ) )

    # Because we only have a finite number of levels, the GOE distribution never completely matches the
    # Wigner semi-circle law (encoded in the GOEDistribution class).  The biggest deviations from the semi-circle
    # law happen on the fringes.  To combat this, we discard paddingNumLevels from either end of the
    # simulated spectrum.  We also have to discard the same region of the semi-circle distribution.
    if paddingNumLevels > 0:
        goeSample = goeSample[paddingNumLevels:-paddingNumLevels]
        goeSingleLevels = GOEDistribution(a, domainMin=goeSample[0], domainMax=goeSample[-1])
        goeSingleLevels = goeSingleLevels.normalize()
        goeSingleLevels = UnivariatePDF(XYs1d=goeSingleLevels)
    else:
        goeSingleLevels = GOEDistribution(a, domainMin=-a, domainMax=a)

    if DOPLOTS:
        goeSingleLevels.plot()
        goeSingleLevels.cdf.plot()

    # The list of x's have the full correlations of the GOE in them, except the gross secular variation.
    # We'll add that back using the real level density.  The "refolds" back in the correct secular variation.
    result = unfoldThenRefoldLevelSequence(
        originalSequence=goeSample,
        originalSequenceSecularVariationDistribution=goeSingleLevels,
        finalSecularVariationFunction=levelDensity,
        offset=E0,
        DOPLOTS=DOPLOTS)

    # This is Monte Carlo, so we can accidentally sample things that are too high.  Let's get rid of them
    if not keepLevelsAboveDomainMax:
        result = [x for x in result if x <= levelDensity.domainMax]

    return result


def unfoldThenRefoldLevelSequence(
        originalSequence,
        originalSequenceSecularVariationDistribution,
        finalSecularVariationFunction,
        offset=0.0,
        DOPLOTS=False):
    """
    Unfold an original energy sequences secular variation, then add back a different secular variation
    This is useful for taking say a GOE level sequence and then stretching it to match a known experimental one.

    The final level sequence then has the large energy scale secular variation of the known experimental one with the
    short range statistical variation of the original sequence.

    :param originalSequence: Original sequence of "energies" whose secular variance is given by
                             originalSequenceSecularVariationDistribution
    :param originalSequenceSecularVariationDistribution: a distribution inheriting from fudge.core.pdf.UnivariatePDF
    :param finalSecularVariationFunction: A level density like object that encodes the final variation.
    :param offset: add an energy offset to the final list of energies
    :param DOPLOTS: Flag to trigger plotting of intermediate distributions (useful for debugging)
    :return:
    """

    # Check types of inputs
    if not isinstance(originalSequenceSecularVariationDistribution, UnivariatePDF):
        raise TypeError("Original sequence's secular variation must be given by an instance of UnivariatePDF")
    if DOPLOTS:
        originalSequenceSecularVariationDistribution.plot()
        originalSequenceSecularVariationDistribution.cdf.plot()

    # Get the original "energies" and remove the secular variation; this is the "unfolding" step
    xList = [originalSequenceSecularVariationDistribution.cdf.evaluate(x) for x in originalSequence]
    #xList = list(map(originalSequenceSecularVariationDistribution.cdf.evaluate, originalSequence))

    # Need to invert the cumulative level distribution for refolding
    invFakeCDF = getCLDInverse(finalSecularVariationFunction)
    if DOPLOTS:
        invFakeCDF.plot()  # title='lookup table')

    # We'll add that back using the real level density.  The "refolds" back in the correct secular variation.
    result = [offset + invFakeCDF.evaluate(x) for x in xList if invFakeCDF.evaluate(x) is not None]
    result.sort()

    return result
