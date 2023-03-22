Cross sections
==============

Types of cross section data
---------------------------

Cross section data can be stored in several forms in GNDS. The simplest form is *XYs1d* with *lin-lin* interpolation,
a single region using the same interpolation rule throughout.
An 'XYs1d' cross-section may use any supported interpolation rule: *lin-lin*, *lin-log*, *log-lin*, *log-log*,
*flat* and in special cases a *chargedParticle* interpolation rule.

Discontinuities are not permitted inside an XYs1d function. Instead, they must be implemented at the boundary
between two different regions of a *regions1d* container. In this case, the regions1d contains a list of 1-d functions (typically XYs1d objects),
each with its own interpolation rule.

For low-energy incident neutrons, the cross section can be represented best using resonance parameters
(see the next section for more detail). In this case, the full cross section includes the sum of the values
computed from resonance plus the tabulated 'background' cross section (which may be 0).
This type of data is stored in GNDS as a *resonancesWithBackground* cross-section form.

Plotting
--------

The plotting tools included with ``Fudge`` require that the cross section first be converted to a single, lin-lin interpolated region.
This can be achieved using the 'toPointwise_withLinearXYs()' function. Thus, if you have a
fudge.reactionData.crossSection.Form instance, you can plot it like this (assuming you have installed matplotlib)::

    >>> from matplotlib import pyplot
    >>> xs, ys = crossSection.toPointwise_withLinearXYs().copyDataToXsAndYs()
    >>> pyplot.plot(xs, ys)

When using the toPointwiseLinear function, you may encounter the following error::

    >>> crossSection.toPointwise_withLinearXYs()
    >>> Exception: resonancesWithBackground cross section has not been reconstructed via reactionSuite.reconstructResonances

See the next section for details on how to reconstruct resonances.

Integral checks
---------------

Not yet implemented.
