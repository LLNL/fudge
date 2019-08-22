Cross sections
==============

Types of cross section data
---------------------------

Cross section data can be stored in several forms in GND. The simplest forms are *pointwise* and *linear*,
which contain a single region that uses the same interpolation rule throughout.
A 'pointwise' cross section may use any supported interpolation rule, whereas a 'linear' cross section
must use linear interpolation on both axes.

Discontinuities are not permitted inside a pointwise region. Instead, they must be implemented at the boundary
between two different regions of a *piecewise* cross section. In this case, each region has its own interpolation and list of X-Y points.

For low-energy incident neutrons, the cross section can be represented best using resonance parameters
(see the next section for more detail). In this case, the full cross section includes the sum of the values
computed from resonance with the tabulated 'background' cross section. This type of data is stored in GND as *resonancesWithBackground*.

Plotting
--------

The plotting tools included with ``Fudge`` require that the cross section first be converted to a single, lin-lin interpolated region.
This can be achieved using the 'toPointwise_withLinearXYs()' function. Thus, if you have a
gnd.reactionData.crossSection instance, you can plot it like this (assuming you have installed Gnuplot.py)::

    >>> crossSection.toPointwise_withLinearXYs().plot()

When using the toPointwiseLinear function, you may encounter the following error::

    >>> crossSection.toPointwise_withLinearXYs()
    >>> Exception: resonancesWithBackground cross section has not been reconstructed via reactionSuite.reconstructResonances

See the next section for details on how to reconstruct resonances.

Integral checks
---------------

Not yet implemented.
