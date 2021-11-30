Outgoing Distributions
======================

Applications of nuclear data are often sensitive not only to the reaction cross section, but also to the angular distributions and energy 
spectra of the reaction products. Each product in GNDS therefore contains a 'distribution' element. The distribution data may be given in 
several different ways.

Angular distributions
---------------------

For reactions producing exactly two outgoing products (such as n + Mn55 -> p + Cr55), the evaluation only needs to include an angular 
distribution for one product. The kinematics of these 'two-body' reactions mean that the angular distribution for the other product 
(and the energy spectra for both) can be computed from the first product.

There are several ways of describing angular distributions in GNDS:

    - Isotropic
    - Legendre coefficients for each incident energy,
    - Pointwise data for each incident energy, where P(mu) is given explicitly for mu from -1 to 1 (mu = cosine of the scattering angle theta)
    - 'Mixed' form, where Legendre coefficients are used at smaller incident energies and pointwise at higher energies (Legendre coefficients 
      can cause trouble with negative probabilities, especially for very forward-peaked distributions).

Energy distributions
--------------------

For a two-body reaction, an energy distribution could be supplied instead of an angular distribution. In practice, however, energy distributions 
are only given for non-2-body reactions. In GNDS, these distributions are used when the outgoing angle and energy are assumed to be uncorrelated. 
There are many ways of storing energy spectra, including:

    - Pointwise
    - Watt spectra
    - Maxwellian
    - and many others

Joint energy-angle distributions
--------------------------------

For non-two-body reactions, the distribution can be best described using a double-differential form that gives the probability as a functon of both outgoing angle *and* energy. Some types of double-differential distribution forms include:

- Legendre coefficients listed as a function of both incident and outgoing energy,
- Pointwise data, where P(mu | E, E' ) is given for each combination of incident energy E and outgoing energy E',
- Kalbach-Mann data, a parameterized form that can efficiently describe the average outgoing spectra.

Multiplicity
------------

Reaction products require a multiplicity. For many common reactions the multiplicity is a constant integer (i.e. n,2n,
n,n', elastic, etc.).  In some special cases (mainly gammas and fission neutrons), an average,
incident-energy-dependent multiplicity is listed instead. Methods for storing the multiplicity include

- constant1d
- XYs1d
- regions1d
- polynomial1d
- reference (i.e. a link to another product's multiplicity)

Where are the discrete gammas?
------------------------------

Discrete gammas can appear in several places in a GNDS evaluation. Where possible, they are located inside
the PoPs particle database, where they are listed as decay modes from excited nuclear states down to lower energy states.
In other cases, the source of the gamma may not be known exactly (or may not have been included in the ENDF evaluation
that was translated to GNDS).  In those cases, discrete gammas are listed as normal reaction products,
with uncorrelated distributions here the energy distribution is listed as ``discreteGamma``.
Primary gammas are listed in a similar fashion, as uncorrelated distributions with a ``primaryGamma`` energy distribution.
