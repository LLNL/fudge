Resonances
==========

Resolved
--------

In the resolved resonance region, each resonance is listed along with partial widths for the various open outgoing
reaction channels. From these parameters, cross sections for each open channel can be reconstructed.

The following formalisms for storing resolved region resonance parameters are supported in GND:
        * BreitWigner (used to store equivalent of LRF 1 and 2 in the ENDF-6 format)
        * RMatrix (LRF 3 and 7)

The Adler-Adler formalism (ENDF LRF 4) is still allowed in ENDF-6 formatted files, but no evaluations using
this parameterization are present anymore in the ENDF/B-VII.1 library. Thus GND does not yet support Adler-Adler
resonance parameters.

Unresolved
----------

In the unresolved resonance region, the density of available excited nuclear states increases and resonances begin to
overlap. Here, average widths and level densities are given as a function of resonance L and J values for each energy
interval. From these parameters, an average cross section for each open channel can be reconstructed.

While the ENDF-6 format supports several ways of storing unresolved resonance parameters, fundamentally they all use
the Single-Level Breit Wigner approximation. GND always stored unresolved parameters in tables, containing partial
widths for each energy interval, for each combination of resonance L and J.

Resonance reconstruction
------------------------

Resonance reconstruction in ``Fudge`` is simple: starting with a reactionSuite instance, we use the
'reconstructResonances' method. For example:

        >>> from fudge.gnd import reactionSuite, styles
        >>> RS = reactionSuite.readXML( "gnd_formatted_file.xml" )
        >>> evalStyle = RS.styles.getEvaluatedStyle()
        >>> reconStyle = styles.crossSectionReconstructed( "reconstructed", derivedFrom=evalStyle.label )
        >>> RS.reconstructResonances( reconStyle, accuracy=0.001, verbose=True )

This method reconstructs the pointwise cross section from resonance parameters, and then adds the background cross
section to obtain the full cross section at all incident energies. The background is typically required since
resonance parameters usually only apply at small incident energies. See the 'resonancesWithBackground' section
(need link to cross section page) for more detail.

Before reconstructing resonances, we must create a new 'style' associated with the reconstructed values.
After reconstruction, cross sections for resonance reactions will have an additional data form with the label
"reconstructed".

The 'reconstructResonances' method takes two optional parameters. If an 'accuracy' is supplied, extra points are added
during reconstruction to ensure that the resulting cross section can be interpolated to the desired
accuracy. The 'verbose' option can be used to print out some details during resonance reconstruction.

'reconstructResonances' is meant as an easy-to-use interface for resonance reconstruction. Under the covers, it uses
the Python module at ``fudge.processing.resonances.fudgeReconstructResonances``. For improved performance
that module uses numpy, multiprocessing, and also some compiled extensions, written in C, that handle the most
computationally intensive tasks.


Angular distributions from resonances
-------------------------------------

Modern resonance parameterizations can be used to describe not only the cross section as a function of incident
energy, but also the outgoing angular distributions for products. For more detail, see section D.1.7.2 of the
ENDF-6 manual (available at https://ndclx4.bnl.gov/gf/project/endf6man).

Angular distribution reconstruction capability has not yet been added to ``Fudge``, although we expect to add
this feature soon.
