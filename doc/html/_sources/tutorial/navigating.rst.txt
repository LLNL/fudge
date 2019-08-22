Navigating
==========

How to find data within a ``reactionSuite``
-------------------------------------------

The Fudge reactionSuite and covarianceSuite classes are meant to reflect the GNDS data hierarchy.
Finding data within them is often just a matter of using the same name as appears in GNDS.  For example,

    >>> from fudge.gnds import reactionSuite
    >>> RS = reactionSuite.readXML("path_to_GNDS_file.xml")
    >>> elastic = RS.getReaction("elastic")
    >>> elastic
    <fudge.gnds.reactions.reaction.reaction object at 0x1020dc590>
    >>> elastic.crossSection
    <fudge.gnds.reactionData.crossSection.component object at 0x10220dcd0>
    >>> elastic.outputChannel.Q
    <fudge.gnds.channelData.Q.component object at 0x10212a0d0>
    >>> elastic.outputChannel.products
    <fudge.gnds.product.products object at 0x102608e10>
    >>> elastic.outputChannel.products[0].multiplicity
    <fudge.gnds.productData.multiplicity.component object at 0x102610310>
    >>> elastic.outputChannel.products[0].distribution
    <fudge.gnds.productData.distributions.distribution.component object at 0x1026103d0>

The 'getReaction' method of reactionSuite provides a simple way to search for a desired reaction.
In addition to strings like "elastic", "fission" or "n + Fe56_e3" (to search for a reaction that emits a neutron
plus an Fe56 atom with the nucleus in its third excited state),  getReaction also accepts ENDF MT numbers.

Within a reaction, data like the cross section, list of products, product multiplicity and distribution information
can all be accessed by using (hopefully) intuitive names. Note, however, that the results are not lists of actual data,
instead they are instances of Fudge classes like
Finding data within them is often just a matter of using the same name as appears in GNDS.  For example,

    >>> from fudge.gnds import reactionSuite
    >>> RS = reactionSuite.readXML("path_to_GNDS_file.xml")
    >>> elastic = RS.getReaction("elastic")
    >>> elastic
    <fudge.gnds.reactions.reaction.reaction object at 0x1020dc590>
    >>> elastic.crossSection
    <fudge.gnds.reactionData.crossSection.component object at 0x10220dcd0>
    >>> elastic.outputChannel.Q
    <fudge.gnds.channelData.Q.component object at 0x10212a0d0>
    >>> elastic.outputChannel.products
    <fudge.gnds.product.products object at 0x102608e10>
    >>> elastic.outputChannel.products[0].multiplicity
    <fudge.gnds.productData.multiplicity.component object at 0x102610310>
    >>> elastic.outputChannel.products[0].distribution
    <fudge.gnds.productData.distributions.distribution.component object at 0x1026103d0>

The 'getReaction' method of reactionSuite provides a simple way to search for a desired reaction.
In addition to strings like "elastic", "fission" or "n + Fe56_e3" (to search for a reaction that emits a neutron
plus an Fe56 atom with the nucleus in its third excited state),  getReaction also accepts ENDF MT numbers.

Within a reaction, data like the cross section, list of products, product multiplicity and distribution information
can all be accessed by using (hopefully) intuitive names. Note, however, that the results are not lists of actual data,
instead they are instances of Fudge classes like
Finding data within them is often just a matter of using the same name as appears in GNDS.  For example,

    >>> from fudge.gnds import reactionSuite
    >>> RS = reactionSuite.readXML("path_to_GNDS_file.xml")
    >>> elastic = RS.getReaction("elastic")
    >>> elastic
    <fudge.gnds.reactions.reaction.reaction object at 0x1020dc590>
    >>> elastic.crossSection
    <fudge.gnds.reactionData.crossSection.component object at 0x10220dcd0>
    >>> elastic.outputChannel.Q
    <fudge.gnds.channelData.Q.component object at 0x10212a0d0>
    >>> elastic.outputChannel.products
    <fudge.gnds.product.products object at 0x102608e10>
    >>> elastic.outputChannel.products[0].multiplicity
    <fudge.gnds.productData.multiplicity.component object at 0x102610310>
    >>> elastic.outputChannel.products[0].distribution
    <fudge.gnds.productData.distributions.distribution.component object at 0x1026103d0>

The 'getReaction' method of reactionSuite provides a simple way to search for a desired reaction.
In addition to strings like "elastic", "fission" or "n + Fe56_e3" (to search for a reaction that emits a neutron
plus an Fe56 atom with the nucleus in its third excited state),  getReaction also accepts ENDF MT numbers.

Within a reaction, data like the cross section, list of products, product multiplicity and distribution information
can all be accessed by using (hopefully) intuitive names. Note, however, that the results are not lists of actual data,
instead they are instances of Fudge classes like ``crossSection.component``, ``Q.component``, ``distribution.component``,
etc.  This reflects the basic design of GNDS: each physical quantity can potentially contain more than one form of data
including the original 'evaluated' data as well as other types that may be processed by heating, grouping, etc.
Each ``component`` class contains a list of one or more ``form`` instances containing the actual data.

The components behave like Python dictionaries:

    >>> elastic.crossSection.keys()
    ['eval', 'recon']

Each of the keys corresponds to a ``style`` of data.

The concept of styles
-------------------------

Evaluated nuclear data is sometimes stored in a form that is not easily plotted or used by application codes. For
example, cross section data may be given as a combination of resonance parameters and a piecewise 'background', but
for plotting it must be converted to pointwise, linearly-interpolable data.

GNDS permits storing both of these data types (or 'styles' of data) together in the same hierarchy, so that the user has easy access to the
original style used by the evaluator as well as to the reconstructed style that was derived from it for easy plotting.
Since multiple styles of data can appear in an evaluation, GNDS needs to be able to keep them straight. A GNDS file contains
a 'styles' section near the top, where each style is defined and given a unique label.  Those same labels are used
throughout the file to indicate what style each data form is associated with.


In the example above, the cross section component contains two keys, each corresponding to a ``style``:
the 'eval' key corresponds to the original evaluated data, and 'recon'
corresponding to the data with resonances reconstructed. More detail about different styles is available in the
``styles`` section of the ``reactionSuite``.

The reconstructed data are most useful for viewing and plotting:

    >>> print( elastic.crossSection['recon'] )
       1.00000000e-05   2.09675600e+00
       6.19409249e+00   2.17147300e+00
       1.23695253e+01   2.25184300e+00
       1.89503588e+01   2.34447500e+00
        ...
       4.50000000e+07   1.43070000e+00
       5.00000000e+07   1.45270000e+00
       5.50000000e+07   1.44300000e+00
       6.00000000e+07   1.40960000e+00

(By the way, the example here is the n-025_Mn_055 candidate evaluation in ENDF-VIII-beta4).
The data has two columns: incident energy in eV, and cross section in b. We can double-check the units by printing
out the axes information:

    >>> print( elastic.crossSection['recon'].axes.toXML() )
    <axes>
      <axis index="1" label="energy_in" unit="eV"/>
      <axis index="0" label="crossSection" unit="b"/></axes>

The data are stored in another Fudge class called ``crossSection.XYs1d`` (the name indicates that this is a 1-dimensional
function stored as a list of X-Y pairs). This is a useful class, with methods for integrating, changing interpolations,
performing mathematical operations on two or more XYs1d instances, etc.  For more information on the capabilities of
XYs1d see the documentation (FIXME need link to xData here).

If desired you can also extract the data into a simple Python list of X-Y pairs:

    >>> xyPairs = elastic.crossSection['recon'].copyDataToXYs()

TO-DO expand discussion to show how to extract data from multiplicities, Q-values, etc.
