``Fudge`` idiosyncrasies
========================

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

Data components
---------------


