``Fudge`` idiosyncrasies
========================

The concept of nativeData
-------------------------

Evaluated nuclear data is sometimes stored in a form that is not easily plotted or used by application codes. For
example, cross section data may be given as a combination of resonance parameters and a piecewise 'background', but
for plotting it must be converted to pointwise, linearly-interpolable data.

GND permits storing both of these data types (or 'forms') together in the same hierarchy, so that the user has easy access to the
original (or 'native') form used by the evaluator and to the linearized form that was derived from it for easy plotting.
In order to differentiate between the two forms, GND uses a special tag called 'nativeData' to designate the original form.

Data components
---------------


