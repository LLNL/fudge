###############################
Sub-packages
###############################

In addition to the *fudge* directory, the **FUDGE** code package includes the following required, parallel directories that
are considered code sub-packages of **FUDGE**. In general, **FUDGE** requires the sub-packages but the sub-packages
do not require **FUDGE**. The sub-packages are included with the *FUDGE* installation.

Sub-packages
============

LUPY
------------------------------------
This sub-package contains functions and classes that are used in various other sub-packages.

pqu
------------------------------------
This sub-package contains classes and functions for representing and manipulating a value with units and uncertainty,
herein called a "Physical Quantity with Uncertainty" or PQU. For example, the following code shows some examples:

    >>> from pqu.PQU import PQU
    >>> length = PQU(10, 'm')
    >>> time = PQU(2.4, 's')
    >>> speed = length / time
    >>> speed.getValueAs('mi / h')
    9.320567883560011

numericalFunctions
------------------------------------
This sub-package contains C functions for storing and operating on a 1-d pointwise function include python wrappers.
In general, the user should not need to interact with this sub-package as it is used by the class
**XYs1d** in the **xData** sub-package and it is the **XYs1d** class that users should use.

xData
------------------------------------
This sub-package contains functions and classes for storing and operating on a 1, 2 and 3-d pointwise functions.

PoPs
------------------------------------
This sub-package contains functions and classes for working with **GNDS** **PoPs** database. A **PoPs** (Properties
of Particles) database contains a list of particles, and for particle, its mass, spin, parity, etc. can be
stored.

crossSectionAdjustForHeatedTarget
------------------------------------
This sub-package contains C functions for "heating" a nuclear cross section for a projectile hitting a material where
the effects of the thermal motion of the atoms in the material are taken into account.

Merced
------------------------------------
This sub-package contains C++ for calculating multi-group transfer matries. The user should not need to 
interact with this sub-package as it is called from within **FUDGE** via the **processMultiGroup** methods.

brownies
------------------------------------

isotopicAbundances
------------------------------------

This sub-package supports isotopic abundance data stored in an *isotopicAbundancesByChemicalElement* node.
Basically, the *isotopicAbundancesByChemicalElement* node is a list of chemical elements where each
chemical element is a list of "natural" occurring isotopes for that chemical element. Each isotope
has its "natural" occurring abundance and with uncertainty. For example,
one evaluation has the following data for the chemical element Oxygen::

    <chemicalElement symbol="O">
      <isotopes>
        <isotope id="O16" atomFraction="0.99757" uncertainty="1.6e-4"/>
        <isotope id="O17" atomFraction="3.8e-4" uncertainty="1e-5"/>
        <isotope id="O18" atomFraction="2.05e-3" uncertainty="1.4e-4"/></isotopes></chemicalElement>

Sub-package dependencies
========================
This section list the dependencies of one sub-package on the other sub-packages.

    * LUPY: none
    * pqu: none
    * numericalFunctions: none
    * xData: numericalFunctions
    * PoPs: LUPY, xData
    * crossSectionAdjustForHeatedTarget: none
    * Merced: none
    * brownies: LUPY, pqu, xData and PoPs
    * isotopicAbundances: LUPY, pqu and xData

bin directory
=============
In general, **FUDGE** should be considered as a toolkit for interaction which nuclear data.
Ergo, the user has to write Python scripts to read, modify, process and write nuclear data.
This directory contains some useful scripts that have been written by the **FUDGE** developers.
