isotopicAbundances
==================

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

It is worth noting that the format for storing isotopic abundance data described in this sub-package
are not part of the **GNDS** specifications. Support for isotopic abundance data was added to **FUDGE**
at the request of **LLNL** users. That said, the format does uses stuff from the **GNDS** specifications.

Various isotopic abundance database exists. For example, both NIST (National Institute of Standards and Technology)
and CIAAW (Commission on Isotopic Abundances and Atomic Weights) have an isotopic abundance database.
To distinguish different database, each database has a **evaluation** attribute.

The top node of an isotopic abundance database is the **isotopicAbundancesByChemicalElement** node. This node
has 4 parts which are:

    * evaluation:       a unique name to distinguish one database from others.
    * format:           specifies the format the data are stored in. Currently, the only defined format is "2.0".
    * documentation:    a **GNDS** documentation node.
    * chemicalElements: a list of **chemicalElement** nodes, one for each chemical element in the database.

The **chemicalElement** node has 3 parts which are:

    * symbol:           the **GNDS** **PoPs** symbol for the chemical element.
    * documentation:    a **GNDS** documentation node.
    * isotopes:         a list of **isotope** nodes, one for each isotope in the chemical element with abundance data.

The **isotope** has 3 parts which are:

    * id:               the **GNDS** **PoPs** id for the isotope.
    * atomFraction:     the atom fraction of the isotope.
    * uncertainty:      the uncertainty of the atom fraction. If uncertainly is unknown, its value will be -1.

The following is an example of what a evaluation looks like in XML::

    <isotopicAbundancesByChemicalElement evaluation="NIST" format="2.0">
      <documentation/>
      <chemicalElements>
        <chemicalElement symbol="H">
          <isotopes>
            <isotope id="H1" atomFraction="0.999885" uncertainty="7e-5"/>
            <isotope id="H2" atomFraction="1.15e-4" uncertainty="7e-5"/></isotopes></chemicalElement>
        <chemicalElement symbol="He">
          <isotopes>
            <isotope id="He3" atomFraction="1.34e-6" uncertainty="3e-8"/>
            <isotope id="He4" atomFraction="0.99999866" uncertainty="3e-8"/></isotopes></chemicalElement>
        <chemicalElement symbol="Li">
            ... <chemicalElement></chemicalElements></isotopicAbundancesByChemicalElement>

To support multiple databases, the **iam** (Isotopic Abundance Map) node was created which stores
information about the available isotopic abundance databases. The **iam** has 5 parts which are:

    * library:              do we need this.
    * format:               specifies the format the data are stored in. Currently, the only defined format is "2.0".
    * checksum:             a checksum calculated by concatenating the checksum of all the entry checksums.
    * algorithm:            the algorithm used to calculate *checksum*.
    * a list of entries:    a list of mixed **iam** **isotopicAbundancesByChemicalElement** and **import** nodes.

The **iam** **isotopicAbundancesByChemicalElement** has has 4 parts which are:

    * path:                 the path to a evaluated isotopic abundance file.
    * evaluation:           name for the evaluation that **path** refers to.
    * checksum:             a checksum calculated by concatenating the checksum of all the entry checksums.
    * algorithm:            the algorithm used to calculate *checksum*.

The **iam** **import** node allow for nesting of **iam** nodes and has has 3 parts which are:

    * path:                 the path to the **iam** file to import
    * checksum:             a checksum calculated by concatenating the checksum of all the entry checksums.
    * algorithm:            the algorithm used to calculate *checksum*.

The following Python code shows how to read in an evaluated isotopic abundance file "NIST.iam",
get the isotopic abundance data for Oxygen:

    >>> from isotopicAbundances import isotopicAbundances as isotopicAbundancesModule
    >>> NIST = isotopicAbundancesModule.read("NIST.iam")
    >>> NIST.evaluation
    'NIST'
    >>> Oxygen = NIST['O']
    >>> for isotope in Oxygen:
    ...     print(isotope.id, isotope.atomFraction, isotope.uncertainty)
    O16 0.99757 0.00016
    O17 0.00038 1e-05
    O18 0.00205 0.00014

The following shows how to get the same data via an **iam** file which has access to four evaluations:

    >>> from isotopicAbundances import iam as iamModule
    >>> iam = iamModule.read('all.iam')iam = iamModule.read('all.iam')
    >>> iam.evaluations()                       # Prints a list of available evaluatioins.
    ['NIST', 'FUDGE', 'COGZA2a', 'CIAAW']
    NIST = iam.findEvaluation('NIST')           # This is the same as the NIST database read in above.
    >>> NIST.evaluation
    'NIST'

.. toctree::
   :maxdepth: 4

isotopicAbundances.iam module
-----------------------------

.. automodule:: isotopicAbundances.iam
    :members:
    :undoc-members:
    :show-inheritance:

isotopicAbundances.isotopicAbundances module
--------------------------------------------

.. automodule:: isotopicAbundances.isotopicAbundances
    :members:
    :undoc-members:
    :show-inheritance:
