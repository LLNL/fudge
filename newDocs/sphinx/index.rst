######################
    FUDGE and GNDS
######################

:Release: |version|
:Date: |today|

Overview
********
The FUDGE (For Updating Data and GEnerating Evaluations) is a python based nuclear data platform that is used
to process, write, modify and view nuclear data (Beck and Mattoon, 2014). It is part of a suite of codes being 
developed at the Lawrence Livermore National Laboratory (LLNL) to support the new General Nuclear Database 
Structure (GNDS) which was developed to mordernize the storage of nuclear data. GNDS (Brown et al., 2019) 
provides representations for cross sections and distributions, particleproduction from any reaction, photo- 
and electro-atomic interaction data, thermal neutron scattering data, and radionuclide production and decay 
data (including fission products).

It was concieved as a nested hierarchy with a `dataset`, an `element` and an `attribute` as the three essential 
building blocks (Mattoon et al., 2012). The table below lists metalanguages that have been used to represent GNDS. 
Although JSON lacks explicit support for the GNDS `attribute`, a special object with the string `attribute` as key 
and an object with the GNDS `attribute` as value. 

.. todo::
   Complete table showing how FUDGE satisfies the GNDS buidling blocks dataset, element and attribute

The table also includes the representation of these building blocks in FUDGE.

    +-----------+-------------+-----------+-----------------+-----------+-------+
    |    GNDS   | File System |     XML   |      JSON       |   HDF5    | FUDGE |
    +===========+=============+===========+=================+===========+=======+
    |  dataset  |    file     |    text   |      array      |  dataset  |       |
    +-----------+-------------+-----------+-----------------+-----------+-------+
    |  element  |  directory  |  element  |     object      +   group   |       |
    +-----------+-------------+-----------+-----------------+-----------+-------+
    | attribute |  meta-data  | attribute | not explicitely | attribute |       |
    |           |             |           | available       |           |       |
    +-----------+-------------+-----------+-----------------+-----------+-------+

The FUDGE documentation explains how the software meets the requirements for accessing GNDS files and the subsequent
processing of nuclear data for downstream applications. This description additionally incudes an overview of the 
relevant nuclear data theory and the GNDS format.

The General Nuclear Database Structure (GNDS)
*********************************************

.. todo::
   Add sections to provide an overview of GNDS

.. toctree::
   :maxdepth: 2

Background Theory
*****************

This section provides a theoretical background to the methods implemented in FUDGE.

.. todo::
   Expand the section as the review of the FUDGE source code continues.

.. toctree::
   :maxdepth: 2

   transferMatrices

The FUDGE Nuclear Data Platform
*******************************
FUDGE is primarily written in Python, with C/C++ extensions for CPU-intensive tasks. It's design allows for usage either
interactively or via batch scripts.

As a nuclear data managing tool, FUDGE includes the following capabilites:
  - Reading and writting GNDS files;
  - Tranlation from the legacy formats ENDF-6 and ENDL;
  - Methods to confirm conformance with the GNDS format and the underlying physics;
  - Plotting;
  - Processing for Monte Carlo and deterministic transport applications.

The sections below describe the FUDGE support for GNDS followed by a description of its functionality beyond GNDS

.. todo::
   Expand sections with more detail for "GNDS Representation in FUDGE" and a section for "Fudge beyond GNDS"

.. toctree::
   :maxdepth: 2

   gndsInFudge

Sample Applications
*******************

.. todo:: 
   Expand list of sample FUDGE applications

.. toctree::
   :maxdepth: 2

   exampleEndf6ToGnds

User Resources
**************

.. todo::
   Expand user resources to include installation instructions, links to more detailed GNDS description, etc.

.. toctree::
   :maxdepth: 2

   installation
   vis/multiplot


Literature References
*********************
Beck, B. R., and C. M. Mattoon. "FUDGE: A Toolkit for Nuclear Data Management and Processing." Transactions 110.1 (2014): 564-567.

Brown, D. A., B. R. Beck, J. L. Conlin, W. Haeck, C. M. Mattoon, and D. Wiarda. WPEC Subgroup-38 Final Report part II: 
Specifications for a new database structure. No. LLNL-TR-774621. Lawrence Livermore National Lab.(LLNL), Livermore, 
CA (United States), 2019.

Cullen, Dermott E. "Nuclear data preparation." Handbook of nuclear engineering. 2010.

Eisberg, Robert, and Robert Resnick. "Quantum physics of atoms, molecules, solids, nuclei, and particles." Quantum Physics of Atoms, 
Molecules, Solids, Nuclei, and Particles, 2nd Edition, by Robert Eisberg, Robert Resnick, pp. 864. ISBN 0-471-87373-X. Wiley-VCH, 
January 1985. (1985): 864.

Hedstrom, Gerald, Bret Beck, and Caleb Mattoon. Merced. No. Merced; 005182MLTPL00. Lawrence Livermore National Lab.(LLNL), 
Livermore, CA (United States), 2016.

Mattoon, C. M., B. R. Beck, N. R. Patel, N. C. Summers, G. W. Hedstrom, and D. A. Brown. Generalized nuclear data: a new 
structure (with supporting infrastructure) for handling nuclear data. Nuclear Data Sheets 113.12 (2012): 3145-3171.

Podgoršak, Ervin B. "Radiation physics for medical physicists." (2006).

Shultis, J. Kenneth, and Richard E. Faw. Radiation shielding. La Grange Park, Illinois: American Nuclear Society, 2000.

Trkov, Andre, and David A. Brown. ENDF-6 formats manual: Data formats and procedures for the evaluated nuclear data files.
No. BNL-203218-2018-INRE. Brookhaven National Lab.(BNL), Upton, NY (United States), 2018.
