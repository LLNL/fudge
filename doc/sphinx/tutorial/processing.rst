Processing GNDS files
=====================

The ``processProtare`` script (in fudge/bin/processProtare.py) is the main driver for processing
nuclear data evaluations. This script can be used to reconstruct resonances (if necessary),
heat reaction cross sections to various temperatures, generate group-averaged cross sections and transfer matrices,
and generate average product energy and momentum deposition data.

The name ``processProtare`` includes the short-hand 'ProTarE', which stands for "Projectile, Target and Evaluation".
Each GNDS ``reactionSuite`` contains this combination: an evaluation of all reactions involving a combination of
projectile and target.

Processed data generated during processProtare are stored internally in Fudge data structures, and once processing
is complete the resulting file is written back out as a GNDS/XML file.

Options for processing
----------------------

processProtare supports many command-line options. Perhaps the most useful is the '-h' or '--help' option:

::

        >python processProtare.py -h
        usage: processProtare.py [-h] [--tag TAG] [-o OUTPUT] [--writeConvertedUnits]
                                 [--energyUnit ENERGYUNIT] [-t TEMPERATURES]
                                 [--temperatureUnit TEMPERATUREUNIT] [--bdfls BDFLS]
                                 [--fluxID FLUXID] [-g GID]
                                 [--legendreMax LEGENDREMAX] [-mg] [-mc] [-up]
                                 [--reconstruct {all,crossSection,angular}]
                                 [--threads THREADS] [-v]
                                 [--reconAccuracy RECONACCURACY]
                                 [--prefixMultiGroup PREFIXMULTIGROUP]
                                 [--prefixMonteCarlo PREFIXMONTECARLO]
                                 [--prefixUpscatter PREFIXUPSCATTER]
                                 [--prefixHeated PREFIXHEATED] [--prefixAep PREFIXAEP]
                                 [--prefixEval PREFIXEVAL] [--prefixRecon PREFIXRECON]
                                 gnds

        Processes all data in a GNDS file.

        positional arguments:
          gnds                   gnds file to process

        optional arguments:
          -h, --help            show this help message and exit
          ...

Sample usage
------------

processProtare takes a GNDS file as input. It also currently requires a LLNL-specific file called 'bdfls', which contains a database of
group boundaries and fluxes.

After translating an ENDF file into GNDS, you can do simple processing by:

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls

This performs basic processing: Doppler-broadens cross sections to 300 K (2.586e-8 MeV/k) and computes average
energy and momentum deposited into reaction products. It writes a new processed file with the extension ".proc.xml" (i.e., if
the original file was named n-001_H_001.xml, this will output n-001_H_001.proc.xml).
This processing adds two new 'style' elements to the GNDS file: ``heated`` corresponding to the Doppler broadening, and ``averageProductData``
corresponding to the energy/momentum calculations.

Below are some examples of more advanced processing options.

Heating cross sections to multiple temperatures:

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls -t 300 -t 600 -t 900
                --temperatureUnit K 

Using eV instead of the default MeV for energies:

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls --energyUnit eV

Heating to two temperatures and generating cumulative distribution functions (CDFs) for faster sampling in Monte Carlo transport codes:

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls -t 300 -t 600
                --temperatureUnit K --energyUnit eV -mc
        
Generating multi-group cross sections and transfer matrices (requires that Merced be compiled):

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls -t 300 -t 600
                --temperatureUnit K --energyUnit eV -mg

By default, processProtare uses LLNL's 230-group energy boundaries for neutrons,
41 groups for photons and 64 groups for light charged particles. Other groups may be selected using the -g option.
For example, to select LLNL's 87-group structure (gid 4) for neutrons:

::

        > python processProtare.py GNDS_file_name.xml --bdfls <path_to_fudge>/fudge/legacy/endl/bdfls -t 300 -t 600
                --temperatureUnit K --energyUnit eV -mg -g n=LLNL_gid_4

Group boundaries are currently defined inside the bdfls file, although we plan to update that to a more user-friendly format soon.

To-do list
----------

processProtare is not yet complete. Some functionality yet to be implemented:

- thermal upscattering corrections for elastic transfer matrices
- processing thermal neutron scattering data including S(alpha,beta)
- generating probability tables from unresolved resonance parameters
- reconstructing detailed angular distributions from resonance parameters
- processing covariances
- and likely other features

