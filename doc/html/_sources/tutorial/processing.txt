Processing GND files
====================

``reconstructResonances()``
``process()``

Options for processing
----------------------

    -h                  help
    -o outFile          specify the output file to be ``outFile``
    -r                  perform resonance reconstruction
    -a                  compute angular distributions from the RRR 
                        (Reich-Moore Limited only)
    -e                  compute average energy and forward linear momentum depositions
    -l                  linearize cross sections and outgoing particle distributions
    -H                  heat target in evaluation to room temperature 
                        (override with --temp option)
    --temp T            override the heating temperature with value T.  
                        T must be a string with format "number units", e.g. "293.15 K".
    -u                  compute URR probability tables
    -x                  compute transfer matrices
    -g                  group all cross sections, outgoing distributions and covariances
    --group groupFile   read the group boundaries from the file ``groupFile``.  The first
                        line of this file is a single string specifying the units of 
                        the group boundaries (e.g. eV).  The group boundaries must be 
                        listed after with one energy per line.
    --outFormat str     specify the output format for the processed data.  May be any of: 
                        [ 'gnd', 'gendf', 'ace', 'pendf', 'endf' ].  Note, not all output
                        format support all processing options, so things may silently get 
                        lost.
