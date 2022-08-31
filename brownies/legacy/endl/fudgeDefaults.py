# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains default parameters used by fudge.

Variables::

    NeedPythonPathToFudge   Fudge attempts to determine if the location of fudge is included
                            in python's search path (e.g., locations set by PYTHONPATH).
                            If fudge determines that it is not in the python search path then
                            it adds its location - as given by DefaultFudgePath/Src - to the
                            python search path. However, if NeedPythonPathToFudge is false 
                            fudge will not add its location to the python search path.
    DefaultFudgePath        Location where fudge scripts and binaries can be found.
                                Scripts are to be in DefaultFudgePath/Src.
                                Binaries are to be in DefaultFudgePath/bin.
    NUCLEAR_DATABASE_DIR    Location of the nuclear data root.
    ENDL_DATABASE_DIR       Location of the default evaluated databases.
    TRANSLATED_DATABASE_DIR Location of the translated database
    ENDL_DEFAULT_DATABASE   Default evaluated database in ENDL_DATABASE_DIR.
                            This can be set with the envirnoment varialbe ENDLPATH.
    MCF_DATABASE_DIR        Location of the default mcf (Monte Carlo transport) processed files.
                            This can be set with the envirnoment varialbe MCFPATH.
    NDF_DATABASE_DIR        Location of the default ndf (deterministic transport) processed files.
                            This can be set with the envirnoment varialbe NDFPATH.

Also see the description of the variables PYTHONPATH and FUDGEPATH in the module fudge.py.
The module fudge.py describes the proper procedure to follow to override these parameters.
"""

import os

NeedPythonPathToFudge = True
DefaultFudgePath = '/usr/apps/fudge/current'
FUDGESRC_DIR = os.path.dirname( os.path.realpath( __file__ ) )
NUCLEAR_DATABASE_DIR = '/usr/gapps/data/nuclear'
ENDL_DATABASE_DIR = os.path.join( NUCLEAR_DATABASE_DIR, 'evaluated' )
if( not os.path.exists( ENDL_DATABASE_DIR ) ) : ENDL_DATABASE_DIR = os.path.join( NUCLEAR_DATABASE_DIR, 'endl_official' )
TRANSLATED_DATABASE_DIR = os.path.join( NUCLEAR_DATABASE_DIR, 'translated' )
ENDL_DEFAULT_DATABASE = "endl99"
MCF_DATABASE_DIR = "/usr/gapps/data/nuclear/processed/endl/endl94/mcf"
NDF_DATABASE_DIR = "/usr/gapps/data/nuclear/processed/current/ndf"
XENSL_DATABASE_DIR = "/usr/gapps/data/nuclear/development/xensl-RIPL-3"
