# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
