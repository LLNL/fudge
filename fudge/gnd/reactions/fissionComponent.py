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
This module contains the fissionComponent class
"""

from fudge.legacy.converting import endfFormats, endf_endl

import fudge
from . import base, reaction

__metaclass__ = type

class fissionComponent( reaction.reaction ) :
    """ special case of reaction, for storing 1st, 2nd, 3rd, etc fission chances when we also have total """

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, documentation = None, attributes = {} ) :
        """Creates a new reaction object of type outputChannel (that is, use 2-body outputChannel for 2-body reaction)."""

        reaction.reaction.__init__( self, outputChannel, label, ENDF_MT, crossSection, documentation, attributes )
        self.moniker = base.fissionComponentToken

    def isBasicReaction( self ) :

        return( False )

    def isCompleteReaction( self ):

        return( False )

    def check( self, info ) :
        warnings = self.__checkCrossSection__( info )
        if not info['crossSectionOnly']:
            warnings += self.__checkOutputChannel__( info )
        return warnings
