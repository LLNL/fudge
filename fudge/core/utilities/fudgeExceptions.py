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
This module contains the fudge exception classes.
"""

class FUDGE_Exception( Exception ) :
    """This class is raised whenever fudge detects an error that is not handled by
    another FUDGE/ENDL exception."""


    def __init__( self, value ) :
        """This constructor should be called as 

>>> raise FUDGE_Exception( "Put exception string here" )."""

        self.value = value

    def __str__( self ) :
        """Returns the string pass to the __init__ method."""

        return repr( self.value )

class ENDL_CheckerException( FUDGE_Exception ) :
    """This class is raised whenever a check method detects bad ENDL data."""

class ENDL_numpyException( Exception ) :
    """This class is raised whenever an error is detected with the numpy module."""

class ENDL_DesignerIsotopeReadException( Exception ) :
    """This class is raised whenever read is called on a designer isotope."""

class ENDL_addFileException_FileExist( FUDGE_Exception ) :
    """This class is raised whenever endlZA.addFile is called and the file already exists."""

class ENDL_fixHeadersException_NoC10Data( FUDGE_Exception ) :
    """This class is raised whenever endlZA.fixHeaders needs missing header data and cannot find the 
    C = 10, I = 0 which it assumes as valid header data."""
