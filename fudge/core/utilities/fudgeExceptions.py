# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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
