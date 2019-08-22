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
This module contains the class fudgeTempFile which simplies the functions in
the tempfile module.
"""

import os
import tempfile

__metaclass__ = type

class fudgeTempFile :
    """This class creates a temporary file using tempfile.mkstemp and supports common functionallity
    for the file (e.g., write). Currently, reading from the temporary file is not supported."""

    def __init__( self, prefix = None, suffix = "", dir = None, deleteOnClose = False ) :
        """Contructor for the fudgeTempFile class. Calls tempfile.mkstemp to open a temporary file.
        If deleteOnClose is 'True', file will be deleted when it is close."""

        if( dir == None ) :
            dir = tempfile.gettempdir( )
        else :
            if( not os.path.exists( dir ) ) : os.makedirs( dir )
        if( prefix == None ) : prefix = tempfile.gettempprefix( )

        self.fd, self.name = tempfile.mkstemp( suffix = suffix, prefix = prefix, dir = dir )
        self.deleted = False
        self.deleteOnClose = deleteOnClose

    def __del__( self ) :
        """If class instance is deleted, the file is proprely closed."""

        if( self.fd != None ) : self.close( )

    def close( self, raiseIfClosed = True ) :
        """Closes the file if still opened. If raiseIfClosed is 'True' and file is already
        closed, a raise is executed."""

        if( self.fd != None ) :
            os.close( self.fd )
            self.fd = None
            if( self.deleteOnClose ) : self.delete( )
        elif( raiseIfClosed ) :
            raise Exception( 'Error from fudgeTempFile.close: file already closed' )

    def delete( self ) :
        """Deletes file if it still exist.  If required, this method calls close first. If file has already
        been deleted, a raise is executed."""

        if( self.deleted ) :
            raise Exception( 'Error from fudgeTempFile.delete: file already deleted' )
        else :
            if( self.fd != None ) : self.close( )
            os.remove( self.name )
            self.deleted = True

    def getName( self ) :
        """Returns self's name which is the full path name of the temporary file."""

        return( self.name )

    def getFileDescriptor( self ) :
        """Returns self's file descriptor."""

        return( self.fd )

    def isOpened( self ) :
        """Returns 'True' if the file is still opened."""

        return( self.fd != None )

    def isDeleted( self ) :
        """Returns 'True' if the file has been deleted."""

        return( self.deleted )

    def write( self, str ) :
        """Write str to file. If file is closed or not all characters were written then a raise is executed."""

        if( self.fd != None ) :
            n = os.write( self.fd, str )
            if( n != len( str ) ) : raise Exception( 'Error from fudgeTempFile.write: only %d of %d characters written' % ( n, len( str ) ) )
        else :
            raise Exception( 'Error from fudgeTempFile.write: attempted to write to a closed file.' )
