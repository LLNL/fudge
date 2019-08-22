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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

        if( dir is None ) :
            dir = tempfile.gettempdir( )
        else :
            if( not os.path.exists( dir ) ) : os.makedirs( dir )
        if( prefix is None ) : prefix = tempfile.gettempprefix( )

        self.fd, self.name = tempfile.mkstemp( suffix = suffix, prefix = prefix, dir = dir )
        self.deleted = False
        self.deleteOnClose = deleteOnClose

    def __del__( self ) :
        """If class instance is deleted, the file is proprely closed."""

        if( self.fd is not None ) : self.close( )

    def close( self, raiseIfClosed = True ) :
        """Closes the file if still opened. If raiseIfClosed is 'True' and file is already
        closed, a raise is executed."""

        if( self.fd is not None ) :
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
            if( self.fd is not None ) : self.close( )
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

        return( self.fd is not None )

    def isDeleted( self ) :
        """Returns 'True' if the file has been deleted."""

        return( self.deleted )

    def write( self, str ) :
        """Write str to file. If file is closed or not all characters were written then a raise is executed."""

        if( self.fd is not None ) :
            n = os.write( self.fd, str )
            if( n != len( str ) ) : raise Exception( 'Error from fudgeTempFile.write: only %d of %d characters written' % ( n, len( str ) ) )
        else :
            raise Exception( 'Error from fudgeTempFile.write: attempted to write to a closed file.' )
