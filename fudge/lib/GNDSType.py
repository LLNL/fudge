# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the function **type** for determining the type (i.e., *reactionSuite*, *covarianceSuite* or *PoPs*) 
of a **GNDS/XML** file and the function **read** for reading into **FUDGE** a *GNDS/XML* file.

If the file is not a *GNDS/XML* file, a raise is executed.
"""

import xml.sax

from fudge.gnds import reactionSuite as reactionSuiteModule
from fudge.gnds import thermalScattering as thermalScatteringModule
from fudge.gnds.covariances import covarianceSuite as covarianceSuiteModule
from PoPs import database as databaseModule

class GNDSTypeException( Exception ) :
    """
    For internal use. The raise executed in **GNDSTypeHandler.startElement** if the file is a valid
    *GNDS/XML* file as there is not other way to stop the parser.
    """

    pass

class GNDSTypeHandler( xml.sax.handler.ContentHandler ) :
    """For internal use. The SAX handler used to determine the *GNDS/XML* file type."""

    def startElement( self, name, attributes ) :
        """
        The SAX handler's startElement. This method always raises on the first start element to stop the parser.
        If the exception GNDSTypeException is raise, a valid tag name for the first element was found. All
        other exceptions are an error.
        """

        self.name = name
        if( name in [ reactionSuiteModule.reactionSuite.moniker, covarianceSuiteModule.covarianceSuite.moniker ] ) :
            self.data = { 'projectile' : attributes['projectile'], 
                          'target'     : attributes['target'], 
                          'evaluation' : attributes['evaluation'] }
        elif( name == databaseModule.database.moniker ) :
            self.data = { 'name' : attributes['name'] }
        elif( name == thermalScatteringModule.thermalScattering.moniker ) :
            self.data = { 'material': attributes['material'] }
        else :
            raise TypeError( 'Invalid XML file with top element = "%s"' % name )
        raise GNDSTypeException( )

def type( fileName, show = False ) :
    """
    This function determines the type of the *GNDS/XML* file given by *fileName* or raises an error if
    the file is not a *GNDS/XML* file. For a *GNDS/XML* file this function returns a tuple of length 2
    items. The first item is the moniker for the *GNDS/XML* type. The second item is a tuple whose contents
    depend of the *GNDS/XML* file type. If the type is *reactionSuite* or *covarianceSuite*, then the tuple
    contains 3 items. The first is the projectile's ID, the second is the target's ID and the third is the 
    evaluation string.  If the type is *PoPs*, the tuple contains 1 item which is the database's name.

    If the file is not a *GNDS/XML* file, a raise is executed.
    """

    parser = xml.sax.make_parser( )
    handler = GNDSTypeHandler( )
    parser.setContentHandler( handler )
    try :
        parser.parse( fileName )
    except GNDSTypeException :
        if( show ) : print('%-16s %s: %s' % (handler.name, fileName, handler.data))
        return( handler.name, handler.data )
    except xml.sax._exceptions.SAXParseException :
        if( show ) :
            print ('%-20s ERROR: %s' % ( "INVALID XML", fileName ))
        else :
            raise
    except :
        if( show ) :
            print ('%-20s ERROR: %s' % ("Oops", fileName))
        else :
            raise

def read( fileName, reactionSuite = None ) :
    """
    This function uses the function **type** to determine the proper **FUDGE** function to call to read 
    in the file *fileName* into **FUDGE**.  It returns the **FUDGE** instance for the type.
    """

    name, dummy = type( fileName )
    if( name == reactionSuiteModule.reactionSuite.moniker ) :
        return( reactionSuiteModule.readXML( fileName ) )
    elif( name == covarianceSuiteModule.covarianceSuite.moniker ) :
        return( covarianceSuiteModule.readXML( fileName, reactionSuite ) )
    elif( name == thermalScatteringModule.thermalScattering.moniker ) :
        return( thermalScatteringModule.readXML( fileName ) )
    return( databaseModule.database.readFile( fileName ) )
