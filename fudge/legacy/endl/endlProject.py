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
The main content of this module is the class endlProject.

To resolve the location of the bdfls file, endlProject does the following,
    1) Uses the bdflsFileName argument if it is None; else,
    2) Uses "./bdfls" if it exists; else,
    3) Uses the bdfls in the database (at the same level as the "ascii" file) if it exists; else,
    4) Uses the bdfls pointed to by the environmental variable BDFLSPATH if it exists; else,
    5) Uses "/usr/gapps/data/nuclear/bdfls".
"""

import os
import tempfile
import string
import sys
import types

from fudge import fudgeDefaults, fudgeParameters
import bdfls
from fudge.core.utilities import brb
import fudgeDocumentationFile
import endlZA
import endl2
import endlmisc
import endl_Z

endlProjectList = []

def exit( ) :
    "Does some house cleaning and then exits python. This is the recommended way to exit fudge."

    n = len( endlProjectList ) - 1
    while ( n >= 0 ) :
        if ( endlProjectList[n].delWorkDirWhenDone ) : del endlProjectList[n]
        n -= 1
    sys.exit( )

def ll( ) :
    "Currently, the same as ls( )."

    ls( )

def ls( ) :
    """List instantiated databases."""

    j = 0
    for i in endlProjectList : print "   endlProjectList[%d] = %s" % ( j, i.database ); j += 1

def removeProject( p, delWorkDirWhenDone = None ) :
    "Removes endlProject p from Fudge's internal list and returns None."

    if ( type( p ) == types.InstanceType ) and ( p.__class__.__name__ == "endlProject" ) :
        if ( delWorkDirWhenDone != None ) : p.delWorkDirWhenDone = delWorkDirWhenDone
        r = range( len( endlProjectList ) )
        for i in r :
            if ( p == endlProjectList[i] ) :
                del endlProjectList[i]
                return None
        raise Exception( "\nError in removeProject: project not found in endlProjectList" )
    raise Exception( "\nError in removeProject: object is not an endlProject instance" )

class endlProject :

    def __init__( self, database = None, projectile = None, yi = None, workDir = None, delWorkDirWhenDone = 0, readOnly = 0, release = 'current',
        bdflsFile = None, cleanWorkDirOfZAs = False ) :
        """Creates a new endlProject with database as the default database to use when reading ZAs. If bdflsFile is not None, then
        it is the bdfls class instance associated with self."""

        if( yi != None ) :
            endlmisc.printWarning( 'Argument "yi" for Constructor of endlProject is deprecated, please use argument "projectile" instead.' )
            if( projectile != None ) : raise Exception( 'Both "yi" and "projectile" argument specified, only use "projectile".' )
            projectile = yi
        if( database == None ) : database = fudgeDefaults.ENDL_DEFAULT_DATABASE
        if( bdflsFile == None ) :
            bdflsFileName = None
            if( os.path.exists( 'bdfls' ) ) : bdflsFileName = './bdfls'    # Else logic handled below.
        else :
            bdflsFileName = bdflsFile.Source
        self.readOnly = readOnly or fudgeParameters.ReadOnly
        self.delWorkDirWhenDone = delWorkDirWhenDone
        dbInDatabaseName, yiInDatabaseName = os.path.split( database )
        if( yiInDatabaseName in [ "yi01", "yi02", "yi03", "yi04", "yi05",  "yi06", "yi07", "yi09" ] ) :
            database = dbInDatabaseName
            if( ( projectile != None ) and ( yiInDatabaseName != endlmisc.incidentParticleTags( projectile )[1] ) ) :
                raise Exception( "\nError in endlProject.__init__: yi in database name and yi differ" )
            projectile = endlmisc.incidentParticleTags( yiInDatabaseName )[5]
        else :
            if( projectile == None ) : projectile = 'n'
        self.yiTags = endlmisc.incidentParticleTags( projectile )
        self.name = self.yiTags[5]
        ( self.database, Status ) = endlmisc.getFullPath( os.path.join( database, 'ascii' ), fudgeDefaults.ENDL_DATABASE_DIR )
        if( Status == 0 ) : ( self.database, Status ) = endlmisc.getFullPath( database, fudgeDefaults.ENDL_DATABASE_DIR )
        if( Status == 0 ) : ( self.database, Status ) = endlmisc.getFullPath( os.path.join( database, 'ascii' ), \
            os.path.join( fudgeDefaults.TRANSLATED_DATABASE_DIR, release ) )
        if ( Status == 0 ) :                                # File or directory does not exist
            if ( self.yiTags[0] not in [ 1, 2, 3, 4, 5, 6, 7, 9 ] ) : raise Exception( "\nError in endlProject.__init__: invalid projectile = %s" % self.name )
            self.database = None
            self.yi = self.yiTags[0]
            if( fudgeParameters.VerboseMode > 0 ) : print self.yiTags[1]
        elif ( Status == 1 ) :
            raise Exception( "\nError in endlProject.__init__: database is not a directory" )
        else :
            self.database = os.path.join( self.database, self.yiTags[1] )
            if ( not os.path.exists( self.database ) ) :
                self.database = None
            else :
                if( bdflsFileName == None ) :
                    bdflsFileName = os.path.join( os.path.dirname( self.database ).split( '/ascii' )[0], 'bdfls' )
                    if( not os.path.exists( bdflsFileName ) ) : bdflsFileName = None
            self.yi = self.yiTags[0]
        self.bdflsFile = bdflsFile
        if( bdflsFile == None ) : self.bdflsFile = bdfls.getBdflsFile( name = bdflsFileName )   # endlProject relies on bdfls class to search BDFLSPATH 
        self.bdflsFileName = endlmisc.getAbsPath_withDefaultDBCheck( self.bdflsFile.Source )[0] # and "/usr/gapps/data/nuclear/bdfls", in that order,
        self.mass = 0                                                                           # for the bdfls file.
        yoZA = endl2.yoToZA( self.yiTags[0] )
        if( self.yiTags[0] != 7 ) : self.mass = self.bdflsFile.mass( abs( yoZA ) )
        if( yoZA < 0 ) : yoZA = 0
        self.Z, self.A = endl2.ZandAFromZA( yoZA )
        if ( fudgeParameters.VerboseMode > 0 ) :
            print 'Using bdfls', self.bdflsFileName
            print self.database
        gappsDir = os.path.realpath( '/usr/gapps/data/nuclear' )
        if self.readOnly :
            self.workDir = None
        else :
            if( workDir == '' ) : workDir = self.database
            if ( workDir == None ) :
                self.workDir = tempfile.mkdtemp( )
            else :
                self.workDir = workDir
            if( gappsDir == os.path.realpath( self.workDir )[:len( gappsDir )] ) : 
                if( 'development' not in os.path.realpath( self.workDir ) ) :
                    raise Exception( 'cannot make workDir point to /usr/gapps/data/nuclear, except in a development sub-directory' )
            if( not os.path.isdir( self.workDir ) ) : os.makedirs( self.workDir )
            self.workDir = endlmisc.appendyiDirIfNeeded( self.workDir, self.yiTags[0] )
            if( not os.path.isdir( self.workDir ) ) : os.makedirs( self.workDir )
        if( not( self.workDir is None ) and cleanWorkDirOfZAs ) : os.system( 'rm -rf %s/za[0-9][0-9][0-9][0-9][0-9][0-9]*' % self.workDir )
        self.zas = []
        self.documentation = None
        if( self.database != None ) :
            documentationFileName = os.path.join( self.database, 'documentation.txt' )
            if( os.path.exists( documentationFileName ) ) : self.documentation = fudgeDocumentationFile.fudgeDocumentationFile( documentationFileName )
        endlProjectList.append( self )

    def __contains__( self, target ) :
        """Returns True if target is in self database, otherwise returns False."""

        strZA = endlmisc.strZASuffix( target )
        return( strZA in self.ZAList( ) )

    def __del__( self ) :
        """For internal use only.
        Does house cleaning on self.  Note, self is not deleted until all references to it are removed."""

        if( hasattr( self, 'zas' ) ) :
            while ( len( self.zas ) > 0 ) : del self.zas[0]
        if( hasattr( self, 'workDir' ) ) :
            if( type( self.workDir ) != type( None ) ) and ( self.delWorkDirWhenDone ) and ( os.path.exists( self.workDir ) ) : 
                os.system( "rm -rf " + self.workDir )

    def __getitem__( self, index ) :
        """Returns the ZA at self's zas[index]."""

        return( self.zas[index] )

        
    def __len__( self ) :
        "Returns the number of ZAs instantiated within self."

        return( len( self.zas ) )

    def addZA( self, ZA, suffix = "" ) :
        """Adds a new empty endlZA class to self and returns a reference to it.  Also see readZA."""

        sZA = endlmisc.strZASuffix( ZA, suffix )
        for z in self.zas :
            if ( z.sZA == sZA ) : raise Exception( "\nError in endlProject.addZA: za = %s already present in this endlProject" % sZA )
        z = endlZA.endlZA( ZA, self.yi, database = None, workDir = self.workDir, suffix = suffix, readOnly = self.readOnly, bdflsFile = self.bdflsFile )
        self.zas.append( z )
        return z

    def getDocumentation( self ) :
        """Return None if there is not documentation file, otherwise, the text of the documentation file is returned."""

        if( self.documentation == None ) : return( None )
        return( self.documentation.getText( ) )

    def getMass( self ) :
        """Returns the mass of the projectile in MeV."""

        return( self.mass )

    def info( self ) :
        """Prints yi, the database directory name, za sub-directory name, work directory name and delWorkDirWhenDone flag for self."""

        print "yi = %2d  symbol = %s  name = %s" % ( self.yi, self.yiTags[2], self.yiTags[3] )
        print "Database directory =", self.database
        print "Working directory =", self.workDir
        print "delWorkDirWhenDone =", self.delWorkDirWhenDone

    def ll( self, Z = None, A = None, suffix = "", w = 5 ) :
        """Prints a table of all ZA directories in self's default source database that match Z and A.  The table is printed with w columns."""

        try :
            zalist = self.ZAList( Z = Z, A = A, suffix = suffix )
        except :
            zalist = [ "No default database" ]
        l = 0;
        for za in zalist : l = max( l, len( za ) )
        Fmt = "%%-%ds" % l
        zal = []
        for za in zalist : zal.append( Fmt % za )
        brb.tylist( zal, w )
        print
        self.ls( )

    def ls( self ) :
        """Prints a list of ZAs instantiated in self."""

        if ( len( self.zas ) > 0 ) :
            print "            ZA        Source"
            print "--------------------------------------------------------------------"
            for i in range( len( self.zas ) ) : print "zas[%3d] %-11.11s  %s" % ( i, `self.zas[i].sZA`, self.zas[i].source )

    def process( self, oldfile, newfile = None, options = "", defines = [], bdflsFile = None, ndfgen = None, mcfgen = None, endep = None, extraZAs = [] ) :
        """Processes the active ZAs in self, using oldfile as a template to 
create newfile.  Active ZAs are ones that have been read-in-to (readZA) 
or added-to (addZA) self.  Uses bdfls pointed to by bdflsFile if not None; otherwise,
used search in the standard algorithm. The options argument sets internal flags (see below).
The defines argument must contain a list of ksh environmental assignments (e.g.,
defines = ( 'greeting="Hello world"', 'NDFGENINPUTSPATH=./ndf.inputs' ), note that
since "Hello world" is 2 words, it is quoted).

Valid options are::
    -i   Creates a ".index" file (this should generally be included),
    -u   Update: newfile is oldfile amended by the active ZAs, 
         otherwise newfile only contains the active ZAs in self,
    -v   Verbose mode (additional information is printed to screen and 
         intermediate files are kept,
    -e   Stop processing after endep has been called (mainly for internal use).

All options can be negated by prepending with 'no_' (e.g., -no_u').
"""

        if( bdflsFile == None ) : bdflsFile = self.bdflsFile.Source
        if ( self.workDir == None ) :
            raise Exception( "\nError in endlProject.process: project does not have a work directory." )
        else :
            if ( not os.path.exists( oldfile ) ) :
                if ( string.rfind( oldfile, "/" ) == -1 ) :
                    f = os.path.join( fudgeDefaults.NDF_DATABASE_DIR, oldfile )
                    if ( not os.path.exists( f ) ) : f = os.path.join( fudgeDefaults.MCF_DATABASE_DIR, oldfile )
                    if ( os.path.exists( f ) ) : oldfile = f
            if ( not os.path.exists( oldfile ) ) :
                raise Exception( "\nError in endlProject.process: old file does not exist (%s)" % oldfile )
            if ( newfile == None ) :
                i = string.rfind( oldfile, "/" ) + 1        # rfind returns -1 if "/" not found. This is great.
                newfile = oldfile[i:] + ".new"
            endlmisc.processDataBase( self.workDir, oldfile, newfile, options, defines = defines, bdflsFile = bdflsFile, ndfgen = ndfgen,
                mcfgen = mcfgen, endep = endep, extraZAs = extraZAs )

    def readZA( self, ZA, database = None, suffix = "" ) :
        """Reads in ZA from database, creating a new endlZA in self.  If database = None, the database from self is used. Also see addZA."""

        sZA = endlmisc.strZASuffix( ZA, suffix )
        for z in self.zas : 
            if ( z.sZA == sZA ) :
                endlmisc.printWarning( "za = %s already present in this endlProject" % sZA )
                return z
        if ( database == None ) : database = self.database
        z = endlZA.endlZA( ZA, self.yi, database = database, workDir = self.workDir, suffix = suffix, readOnly = self.readOnly, bdflsFile = self.bdflsFile )
        self.zas.append( z )
        return z

    def read( self, target, database = None ) :

        if( target[:18] == 'FissionProductENDL' ) :
            ZAs = target[18:]
        else :
            Z, A, suffix, ZA = endl2.gndNameToEndlZ_A_Suffix( target )
            ZAs = endlmisc.strZASuffix( 1000 * Z + A, suffix )
        return( self.readZA( ZAs, database = database ) )

    def remove( self, target ) :
        """Removes self's reference to target and deletes target's directory from self's work 
        directory (i.e., if process( ) is called this target will not be processed).  Also see unread."""

        for i in range( len( self.zas ) ) :
            if ( self.zas[i] == target ) :
                if ( self.zas[i].workDir != None ) : os.system( "rm -rf " + self.zas[i].workDir )
                del self.zas[i]
                return
        endlmisc.printWarning( "Warning in endlProject.remove: target = %s is not a part of project" % target.name )

    def removeZA( self, ZA, suffix = "" ) :
        """Removes self's reference to ZA and deletes ZA's directory from self's work 
        directory (i.e., if process( ) is called this ZA will not be processed).  Also see unreadZA."""

        sZA = endlmisc.strZASuffix( ZA, suffix )
        for i in range( len( self.zas ) ) :
            if ( self.zas[i].sZA == sZA ) :
                if ( self.zas[i].workDir != None ) : os.system( "rm -rf " + self.zas[i].workDir )
                del self.zas[i]
                return
        endlmisc.printWarning( "Warning in endlProject.removeZA: za = %s not in project" % sZA )

    def save( self ) :
        "Saves all ZAs instantiated within self."

        for z in self.zas : z.save( )

    def saveProcess( self, oldfile, newfile = None, options = " -i -u " ) :
        """Saves all ZAs instantiated within self and then calls process( oldfile, newfile, options )."""

        self.save( )
        self.process( oldfile, newfile, options )


    def targetList( self, Symbol = None, A = None, suffix = None ) :
        """Returns a list of target sub-directories in self's database matching Symbol, A and suffix."""

        if( self.database == None ) : raise Exception( "No default database defined." )
        targetList = []
        targets = os.listdir( self.database )
        targets.sort( )
        if( ( type( Symbol ) == type( None ) ) or ( type( Symbol ) == type( '' ) ) ) : Symbol = [ Symbol ]
        if( ( type( A ) == type( None ) ) or ( type( A ) == type( 0 ) ) ) : A = [ A ]
        if( ( type( suffix ) == type( None ) ) or ( type( suffix ) == type( '' ) ) ) : suffix = [ suffix ]
        for target in targets :
            if ( endlmisc.validZADirectoryName_( target ) ) :
                S_ = None
                if( Symbol != [ None ] ) : S_ = endl_Z.endl_ZSymbol( int( target[2:5] ) )
                A_ = None
                if( A != [ None ] ) : A_ = int( target[5:8] )
                s_ = None
                if( suffix != [ None ] ) : s_ = target[8:]
                if( ( S_ in Symbol ) and ( A_ in A ) and ( s_ in suffix ) ) :
                    d = os.path.join( self.database, target )
                    if ( os.path.exists( d ) ) : targetList.append( endl2.endlToGNDName( target ) )
        return targetList

    def unread( self, target ) :
        """Removes self's reference to target, but does not delete target's directory from self's work 
        directory (i.e., if process( ) is called this target will be processed). Also see remove."""

        for i in range( len( self.zas ) ) :
            if ( self.zas[i] == target ) :
                del self.zas[i]
                return
        endlmisc.printWarning( "Warning in endlProject.unread: target = %s is not a part of project" % target )

    def unreadZA( self, ZA, suffix = "" ) :
        """Removes self's reference to ZA, but does not delete ZA's directory from self's work 
        directory (i.e., if process( ) is called this ZA will be processed). Also see removeZA."""

        sZA = endlmisc.strZASuffix( ZA, suffix )
        for i in range( len( self.zas ) ) :
            if ( self.zas[i].sZA == sZA ) :
                del self.zas[i]
                return
        endlmisc.printWarning( "Warning in endlProject.unreadZA: za = " + sZA + " not in project" )

    def ZA( self, ZA, suffix = "" ) :
        """Returns a reference to the endlZA matching ZA in self.  ZA is passed to endlmisc.intZA 
        and its return value to used for matching ZAs. If ZA has not been instantiated then None is returned."""

        sZA = endlmisc.strZASuffix( ZA, suffix )
        for z in self.zas :
            if ( z.sZA == sZA ) : return z 
        return None

    def ZAList( self, Z = None, A = None, suffix = "" ) :
        """Returns a list of ZA sub-directories in self's database matching Z and A."""

        if( self.database == None ) : raise Exception( "\nError in endlProject.ZAList: no default database defined." )
        Zb = 2
        Ze = 2
        ZStr = ""
        Ab = 5
        Ae = 5
        AStr = ""
        if ( Z != None ) :
            ZStr = "%3.3d" % Z
            Ze = 5
        if ( A != None ) :
            AStr = "%3.3d" % A
            Ae = 8
        zalist = []
        ZAs = os.listdir( self.database )
        ZAs.sort( )
        for za in ZAs :
            if ( endlmisc.validZADirectoryName_( za ) ) :
                if ( ZStr == za[Zb:Ze] ) and ( AStr == za[Ab:Ae] ) and ( ( suffix == "" ) or ( suffix == za[8:] ) ) :
                    d = os.path.join( self.database, za )
                    if ( os.path.exists( d ) ) : zalist.append( za )
        return zalist
