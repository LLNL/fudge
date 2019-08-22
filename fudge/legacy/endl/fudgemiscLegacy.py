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

def databaseList( yi = None, path = None ) :
    """Returns a list of the default databases. The list can be narrowed to a single incident particle
    type if yi is specified. Also, the location where databases are to be searched for can be set with path.
    This allows a user to have their own directory containing databases."""

    def findYis( yiPath, yi, list ) :
        """For internal use only."""

        yis = []
        if ( os.path.isdir( yiPath ) ) :
            ls = os.listdir( yiPath )
            ls.sort( )
            for f in ls :
                if( re.match( r"yi\d\d$", f ) ) :
                    y = os.path.join( yiPath, f )
                    if( os.path.isdir( y ) ) :
                        iyi = int( f[2:] )
                        if( iyi in [ 0, 1, 2, 3, 4, 5, 6, 7, 9 ] ) :
                            if( ( yi == None ) or ( yi == iyi ) ) : yis.append( f )
        if( len( yis ) > 0 ) : list.append( [ db, yis ] )

    listTranslated = path == None
    if( yi != None ) : yi = endlmisc.incidentParticleTags( yi )[0]
    if ( path == None ) : path = fudgeDefaults.ENDL_DATABASE_DIR
    l = os.listdir( path )
    l.sort( )
    ll = []
    for db in l :
        d = os.path.join( path, db )
        a = os.path.join( d, 'ascii' )
        if( os.path.exists( a ) ) : d = a
        findYis( d, yi, ll )
    if( listTranslated ) :
        l = os.listdir( os.path.join( fudgeDefaults.TRANSLATED_DATABASE_DIR, 'current' ) )
        for db in l :
            translated = os.path.join( fudgeDefaults.TRANSLATED_DATABASE_DIR, 'current', db, 'ascii' )
            findYis( translated, yi, ll )
    return ll

def projectileList( path = None ) :
    """Returns a list of the projectiles in path."""

    def getYis( database ) : 

        yis = []
        ls = os.listdir( database )
        ls.sort( )
        for f in ls :
            if( re.match( r"yi\d\d$", f ) ) :
                y = os.path.join( database, f )
                if( os.path.isdir( y ) ) :
                    iyi = int( f[2:] )
                    if( iyi in [ 0, 1, 2, 3, 4, 5, 6, 7, 9 ] ) : yis.append( f )
        return( yis )

    if( path == None ) : path = fudgeDefaults.ENDL_DEFAULT_DATABASE
    database, Status = endlmisc.getFullPath( path, fudgeDefaults.ENDL_DATABASE_DIR )
    yis = getYis( database )
    if( len( yis ) == 0 ) : yis = getYis( os.path.join( database, "ascii" ) )
    return yis

def info( yi = None, w = 2 ) :
    """List the default database path and all databases in it. The list of databases is printed with w columns."""

    bdflsFile = bdfls.getDefaultBdfls( )
    print "bdfls file =", bdflsFile.Source
    print "Default database path =", fudgeDefaults.ENDL_DATABASE_DIR
    dbs = databaseList( yi = yi )
    lll = []
    for d in dbs :
        for yi in d[1] : lll.append( d[0] + "/" + yi )
    if ( len( lll ) > 0 ) : brb.tlist( lll, w )
