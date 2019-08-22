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
