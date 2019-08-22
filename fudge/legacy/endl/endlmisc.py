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
This module contains miscellaneous functions related to ENDL.
"""

import os
import sys
import copy
import string
import math
import re
from fudge import fudgeDefaults

def copyEndlHeader( Desc, Src ) :
    "This needs to be documented."

    if( not hasattr( Src, "h" ) ) :
        raise Exception( "\nError from copyEndlHeader: Source does not have a header" )
    Desc.h = copy.deepcopy( Src.h )
    if( Desc.I >= 0 ) :
        Desc.h[1] = "%s%3d%s" % ( Desc.h[1][:2], Desc.I, Desc.h[1][5:] )
    Desc.C = Src.C
    Desc.Q = Src.Q
    Desc.S = Src.S
    Desc.X1 = Src.X1
    Desc.X2 = Src.X2
    Desc.X3 = Src.X3
    Desc.X4 = Src.X4
    Desc.yo = Src.yo

def incidentParticleTags( yi ) :
    """Returns a list of [ iYi, strYi, symbol, name ] for projectile yi where yi 
    can be any of iYi, strYi, label or name. iYi is the projectile integer id, strYi
    is a string of the form "yi##" where ## is iYi converted to a two digit string,
    and the last two (symbol and name) are the projectiles symbol and name. For
    example, incidentParticleTags( "yi02" ) would return [ 2, "yi02", "p", "proton" ]."""

    yis = [ [ 0, "yi00", "u",  "unknown",  "unknown", "unknown"  ], [ 1, "yi01", "n",   "neutron",  "n",   "n"      ],
            [ 2, "yi02", "p",  "proton",   "H1",     "H1"      ], [ 3, "yi03", "d",   "deuteron", "H2",   "H2"      ],
            [ 4, "yi04", "t",  "triton",   "H3",     "H3"      ], [ 5, "yi05", "He3", "helium3",  "He3",  "He3"     ],
            [ 6, "yi06", "a",  "alpha",    "He",      "He4"     ], [ 7, "yi07", "g",   "gamma",    "gamma", "gamma"    ],
            [ 8, "yi08", "e+", "positron", "b+",      "positron" ], [ 9, "yi09", "e-",  "electron", "b-",    "electron" ] ]
    for y in yis :
        if( yi in y ) : return( y )
    raise Exception( "\nError from incidentParticleTags: invalid incident particle tag = %s" % `yi` )

def outgoingParticleTags( yo ) :
    """Calls incidentParticleTags but allows for yo to be integers in the range 10 to 19 inclusive."""

    yo_ = yo
    if( type( yo ) == type( 1 ) ) :
        if( ( yo > 9 ) and ( yo < 20 ) ) : yo -= 10
    return( [ yo_, ] + incidentParticleTags( yo ) )

def intYoCIS( fileName ) :
    """Returns the list of integers ( yo, C, I, S ) from fileName where fileName must be of the from 'yo##c##i###s###'."""

    if ( len( fileName ) != 15 ) : raise Exception( "\nError from intYoCIS: bad file name = %s" % `fileName` )
    return ( int( fileName[2:4] ), int( fileName[5:7] ), int( fileName[8:11] ), int( fileName[12:15] ) )

def intZASuffix( ZA_ii ) :
    """ZA_ii can be an integer (e.g., 92235) or a string of the from "zaZZZAAA_suffix" 
    (e.g., "za095242m"). Returns the tuple integer ZA and _suffix from ZA_i. For 
    examples, ZA_i = 92235 would return ( 92235, "" ) and ZA_i = "za095242m" would
    return ( 95242, "m" )."""

    try :
        ZA_i = int( ZA_ii )
    except :
        ZA_i = ZA_ii
    if ( type( ZA_i ) == type( 1 ) ) : return ( ZA_i, "" )
    ZA = None
    Suffix = ""
    if ( ( type( ZA_i ) == type( "" ) ) or ( type( ZA_i ) == type( u"" ) ) ) :
        if( ZA_i[:2] == "za" ) and ( len( ZA_i ) > 7 ) :
            ZA = int( ZA_i[2:8] )
            if ( len( ZA_i ) > 8 ) and ( ZA_i[8] != "/" ) : Suffix = ZA_i[8:]
    elif ( type( ZA_i ) == type( 1 ) ) :
        ZA = int( ZA_i )
    if( ZA == None ) : raise Exception( "\nError in intZA: invalid ZA = %s" % `ZA_i` )
    return ( ZA, Suffix ) 

def intZA( ZA_i ) :
    "Returns the ZA part of intZASuffix( ZA_i )."

    return intZASuffix( ZA_i )[0]

def strZASuffix( ZA, suffix = "" ) :
    """Converts ZA and suffix into 'za%6.6d%s' % ( ZA, suffix ). If ZA is already of the form 'za%6.6d%s', it is returned
    without any change."""

    if( suffix == 'm1' ) : suffix = 'm'
    if( type( ZA ) == type( 1 ) ) :
        strZA = "za%6.6d%s" % ( ZA, suffix )
    elif( type( ZA ) == type( "" ) ) :
        try :
            ZA = "za%6.6d" % int( ZA )
        except :
            pass
        strZA = "%s%s" % ( ZA, suffix )
    else :
        raise Exception( "ZA must be an integer or a string of the form 'zaZZZAAA'" )
    return strZA

def getmodestring_( I, callingRoutine ) :
    "For internal use only."

    if ( I in [ 0, 7, 9, 10, 11, 12, 13, 80, 89, 90, 91, 92, 941, 942, 951 ] ) :
        s = "2"
    elif ( I in [ 1, 21, 22, 81, 84, 952 ] ) :
        s = "3"
    elif ( I in [ 3, 4, 20 ] ) :
        s = "4"
    else :
        raise Exception( "\nError in %s: unsupported I-value = %d" % ( callingRoutine, I ) )
    return s

def getNumberOfColumns_( I, callingRoutine ) :
    "For internal use only."

    return int( getmodestring_( I, callingRoutine )[0] )

def processDataBase( database, oldFile, newFile, options = "", defines = [], bdflsFile = None, ndfgen = None, mcfgen = None, endep = None, extraZAs = [] ) :
    """For internal use only."""

    ndfgenFile = ""
    mcfgenFile = ""
    endepFile = ""
    if( bdflsFile != None ) :
        bdflsName = os.path.join( database, 'bdfls' )
        if( os.path.exists( bdflsName ) ) :
            if( not os.path.islink( bdflsName ) ) : raise Exception( '\nError in processDataBase: bdfls file already exists.' )
            os.system( 'rm -f ' + bdflsName )
        os.system( 'ln -s %s %s' % ( os.path.realpath( bdflsFile ), bdflsName ) )
    if( ndfgen != None ) :
        if( not os.path.exists( ndfgen ) ) : raise Exception( 'Error in processDataBase: ndfgen does not exists at %s.' % ndfgen )
        ndfgenFile = ' -ndfgen %s' % os.path.realpath( ndfgen )
    if( mcfgen != None ) :
        if( not os.path.exists( mcfgen ) ) : raise Exception( 'Error in processDataBase: mcfgen does not exists at %s.' % mcfgen )
        mcfgenFile = ' -mcfgen %s' % os.path.realpath( mcfgen )
    if( endep != None ) :
        if( not os.path.exists( endep ) ) : raise Exception( 'Error in processDataBase: endep does not exists at %s.' % endep )
        endepFile = ' -endep %s' % os.path.realpath( endep )
    ds = ""
    for d in defines : ds += 'export %s; ' % d
    if( "-no_u" not in options ) : options += " -u "
    if( "-no_i" not in options ) : options += " -i "
    if( len( extraZAs ) != 0 ) : 
        f = open( os.path.join( database, 'extraZAs' ), 'w' )
        for zas in extraZAs :
            ZA, suffix = intZASuffix( zas )
            f.write( "za%.6d%s\n" % ( ZA, suffix ) )
        f.close( )
    s = '$FUDGEPATH/Src/endlmod.com ' + options + ndfgenFile + mcfgenFile + endepFile + ' ' + database + ' ' + oldFile + ' ' + newFile
    if( ds != "" ) : print ds
    print s
    status = os.system( ds + s )
    if( status != 0 ) : raise Exception( "\nError in processDataBase: endlmod.com failed" )

def headerFunkyDouble2String( d ) :
    "For internal use only."

    nines = 9.999                       # Must be less than 15 nines otherwise log10 may not return as excepted.
    maxValue, abs_d = nines * 1e10, abs( d )
    if( d < 0 ) : maxValue /= 10
    if( 1e-2 < abs_d < maxValue ) :
        precisionMax = 8
        if( d > 0 ) : precisionMax += 1
        precision = max( 0, min( precisionMax, precisionMax - int( math.log10( abs_d ) ) ) )
        s = ( "%%11.%df" % precision ) % d
    else :
        precision = 6
        eOffset = 8
        if( d < 0 ) : precision -= 1
        if( ( abs_d < ( nines * 1e-100 ) ) or ( abs_d > ( nines * 1e99 ) ) ) :
            if( abs_d != 0 ) : eOffset -= 1
            precision -= 1
        s = ( "%%12.%de" % precision ) % d
        s = s[:eOffset] + s[eOffset+1:]
    return s

def headerString2FunkyDouble( s, i, callingRoutine, newMethod = True ) :
    "For internal use only."

    if( newMethod ) :
        fd = s[i:i+11]
        for i in xrange( 11 ) :
            if( fd[i] != ' ' ) : break
        if( ( i < 11 ) and ( fd[i] in '-+' ) ) : i += 1
        for i in xrange( i, 11 ) :
            if( fd[i] not in string.digits ) : break
        if( ( i < 11 ) and ( fd[i] == '.' ) ) :
            i += 1
            for i in xrange( i, 11 ) :
                if( fd[i] not in string.digits ) : break
        f = fd[:i]
        if( ( i < 11 ) and ( fd[i] in '-+' ) ) : f = f + 'e'
        for j in xrange( i, 11 ) :
            if( fd[j] == ' ' ) :
                f = f + '0'
            else :
                f = f + fd[j]
    else :
        f = s[i:i+11]
        try :
            d = float( f )
        except :
            if ( ( s[i+7] == "e" ) or ( s[i+7] == "E" ) ) :
                f = s[i:i+11]
            else :
                f = s[i:i+8] + "e" + s[i+8:i+11]
                if ( ( f[1] == " " ) or ( f[1] == "-" ) ) : f = f[1:8] + "0" + f[8:]
                if ( f[10] == " " ) : f = f[:10] + "0" + f[11:]
    try :
        d = float( f )
    except :
        raise Exception( "\nError in %s: invalid FunkyDouble\n<%s>" % ( callingRoutine, s[i:] ) )
    return( d )

def getFullPath( file, path = None ) :
    "For internal use only."

    if ( type( file ) != type( "" ) ) : raise Exception( "\nError in getFullPath: file name is not a string" )
    if ( len( file ) < 1 ) : raise Exception( "\nError in getFullPath: file name is an empty string" )
    if ( os.path.isabs( file ) ) :                  # Absolute path.
        f = file
    elif ( file[0] == "~" ) :                       # Users home directory.
        f = os.path.expanduser( file )
    elif ( os.path.exists( file ) ) :
        f = os.path.abspath( file )
    elif ( ( file[:2] != "./" ) and ( path != None ) ) :
        if ( type( path ) != type( "" ) ) : raise Exception( "\nError in getFullPath: path name is not a string" )
        if ( len( path ) < 1 ) : raise Exception( "\nError in getFullPath: path name is an empty string" )
        f = os.path.join( path, file )
    else :
        f = os.path.abspath( file )
    f = os.path.normpath( f )
    if ( not os.path.exists( f ) ) : return ( f, 0 )        # File or directory does not exist
    if (     os.path.isfile( f ) ) : return ( f, 1 )        # f is a file.
    if (     os.path.isdir( f ) )  : return ( f, 2 )        # f is a directory.
    raise Exception( "\nError in getFullPath: %s exists but is not a file or directory" )

def getAbsPath_withDefaultDBCheck( file ) :
    """This routine returns the absolute path for file using getFullPath and then os.path.abspath. If the file is in 
    fudgeDefaults.NUCLEAR_DATABASE_DIR, then the name will contain this even if fudgeDefaults.NUCLEAR_DATABASE_DIR
    is itself a link."""

    path, Status = getFullPath( file )
    path = os.path.realpath( path )
    DB_DIR = fudgeDefaults.NUCLEAR_DATABASE_DIR
    real_DB_DIR = os.path.realpath( DB_DIR )
    if( path[:len( real_DB_DIR )] == real_DB_DIR ) : path = os.path.join( DB_DIR, path[len( real_DB_DIR ) + 1:] )
    return( path, Status )

def appendyiDirIfNeeded( path, yi ) :
    """Appends the sub-directory yi## if it is not the last sub-directory of
    path.  ## is yi convert to a 2-digit integer."""

    Stryi = incidentParticleTags( yi )[1]
    if ( Stryi != path[-4:] ) : path = os.path.join( path, Stryi )
    return path

def validZADirectoryName_( name ) :
    "For internal use only."

    if ( ( re.match( r"za\d\d\d\d\d\d$", name ) ) or ( re.match( r"za\d\d\d\d\d\dm", name ) ) ) : return 1
    return 0

def check_mapMuEpToGrid( data ) :
    "For internal use only."

    E = None
    for d in data.data :
        if ( E != d[0] ) :
            if ( E != None ) :  print "E =", E, "  Sum =", s / 2.
            E = d[0]
            Ep2 = d[1]
            P2 = d[2]
            s = 0
        else :
            Ep1 = Ep2
            Ep2 = d[1]
            P1 = P2
            P2 = d[2]
            s += ( P2 + P1 ) * ( Ep2 - Ep1 )
    if ( E != None ) : print "E =", E, "  Sum =", s / 2.

def simplifyDataLine( l, n, callingRoutine, comment ) :
    "For internal use only."

    ls = string.split( l, comment )[0]
    ls = string.replace( ls, ",", " " )
    ls = string.replace( ls, "\t", " " )
    if ( ls == "" ) : return None
    v = string.split( ls )
    if( ( n != None ) and ( len( v ) != n ) ) :
        raise Exception( "\nError in %s: line has %d values, expected %d\n%s" % ( callingRoutine, len( v ), n, l ) )
    return v

def getFormat( self ) :
    """Checks if self has a getFormat method. If it does not, then None is returned. If it does then: 
    1) if getFormat returns a string or None then it is returned else getFormat must return an 
    integer, n1. If n1 is 12 (default for legacy ENDL) then None is returned else returns string 
    '%n1.n2e' where n2 = n1 - 7."""

    try :
        format = self.getFormat( )
        if( ( format == None ) or ( type( format ) == type( '' ) ) ) : return( format )
        if( format == 12 ) :
            format = None
        else :
            format = '%%%d.%de' % ( format, format - 7 )
    except :
        format = None
    return( format )

def read1dDataFile( fn, comment = "#", i0 = 0 ) :
    """Opens file with name fn, reads in the data, closes the file and then
    calls translate1dStringData with the data.  Returns the results of translate1dStringData."""

    f = open( fn, "r" )
    ls = f.readlines( )
    f.close( )
    return translate1dStringData( ls, comment = comment, i0 = i0 )

def translate1dStringData( ls, comment = "#", i0 = 0 ) :
    """
Translates 1d data from the list of strings ls and returns ::
    list[ number ].
This list can be used as input when creating endl1dmath instances. All strings starting 
with comment are ignored and all data in a string after comment are ignored.  
After comments are removed, the strings can only contain one numbers per line 
else a raise is issued.  The x-component is column i0."""

    if ( i0 != 0 ) :
        raise Exception( "\nError in endlmisc.translate1dStringData: Invalid i0 = %s value" % ( `i0` ) )
    d1 = []
    n = None
    for l in ls :
        sl = simplifyDataLine( l, n, "translate1dStringData", comment )
        if( sl != None ) :
            if( n == None ) : n = len( sl )
            d1.append( float( sl[i0] ) )
    return d1

def read2dDataFile( fn, comment = "#", i0 = 0, i1 = 1 ) :
    """Opens file with name fn, reads in the data, closes the file and then
    calls translate2dStringData with the data.  Returns the results of translate2dStringData."""

    f = open( fn, "r" )
    ls = f.readlines( )
    f.close( )
    return translate2dStringData( ls, comment = comment, i0 = i0, i1 = i1 )

def translate2dStringData( ls, comment = "#", i0 = 0, i1 = 1 ) :
    """
Translates 2d data from the list of strings ls and returns ::
    list[ number, number ].
This list can be used as input when creating endl2dmath instances. All strings starting 
with comment are ignored and all data in a string after comment are ignored.  
After comments are removed, the strings can only contain two numbers per line 
else a raise is issued.  The x-component is column i0 of the data and the
y-component is column i1 (i0 and i1 can only be 0 or 1)."""

    i = [ i0, i1 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) :
        raise Exception( "\nError in endlmisc.translate2dStringData: Invalid i0 = %s and/or i1 = %s value(s)" % ( `i0`, `i1` ) )
    ConvertDs = 0
    if ( len( ls ) > 0 ) and ( len( ls[0] ) > 25 ) :
        if ( ls[0][10] == "D" ) and ( ls[0][22] == "D" ) : ConvertDs = 1
    d2 = []
    n = None
    for l in ls :
        if ( ConvertDs ) : l = re.sub( r'[D]', 'E', l )
        sl = simplifyDataLine( l, 2, "translate2dStringData", comment )
        if( sl != None ) : 
            if( n == None ) : n = len( sl )
            d2.append( [ float( sl[i0] ), float( sl[i1] ) ] )
    return d2

def read3dDataFile( fn, comment = "#", i0 = 0, i1 = 1, i2 = 2 ) :
    """Opens file with name fn, reads in the data, closes the file and then
    calls translate3dStringData with the data.  Returns the results of translate3dStringData."""

    f = open( fn, "r" )
    ls = f.readlines( )
    f.close( )
    return translate3dStringData( ls, comment = comment, i0 = i0, i1 = i1, i2 = i2 )

def translate3dData( data, i0 = 0, i1 = 1, i2 = 2 ) :
    """
Translates 3d data from the list of (x,y,z) data and returns ::
    list[ number, list[ number, number ] ].
The returned list can be used as input when creating endl3dmath instances.
The x-component is column i0 of the data, the y-component is column i1 and
the z-component is column i2 (i0, i1 and i2 can only be 0, 1 or 2)."""

    i = [ i0, i1, i2 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) :
        raise Exception( "\nError in endlmisc.translate3dData: Invalid i0 = %s, i1 = %s and/or i2 = %s value(s)" % ( `i0`, `i1`, `i2` ) )
    x0 = None
    d3 = []
    for d in data :
        x = d[i0]
        if ( x != x0 ) :
            if ( x0 != None ) : d3.append( [ x0, d2 ] )
            d2 = []
            x0 = x
        d2.append( [ d[i1], d[i2] ] )
    if ( x0 != None ) : d3.append( [ x, d2 ] )
    return d3

def translate3dStringData( ls, comment = "#", i0 = 0, i1 = 1, i2 = 2 ) :
    """
Translates 3d data from the list of strings ls and returns ::
    list[ number, list[ number, number ] ]. 
The returned list can be used as input when creating endl3dmath instances. All strings
starting with comment are ignored and all data in a string after comment are ignored.
After comments are removed, the strings can only contain three numbers per line
else a raise is issued.  The x-component is column i0 of the data, the y-component
is column i1 and the z-component is column i2 (i0, i1 and i2 can only be 0, 1 or 2)."""

    i = [ i0, i1, i2 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) :
        raise Exception( "\nError in endlmisc.translate3dStringData: Invalid i0 = %s, i1 = %s and/or i2 = %s value(s)" % ( `i0`, `i1`, `i2` ) )
    x0 = None
    d3 = []
    n = None
    for l in ls :
        sl = simplifyDataLine( l, 3, "translate3dStringData", comment )
        if( sl != None ) :
            if( n == None ) : n = len( sl )
            x = float( sl[i0] )
            if ( x != x0 ) :
                if ( x0 != None ) : d3.append( [ x0, d2 ] )
                d2 = []
                x0 = x
            d2.append( [ float( sl[i1] ), float( sl[i2] ) ] )
    if ( x0 != None ) : d3.append( [ x, d2 ] )
    return d3

def print3dData( d, i0 = 0, i1 = 1, i2 = 2, fmt0 = None, fmt1 = None, fmt2 = None ) :
    """Calls string3dData and print its results."""

    s = string3dData( d, i0 = i0, i1 = i1, i2 = i2, fmt0 = fmt0, fmt1 = fmt1, fmt2 = fmt2 )
    for i in s : print i

def string3dData( d, i0 = 0, i1 = 1, i2 = 2, fmt0 = None, fmt1 = None, fmt2 = None ) :
    """Returns a list of strings for data d.  D must be list[ number, list[ number, number ] ]."""

    from fudge.legacy.endl import endl3dmathmisc
    i = [ i0, i1, i2 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) :
        raise Exception( "\nError in endlmisc.string3dData: Invalid i0 = %s, i1 = %s and/or i2 = %s value(s)" % ( `i0`, `i1`, `i2` ) )
    if( fmt0 == None ) : fmt0 = endl3dmathmisc.endl3d_repr_xFormat
    if( fmt1 == None ) : fmt1 = endl3dmathmisc.endl3d_repr_yFormat
    if( fmt2 == None ) : fmt2 = endl3dmathmisc.endl3d_repr_zFormat
    fmt = [ fmt0, fmt1, fmt2 ]
    s = []
    sa = [ "", "", "" ]
    for d0 in d :
        sa[i0] = fmt[i0] % d0[0]
        for d1 in d0[1] :
            sa[i1] = fmt[i1] % d1[0]
            sa[i2] = fmt[i2] % d1[1]
            s.append( "%s %s %s" % ( sa[0], sa[1], sa[2] ) )
    return s

def read4dDataFile( fn, comment = "#", i0 = 0, i1 = 1, i2 = 2, i3 = 3 ) :
    """Opens file with name fn, reads in the data, closes the file and then calls
    translate4dStringData with the data.  Returns the results of translate4dStringData."""

    f = open( fn, "r" )
    ls = f.readlines( )
    f.close( )
    return translate4dStringData( ls, comment = comment, i0 = i0, i1 = i1, i2 = i2, i3 = i3 )

def translate4dData( data, comment = "#", i0 = 0, i1 = 1, i2 = 2, i3 = 3 ) :
    """
Translates 4d data from the list of (x,y,z,zz) data and returns ::
    list[ number, list[ number, list[ number, number ] ] ].
The returned list that can be used as input when creating endl4dmath instances.
The x-component is column i0 of the data, the y-component is column i1, the
z-component is column i2 and the zz-component is column i3
(i0, i1, i2 and i3 can only be 0, 1, 2 or 3)."""

    i = [ i0, i1, i2, i3 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) or ( i[3] != 3 ) :
        s = "Invalid i0 = %s, i1 = %s, i2 = %s and/or i3 = %s value(s)" % ( `i0`, `i1`, `i2`, `i3` )
        raise Exception( "\nError in endlmisc.translate4dData: %s" % s )
    x0 = None
    y0 = None
    d4 = []
    for d in data :
        x = d[i0]
        y = d[i1]
        if ( x != x0 ) :
            if ( y0 != None ) : d3.append( [ y0, d2 ] )
            if ( x0 != None ) : d4.append( [ x0, d3 ] )
            d2 = []
            d3 = []
            x0 = x
            y0 = y
        elif ( y != y0 ) :
            d3.append( [ y0, d2 ] )
            d2 = []
            y0 = y
        d2.append( [ d[i2], d[i3] ] )
    if ( y0 != None ) : d3.append( [ y0, d2 ] )
    if ( x0 != None ) : d4.append( [ x0, d3 ] )
    return d4

def translate4dStringData( ls, comment = "#", i0 = 0, i1 = 1, i2 = 2, i3 = 3 ) :
    """
Translates 4d data from the list of strings ls and returns ::
    list[ number, list[ number, list[ number, number ] ] ].

This list that can be used as input when creating endl4dmath instances. All strings
starting with comment are ignored and all data in a string after comment are ignored.
After comments are removed, the strings can only contain four numbers per line
else a raise is issued.  The x-component is column i0 of the data, the y-component
is column i1, the z-component is column i2 and the zz-component is column i3 
(i0, i1, i2 and i3 can only be 0, 1, 2 or 3)."""

    i = [ i0, i1, i2, i3 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) or ( i[3] != 3 ) :
        s = "Invalid i0 = %s, i1 = %s, i2 = %s and/or i3 = %s value(s)" % ( `i0`, `i1`, `i2`, `i3` )
        raise Exception( "\nError in endlmisc.translate4dStringData: %s" % s )
    x0 = None
    y0 = None
    d4 = []
    n = None
    for l in ls :
        sl = simplifyDataLine( l, 4, "translate4dStringData", comment )
        if( sl != None ) :
            if( n == None ) : n = len( sl )
            x = float( sl[i0] )
            y = float( sl[i1] )
            if ( x != x0 ) :
                if ( y0 != None ) : d3.append( [ y0, d2 ] )
                if ( x0 != None ) : d4.append( [ x0, d3 ] )
                d2 = []
                d3 = []
                x0 = x
                y0 = y
            elif ( y != y0 ) :
                d3.append( [ y0, d2 ] )
                d2 = []
                y0 = y
            d2.append( [ float( sl[i2] ), float( sl[i3] ) ] )
    if ( y0 != None ) : d3.append( [ y0, d2 ] )
    if ( x0 != None ) : d4.append( [ x0, d3 ] )
    return d4

def print4dData( d, i0 = 0, i1 = 1, i2 = 2, i3 = 3, fmt0 = None, fmt1 = None, fmt2 = None, fmt3 = None ) :
    """Calls string4dData and print its results."""

    s = string4dData( d, i0 = i0, i1 = i1, i2 = i2, i3 = i3, fmt0 = fmt0, fmt1 = fmt1, fmt2 = fmt2, fmt3 = fmt3 )
    for i in s : print i

def string4dData( d, i0 = 0, i1 = 1, i2 = 2, i3 = 3, fmt0 = None, fmt1 = None, fmt2 = None, fmt3 = None ) :
    """Returns a list of strings for data d.  D must be list[ number, list[ number, list[ number, number ] ] ]."""

    from fudge.legacy.endl import endl4dmathmisc
    i = [ i0, i1, i2, i3 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) or ( i[3] != 3 ) :
        raise Exception( "\nError in endlmisc.string4dData: %s" % "Invalid i0 = %s, i1 = %s, i2 = %s, and/or i3 = %s" % (`i0`,`i1`,`i2`,`i3`) )
    if( fmt0 == None ) : fmt0 = endl4dmathmisc.endl4d_repr_tFormat
    if( fmt1 == None ) : fmt1 = endl4dmathmisc.endl4d_repr_xFormat
    if( fmt2 == None ) : fmt2 = endl4dmathmisc.endl4d_repr_yFormat
    if( fmt3 == None ) : fmt3 = endl4dmathmisc.endl4d_repr_zFormat
    fmt = [ fmt0, fmt1, fmt2, fmt3 ]
    s = []
    sa = [ "", "", "", "" ]
    for d0 in d :
        sa[i0] = fmt[i0] % d0[0]
        for d1 in d0[1] :
            sa[i1] = fmt[i1] % d1[0]
            for d2 in d1[1] :
                sa[i2] = fmt[i2] % d2[0]
                sa[i3] = fmt[i3] % d2[1]
                s.append( ' '.join( sa ) )
    return s

def stringWithPrefixSuffix( list, Prefix = "", Suffix = "" ) :
    "For internal use only."

    PS = "\n"
    if( len( list ) > 0 ) :
        list[0] = Prefix + list[0]
        list[-1] = list[-1] + Suffix + "\n"
        PS = Suffix + "\n" + Prefix
    list = PS.join( list )
    return list

class endlCheckerObject :

    sortOrder = [ 'database', 'yi', 'ZA', 'suffix', 'C', 'S', 'X1', 'yo', 'Q', 'I', 'X2', 'X3', 'X4', 'message' ]
    
    def __init__( self, data = None, database = None, yi = None, ZA = None, suffix = None, yo = None, C = None, I = None, S = None, X1 = None, \
        X2 = None, X3 = None, X4 = None, Q = None, message = '' ) :

        self.database = database
        self.yi = yi
        self.ZA = ZA
        self.suffix = suffix
        self.yo = yo
        self.C = C
        self.I = I
        self.S = S
        self.X1 = X1
        self.X2 = X2
        self.X3 = X3
        self.X4 = X4
        self.Q = Q
        self.message = message
        self.indentation = '  '
        if( data != None ) :
            self.yi = data.yi
            self.ZA = data.ZA
            self.yo = data.yo
            self.C = data.C
            self.I = data.I
            self.S = data.S
            self.X1 = data.X1
            self.X2 = data.X2
            self.X3 = data.X3
            self.X4 = data.X4
            self.Q = data.Q

    def __cmp__( self, other ) :

        def diffWithNone( v1, v2 ) :

            if( v1 == v2 ) : return( 0 )
            if( v1  > v2 ) : return( 1 )
            return( -1 )

        if( type( self ) != type( other ) ) : raise Exception( 'Error in endlmisc.endlCheckerObject.__cmp__: other not a endlCheckerObject class object' )
        for attr in endlCheckerObject.sortOrder :
            if( getattr( self, attr ) != getattr( other, attr ) ) : return( diffWithNone( getattr( self, attr ), getattr( other, attr ) ) )
        return( 0 )
        
    def __hash__( self ):
        h = ( hash( self.X2 ) + hash( self.X3 ) + hash( self.X4 ) + hash( self.Q ) ) \
            + hash( self.C ) * 100 \
            + ( hash( self.S ) + hash( self.X1 ) ) * 1000 \
            + hash( self.yo ) * 10000 \
            + hash( self.I ) * 100000 \
            + ( hash( self.ZA ) + hash( self.suffix ) ) * 1000000 \
            + ( hash( self.yi ) ) * 10000000 \
            + ( hash( self.database ) ) * 100000000
        if type( self.message ) == str: h += hash( self.message )
        elif type( self.message ) == list: h += sum( map( hash, self.message ) )
        return h

    def __repr__( self ) :

        return( '\n'.join( self.toListOfStrings( indentation = self.indentation, subIndentation = self.indentation ) ) )

    def toString( self, selections = [], returnNoneOnNoMatches = True ) :

        slos = self.toSelectedListOfStrings( indentation = self.indentation, subIndentation = self.indentation, selections = selections )
        if( ( len( slos ) == 1 ) and returnNoneOnNoMatches ) : return( None )
        return( '\n'.join( slos ) )

    def toSelectedListOfStrings( self, indentation = None, subIndentation = None, selections = [] ) :

        los = self.toListOfStrings( indentation = self.indentation, subIndentation = self.indentation )
        if( len( selections ) > 0 ) :
            slos = los[:1]
            for lo in los[1:] :
                for selection in selections :
                    if( re.search( selection, lo ) ) :
                        slos.append( lo )
                        break
            los = slos
        return( los )

    def toListOfStrings( self, indentation = None, subIndentation = None ) :

        def toInt( v, n ) :

            abs_n = abs( n )
            try :
                s = '%d' % v
            except :
                s = abs_n * '?'
            sn = len( s )
            filler = '0'
            if( n < 0 ) : filler = ' '
            if( sn < abs_n ) : s = ( abs_n - sn ) * filler + s
            return( s )

        def toFloat( v ) :

            try :
                s = '%e' % v
            except :
                s = '????????????'
            return( s )

        def toStr( s ) :
            
            if s is None : s = ''
            s = str( s )
            return( s )

        if( indentation == None ) : indentation = self.indentation
        if( subIndentation == None ) : subIndentation = self.indentation
        s = []
        if( self.message != '' ) :
            yi = toInt( self.yi, -2 )
            ZA = toInt( self.ZA, -6 )
            suffix = toStr( self.suffix )
            yo = toInt( self.yo, 2 )
            C = toInt( self.C, 2 )
            I = toInt( self.I, 3 )
            S = toInt( self.S, 3 )
            X1 = toFloat( self.X1 )
            X2 = toFloat( self.X2 )
            X3 = toFloat( self.X3 )
            X4 = toFloat( self.X4 )
            Q = toFloat( self.Q )
            fileName = 'yo%sc%si%ss%s' % ( yo, C, I, S )
            s = [ indentation + 'Error for yi = %s, ZA = %s%s, %s, X1 = %s, X2 = %s, X3 = %s, X4 = %s, Q = %s' % ( yi, ZA, suffix, fileName, X1, X2, X3, X4, Q ) ]
            s += checkMessagesToString( self.message, indentation = indentation + subIndentation )
        return( s )

def checkMessagesToString( message, indentation = '  ', subIndentation = '  ' ) :

    s = []
    if( type( message ) == type( '' ) ) :
        s += [ indentation + message ]
    else :
        for m in message :
            if( type( m ) == type( '' ) ) :
                s += [ indentation + m  ]
            elif( ( type( m ) == type( [] ) ) or ( type( m ) == type( [] ) ) ) :
                if( len( m ) != 2 ) : raise Exception( 'Error in endlmisc.checkMessagesToString: len( m ) = %d != 2 )' % len( m ) )
                s += [ indentation + m[0] ]
                s += checkMessagesToString( m[1], indentation = indentation + subIndentation  )
            else :
                raise Exception( 'Error in endlmisc.checkMessagesToString: invalid message type = %s' % type( m ) )
    return( s )

def checkCulling( errs, cullStrings ) :
    """
    Removes all messages containing sub-strings listed in cullStrings. cullStrings can be either a string or a
    list of strings. If as list of strings, each string must be a sub-string in a message for the message to
    be culled.
    """

    def checkCullingMatch( message, cullStrings ) :

        found = True
        for cullString in cullStrings : found = found and ( cullString in message )
        return( found )

    def checkCulling2( message, cullStrings, level = 0 ) :

        if( isinstance( message, list ) ) :
            messages = []
            for msg in message :
                msg1 = checkCulling2( msg, cullStrings, level + 1 )
                if( msg1 is not None ) : messages.append( msg1 )
            if( len( messages ) < 2 ) : messages = None
            return( messages )
        else :
            if( checkCullingMatch( message, cullStrings ) ) : return( None )
            return( message )

    if( isinstance( cullStrings, str ) ) : cullStrings = [ cullStrings ]
    errs2 = []
    for err in errs :
        messages = []
        if( isinstance( err.message, str ) ) :
            if( not( checkCullingMatch( err.message, cullStrings ) ) ) : errs2.append( err )
        else :
            for message in err.message :
                message = checkCulling2( message, cullStrings )
                if( message is not None ) :
                    messages.append( message )
            if( len( messages ) > 0 ) :
                err.message = messages
                errs2.append( err )
    return( errs2 )

def printWarning( s ) :

    sys.stderr.write( s )
    if( s[-1:] != '\n' ) : sys.stderr.write( '\n' )

def endlToFudgeInterpolation( interpolation ) :
    """This function converts an endl interpolation value (0 or 2 is lin-lin, 3 is log-lin, 4 is lin-log and
    5 is log-log) into a fudge interpolation value (0 is lin-lin, 1 is log-lin, 2 is lin-log and 3 is log-log)."""

    if( ( interpolation < 0 ) or ( interpolation > 5 ) or ( interpolation == 1 ) ) : raise Exception( "Invalid ENDL interpolation value = %d" % interpolation )
    return( ( 0, None, 0, 1, 2, 3 )[interpolation] )

def fudgeToEndlInterpolation( interpolation ) :
    """This function converts a fudge interpolation value into an end interpolation value (see endlToFudgeInterpolation
    for valid interpolation vales."""

    if( ( interpolation < 0 ) or ( interpolation > 3 ) ) : raise Exception( "Invalid FUDGE interpolation value = %d" % interpolation )
    return( ( 0, 3, 4, 5 )[interpolation] )
