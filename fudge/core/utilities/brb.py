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
This module contains a set of functions that are hopefully useful when using python.
"""
import sys
import os
import re
import types
import glob

def uniquify(seq):
    ''' 
    Fast implimentation of a function to strip out non-unique entries in a Python list, but preserving list order.
    Usage:
        >>> print uniquify( [1,4,2,4,7,2,1,1,1,'s','a',0.1] ) 
            [1, 4, 2, 7, 's', 'a', 0.10000000000000001]
    '''
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)] 

def banner( s, justification = 'c' ):
    '''A flower-box ascii banner generator.  The input string can have line feeds.  Each line will get centered in the box.'''

    if( justification not in [ 'l', 'c', 'r' ] ) : raise Exception( 'Invalid justification = "%s", must be "l", "c", "r"' % justification )
    ss = s.split( '\n')
    m = 0
    for l in ss: m = max( m, len( l ) )
    ans = '+-' + m*'-' + '-+\n'
    if(   justification == 'l' ) :
        for l in ss: ans += '| ' + l.ljust( m ) + ' |\n'
    elif( justification == 'c' ) :
        for l in ss: ans += '| ' + l.center( m ) + ' |\n'
    else :
        for l in ss: ans += '| ' + l.rjust( m ) + ' |\n'
    ans += '+-' + m*'-' + '-+'
    return ans

# Small banner, with wings
def winged_banner( x, wingsize=10 ):
    return wingsize*'*'+' '+x.replace( '\n', '; ' )+' '+wingsize*'*'

def Pause( Prompt = "Enter <RET> to continue : ", ExtraStr = "" ) :
    "Prompts and waits for user to enter a line. The line is returned. ExtraStr is added to the prompt string."

    if( ExtraStr != "" ) : Prompt = Prompt + "(%s)" % ExtraStr
    try :
        s = raw_input( Prompt )
    except EOFError :
        print
        s = None
    except :
        raise
    return( s )

def tlist( l, w = 5, sep = None, rightJustify = 0 ) :
    "Calls tylist with arguments passed."

    tylist( l, w, sep = sep, rightJustify = rightJustify )

def txlist( l, w = 5, sep = None, rightJustify = 0 ) :
    """Prints the list given by the first argument in a tabled format. The purpose of this 
    function is to print a python list in a more readable format. The number of items from list
    printed per line is set by w.  The argument sep, if not None, is printed as a separator
    between each item on a line. The items are left justified unless the argument rightJustify
    is true. This function prints each w consecutive items horizontally. Also see tylist."""

    if ( w < 0 ) : w = 1
    s = tlistMaxLen( l, rightJustify )
    if s == "" :
        ll = []
        for i in l : ll.append( `i` )
        s = tlistMaxLen( ll, rightJustify )
    else :
        ll = l
    Counter = 0
    for i in ll :
        print s % i,
        Counter = Counter + 1
        if( Counter == w ) :
            Counter = 0
            print ""
        else :
            if( sep is not None ) : print sep,
    if( Counter > 0 ) : print ""

def tylist( l, w = 5, sep = None, rightJustify = 0 ) :
    """Prints the list given by the first argument in a tabled format. The purpose of this 
    function is to print a python list in a more readable format. The number of items from list
    printed per line is set by w.  The argument sep, if not None, is printed as a separator
    between each item on a line. The items are left justified unless the argument rightJustify
    is true. This function prints each consecutive item vertically while maintaining w columns. 
    Also see txlist."""

    if ( w < 0 ) : w = 1
    s = tlistMaxLen( l, rightJustify )
    if s == "" :
        ll = []
        for i in l : ll.append( `i` )
        s = tlistMaxLen( ll, rightJustify )
    else :
        ll = l
    r = ( len( ll ) + w - 1 ) / w
    for ir in range( r ) :
        iw = 0
        for i in range( ir, len( ll ), r ) :
            print s % ll[ i ],
            iw += 1
            if( ( sep is not None ) and ( iw < w ) ) : print sep,
        print ""

def tlistMaxLen( l, rightJustify ) :
    """For internal use only."""

    MaxLen = 2
    for i in l : 
        if ( type( i ) == type( 1 ) ) : MaxLen = max( MaxLen, 8 )
        elif ( type( i ) == type( "" ) ) : MaxLen = max( MaxLen, len( i ) + 1 )
        else : return ""
    justify = "-"
    if( rightJustify ) : justify = ""
    return "%%%s%ds" % ( justify, MaxLen )

def tdir( a = None, w = 5, pattern = None, wpattern = None ) :
    """This function is an attempt at improving the output of python's dir function. It uses
    tlist to print the output in a more readable format. The argument w is passed to tlist.
    The argument pattern, if not None, is used to restrict the items displayed to those that match it.
    Pattern can be any pattern understood by the re module (e.g., to display all items
    containing the two consecutive letters "2d" set pattern to ".*2d.*"). If pattern is None
    and wpattern is NOT None then pattern is set to ".*" + wpattern + ".*"."""

    if( ( pattern is None ) and ( wpattern is not None ) ) : pattern = ".*" + wpattern + ".*"
    if( ( ( type( a ) == type( [] ) ) or ( type( a ) == type( () ) ) ) and ( len( a ) > 0 ) ) :
        l = a
    elif ( type( a ) == type( {} ) ) :
        l = a.keys( )
    else :
        l = dir( a )
    if( pattern is not None ) :
        m = []
        p = re.compile( pattern )
        for i in l :
            if( p.match( i ) ) : m.append( i )
        l = m
    tlist( l, w )

def doc( a ) :
    """This function simply executes "print a.__doc__" where "a" is the argument."""

    print a.__doc__

def docmethods( o, Full = 1 ) :
    """This function attempts to print all methods and their documentation contained in the first argument."""

    def docmethods2( o, typeslist ) :
        for i in dir( o ) :
            a = getattr( o, i )
            if ( type( a ) in typeslist ) :
                print "\n------------------ %-30s ------------------" % i
                print a.__doc__

    typeslist = [ types.FunctionType, types.MethodType, types.UnboundMethodType ]
    if ( Full > 0 ) : typeslist.append( types.ModuleType )
    if ( Full > 5 ) :
        typeslist.append( types.BuiltinFunctionType )
        typeslist.append( types.BuiltinMethodType )
    if   ( type( o ) == types.InstanceType ) : docmethods2( o.__class__, typeslist )    # Instance of a class.
    elif ( type( o ) == types.ClassType ) : docmethods2( o, typeslist )                 # Class
    elif ( type( o ) == types.ModuleType ) : docmethods2( o, typeslist )                # Function

def py( depth = 1, start = 0 ) :
    """This routine prints all python scripts in sys.path[start:start+depth] in a nice format."""

    Counter = 0
    for d in sys.path :
        if ( not d == "" ) :
            if( Counter >= start ) :
                print "\n", d
                p = os.path.join( d, "*.py" )
                fs = glob.glob( p )
                l = []
                for f in fs : l.append( os.path.split( f )[1] )
                txlist( l )
            Counter = Counter + 1
            if( Counter == ( start + depth ) ) : return

def objectoutline( o, MaxLevel = 1, Full = 1 ) :
    """This routine attempts to print the objects - methods, members and possible other things - that are a part
    of the first argument. It will recursively traverse each object found down to level MaxLevel."""

    objectoutline2( o, MaxLevel, 0, "", Full )

def objectoutline2( o, MaxLevel, level, name, Full ) :
    "For internal use only."

    if( type( name ) != type( "" ) ) : name = str( name )
    if ( level == ( MaxLevel + 1 ) ) : return
    DoNextLevel = not ( level == MaxLevel )
    nextlevel = level + 1
    sm = " " * ( 4 * level )
    s = sm
    isInstance = ( type( o ) == types.InstanceType )
    if( hasattr( o, '__class__' ) ) : isInstance |= ( type( o ) == o.__class__ )    # New style class instance.
    
    try :
        if ( name != "" ) : s = s + name + ": "
    except :
        print '<%s> <%s>, %s, %s' % ( `s`, `name`, `type( s )`, `type( name )` )
        raise

    if ( type( o ) == types.ModuleType ) :                    # Function
        print "%smethod %s" % ( s, o.__name__ )
        for i in dir( o ) : objectoutline2( getattr( o, i ), MaxLevel, nextlevel, i, Full )
    elif ( type( o ) == type( [] ) ) :                          # List of objects
        print "%slist of len %d" % ( s, len( o ) )
        if ( ( Full ) and ( CheckListForSameTypeOfObjects( o, [] ) ) ) :
            if ( DoNextLevel ) :
                for i in o : objectoutline2( i, MaxLevel, nextlevel, "", Full )
        else :
            if ( o != [] ) : objectoutline2( o[0], MaxLevel, nextlevel, "", Full )
    elif ( type( o ) == type( () ) ) :                          # Tuple of objects
        print "%stuple of len %d" % ( s, len( o ) )
        if ( Full ) :
            if ( DoNextLevel ) :
                for i in o : objectoutline2( i, MaxLevel, nextlevel, "", Full )
        else :
            if ( o != () ) : objectoutline2( o[0], MaxLevel, nextlevel, "", Full )
    elif ( type( o ) == type( {} ) ) :                          # Dictionary
        print "%sdictionary of %d items" % ( s, len( dir( o ) ) )
        for i in o.keys( ) : objectoutline2( o[i], MaxLevel, nextlevel, i, Full )
    elif( ( type( o ) == types.ClassType ) or ( type( o ) == types.TypeType ) ) :                     # Class old and new style
        print "%sclass %s" % ( s, o.__name__ )
        if ( DoNextLevel ) :
            if ( len( o.__bases__ ) > 0 ) :
                print "%s  -- base(s) --" % sm
                for i in o.__bases__ : objectoutline2( i, MaxLevel, nextlevel, i.__name__, Full )
            print "%s  -- methods --" % sm
            for i in dir( o ) : objectoutline2( getattr( o, i ), MaxLevel, nextlevel, i, Full )
        else :
            if ( len( o.__bases__ ) > 0 ) :
                print "%s  -- base(s) --" % sm
                for i in o.__bases__ : objectoutline2ClassesOnly( i, nextlevel, i.__name__ )
    elif ( type( o ) == type( 1 ) ) or ( type( o ) == type( 1. ) ) or ( type( o ) == type( "" ) ) or ( type( o ) == type( u"" ) ) or ( type( o ) == type( True ) ) :
        print "%s%s = %s" % ( s, `type( o )`, `o` )
    elif ( type( o ) == types.MethodType ) :
        print "%s %s" % ( s, o.__class__.__name__ )
    elif( isInstance ) :                  # Instance of a class.
        print "%sinstance of class %s" % ( s, o.__class__.__name__ )
        if ( DoNextLevel ) :
#            if ( len( o.__class__.__bases__ ) > 0 ) :
#                print "%s  -- bases --" % sm
#                for i in o.__class__.__bases__ : objectoutline2( i, MaxLevel, nextlevel, i.__name__, Full )
            print "%s  -- member --" % sm
            for i in dir( o ) : 
                if ( type( getattr( o, i ) ) != types.MethodType ) : objectoutline2( getattr( o, i ), MaxLevel, nextlevel, i, Full )
            print "%s  -- methods --" % sm
            for i in dir( o ) :
                if ( type( getattr( o, i ) ) == types.MethodType ) : objectoutline2( getattr( o, i ), MaxLevel, nextlevel, i, Full )
    else :                                                      # Integer, Float, etc...
        print "%s%s" % ( s, `type( o )` )

def objectoutline2ClassesOnly( o, level, name ) :
    "For internal use only!"

    sm = " " * ( 4 * level )
    if( ( type( o ) == types.ClassType ) or ( type( o ) == types.TypeType ) ) :                     # Class old and new style
        print "%sclass %s" % ( sm, o.__name__ )
        if ( len( o.__bases__ ) > 0 ) :
            if( len( o.__bases__[0].__bases__ ) == 0 ) : return
            print "%s  -- base(s) --" % sm
            for subo in o.__bases__ : objectoutline2ClassesOnly( subo, level + 1, subo.__name__ )

def CheckListForSameTypeOfObjects( l, e ) :
    """For internal use only."""

    if( l == e ) : return 0
    t = type( l[0] )
    for i in l :
        if ( t != type( i ) ) : return 1
    return 0

def objectvalues( o ) :
    """Prints a list of all objects return by dir( o ) and their data (or object type) where "o" 
    is the first argument."""

    if ( type( o ) == types.InstanceType ) :                  # Instance of a class.
        print "class %s" % o.__class__.__name__
        print "  -- member --"
        for i in dir( o ) : objectvalues2( getattr( o, i ), i )
    else :
        objectvalues2( o, "" )

def objectvalues2( o, name ) :
    """For internal use only."""
    if   ( type( o ) == types.InstanceType ) :                  # Instance of a class.
        return
    elif ( type( o ) == types.ModuleType ) :                    # Function
        return
    if ( type( o ) == type( [] ) ) :                            # List of objects
        print "  %s: list of len %d" % ( name, len( o ) )
    elif ( type( o ) == type( () ) ) :                          # Tuple of objects
        print "  %s: tuple of len %d" % ( name, len( o ) )
    elif ( type( o ) == type( {} ) ) :                          # Dictionary
        print "  %s dictionary of %d items" % ( name, len( dir( o ) ) )
        for i in o.keys( ) : print "    %s = %s" % ( i, o[i] )
    elif ( type( o ) == types.ClassType ) :                     # Class
        return
    else :                                                      # Integer, Float, etc...
        print "  %s: %s = %s" % ( name, `type( o )`, `o` )

def getType( o ) :

    try :
        return( "%s %s" % ( str( o.im_class ), str( o.im_func ) ) )
    except :
        try :
            return( str( o.__class__ ) )
        except :
            return( type( o ) )

def limitObjectToString( object ) :

    if( isinstance( object, str ) ) :
        s1 = object
    else :
        s1 = repr( object )
    if( len( s1 ) > 64 ) : s1 = s1[:61] + '...'
    return( s1 )
