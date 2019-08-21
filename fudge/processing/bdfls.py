# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains classes for dealing with a bdfls file and its contents.
"""

import os
import string
from fudge.core import fudgemisc
from fudge.core.utilities import fudgeFileMisc, fudgeZA, subprocessing
from fudge.vis.gnuplot import plotbase

default_bdfls = None

# set default bdfls file path
_BDFLS_FILE_ = "/usr/gapps/data/nuclear/bdfls"
bdfls_EOD =    "                                                                        1\n"
bdfls_mh_EOD = "    999999                                                              1\n"

def targetToZA( za ) :
    """This routine converts an integer or string into a ZA (e.g., "FissionProduct_ENDL99120" to 99120)."""

    try :
        iza = int( za )
    except :
        if( type( za ) == ( type( '' ) ) ) :
            if( za[:2] == 'za' ) :
                iza = int( za[2:8] )
            elif( za[:19] == 'FissionProduct_ENDL' ) :
                iza = int( za[19:] )
            else :
                if( '_' in za ) :
                    ZA, suffix = za.split( '_' )
                else :
                    ZA, suffix = za, None
                if( suffix == 'natural' ) :
                    Z = ZA
                    A = 0
                else :
                    for i, s in enumerate( za ) :
                        if( s.isdigit( ) ) : break
                    Z, A = za[:i], za[i:]
                foundZA = False
                for iZ in xrange( fudgeZA.nZs ) :
                    if( Z == fudgeZA.ZToSymbol( iZ ) ) : 
                        iza = 1000 * iZ + int( A )
                        foundZA = True
                        break
                if( not foundZA ) : raise Exception( 'Could not handle za = %s.' % `za` )
        else :
            raise Exception( 'za = "%s" cannot be converted into an integer' % `za` )
    return( iza )

class bdfls :
    """This class reads in a bdfls file."""

    def __init__( self, name = None, template = None ) :
        """Reads in the bdfls file give by template. If self's save is called with name = None,
        then this methods name is used as the output bdfls file's name."""

        if ( template == None ) :
            if ( os.path.exists( "./bdfls" ) ) :
                template = "./bdfls"
            else :
                if( 'BDFLSPATH' in os.environ ) :
                    template = os.environ[ 'BDFLSPATH' ]
                else :
                    template = _BDFLS_FILE_
        if ( not os.path.exists( template ) ) : raise Exception( "\nError in bdfls.__init__: bdfls file does not exist\n%s" % template )
        self.name = name
        if ( name == None ) : self.name = "./" + template.split( "/" )[-1]
        self.template = template
        self.g = []
        self.f = []
        self.m = []
        self.massFormat = 'old'
        self.h = []
        self.c = []
        self.t = []
        self.s = []
        self.read( template )

    def __repr__( self ) :
        """Returns a string that summarizes the information in this class."""

        s = "File name = %s\n" % self.name
        s += "Template file = %s\n" % self.template
        s += "  Group data\n"
        for g in self.g : s += "    id = %4d  groups = %4d  boundaries = %d\n" % ( g.id, g.n - 1, g.n )
        s += "  Flux data\n"
        for f in self.f : s += "    id = %4d  lMax = %d\n" % ( f.id, f.lMax )
        s += "  There are %s masses\n" % len( self.m )
        s += "  There are %s halflifes\n" % len( self.h )
        s += "  There are %s constants\n" % len( self.c )
        s += "  Temperature data\n"
        for t in self.t : s = s + "    id = %4d  n = %d\n" % ( t.id, t.n )
        s += "  There are %s subshell designators\n" % len( self.s )
        return s

    def addGroup( self, ebs, id, label = "", overwrite = 0 ) :
        """Adds the group with id to self's list of groups. ebs is a python list
        of numbers and label is the groups label (see class bdfls_group). If overwrite
        is false, a raise is triggered if self already contains group with id, else
        the old one is overwritten."""

        pos = None
        for i in range( len( self.g ) ) :
            if( self.g[i].id == id ) :
                if( overwrite ) :
                    del self.g[i]
                    pos = i
                    break
                else :
                    raise Exception( "\nError in bdfls.addGroup: group with id = %d already exist" % id )
            elif( self.g[i].id > id ) :
                pos = i
                break
        g = bdfls_group( ebs )
        g.n = len( ebs )
        g.id = id
        g.label = label
        g.name = 'LLNL_gid_%d' % g.id
        if( pos == None ) :
            self.g.append( g )
        else :
            self.g.insert( pos, g )

    def addFlux( self, EF_ls, id, label = "", overwrite = 0 ) :
        """Adds the flux with id to self's list of fluxes. EF_ls is a python list
        with the following structure:
            [
                [ [E0, f0], [E1, f1], ... ], # for l = 0 
                [ [E0, f0], [E1, f1], ... ], # for l = 1 
                ...
            ]
        and label is the flux's label (see class bdfls_flux). If overwrite
        is false, a raise is triggered if self already contains group with id, else
        the old one is overwritten."""

        pos = None
        for i in range( len( self.f ) ) :
            if( self.f[i].id == id ) :
                if( overwrite ) :
                    del self.f[i]
                    pos = i
                    break
                else :
                    raise Exception( "\nError in bdfls.addFlux: flux with id = %d already exist" % id )
            elif( self.f[i].id > id ) :
                pos = i
                break
        f = bdfls_flux( EF_ls )
        f.id = id
        f.label = label
        f.name = 'LLNL_fid_%d' % f.id
        if( pos == None ) :
            self.f.append( f )
        else :
            self.f.insert( pos, f )

    def addZA( self, za, mass, halflife, warning = 1 ) :
        """Calls addZAHalflife with halflife and addZAMass with mass."""

        self.addZAHalflife( za, halflife, warning )
        self.addZAMass( za, mass, warning )

    def addZAHalflife( self, za, halflife, warning = 1 ) :
        """Adds or changes the halflife for ZA. If warning is true and ZA is already in self's halflife list,
        then a warning message is printed."""

        i = 0
        for h in self.h :
            if ( h.za == za ) :
                if ( warning ) : fudgemisc.printWarning( "changing halflife for za = %d from %s to %s" % ( za, `h.halflife`, `halflife` ) )
                h.halflife = halflife
                return
            elif ( h.za > za ) :
                self.h.insert( i, bdfls_halflife( za, halflife ) )
                return
            i += 1

    def addZAMass( self, za, mass, warning = 1 ) :
        """Adds or changes the mass for ZA. If warning is true and ZA is already in self's mass list,
        then a warning message is printed."""

        i = 0
        for m in self.m :
            if ( m.za == za ) :
                if ( warning ) : fudgemisc.printWarning( "changing mass for za = %d from %s to %s" % ( za, `m.mass`, `mass` ) )
                m.mass = mass
                return
            elif ( m.za > za ) :
                self.m.insert( i, bdfls_mass( za, mass, format = self.massFormat ) )
                return
            i += 1

    def flux( self, id ) :
        "Returns the flux that matches id."

        if( type( id ) == type( 1 ) ) :
            for f in self.f :
                if( f.id == id ) : return f
        else :
            for f in self.f :
                if( f.name == id ) : return f
        raise Exception( "\nError in bdfls.flux: id does not exist (%s)" % id )

    def group( self, id ) :
        "Returns the group that matches id."

        if( type( id ) == type( 1 ) ) :
            for g in self.g :
                if( g.id == id ) : return g
        else :
            for g in self.g :
                if( g.name == id ) : return g
        raise Exception( "\nError in bdfls.group: id does not exist (%s)" % id )

    def subGroups( self, id ) :
        "Returns the sub-group ids for group id."

        subgs = []
        for g in self.g :
            if( g.id == id ) : continue
            if( self.isSubGroup( g.id, id ) ) : subgs.append( g.id )
        subgs.sort( )
        return( subgs )

    def isSubGroup( self, sid_, id_, eps = 0. ) :
        """Returns true if sid is a sub-group of id, false otherwise."""

        sid = None
        id = None
        for g in self.g :
            if(   sid_ == g.id ) : sid = g
            elif( id_  == g.id ) : id = g
        if( id == None ) : raise Exception( "isSubGroup: group with id = %d does not exist" )
        if( sid == None ) : raise Exception( "isSubGroup: sub group with id = %d does not exist" )
        i = 0
        n = id.n
        for sg in sid.gb :
            while( 1 ) :
                if( n == i ) : return( 0 )
                d = sg - id.gb[i]
                if( abs( d ) <= eps ) :
                    i += 1
                    break
                if( d < 0. ) : return( 0 )
                i += 1
        return( 1 )

    def halflife( self, za ) :
        "Returns the halflife of the za if in database, else returns None."

        if( za == 'gamma' ) : return( 1e50 )
        za_ = targetToZA( za )
        iMin = 0
        iMax = len( self.h ) - 1
        while( True ) :
            iMid = ( iMin + iMax ) / 2
            h = self.h[iMid]
            if( h.za == za_ ) : return h.halflife
            if( iMin == iMax ) : break
            if( h.za > za_ ) :
                iMax = iMid
            else :
                iMin = iMid + 1
        return None

    def mass( self, za ) :
        "Returns the mass of the za if in database, else returns None (the mass of a za is for a neutral atom, not the nucleus)."

        if( za == 'gamma' ) : return( 0. )
        za_ = targetToZA( za )
        iMin = 0
        iMax = len( self.m ) - 1
        while( True ) :
            iMid = ( iMin + iMax ) / 2
            m = self.m[iMid]
            if( m.za == za_ ) : return m.mass
            if( iMin == iMax ) : break
            if( m.za > za_ ) :
                iMax = iMid
            else :
                iMin = iMid + 1
        return None

    def constant( self, i ) :
        "Returns the (i+1)^th constant from self's constant list."

        return( self.c[i].value )

    def printza( self, za, im = 0, ih = 0 ) :
        "Prints the za, mass and halflife for za."

        ( za, m, h ) = self.za( za, im, ih )
        if ( za == None ) : return
        sm = 15 * " "
        if( m != None ) : sm = "%15.8e" % m
        sh = 12 * " "
        if ( h != None ) : sh = "%12e" % self.h[ih].halflife
        print "%6d %s %s" % ( za, sm, sh )

    def read( self, name ) :
        """Reads in the bdfls file given by name. For internal use only."""

        class fileInfo :

            def __init__( self, name ) :
                "For internal use only."

                self.f = open( name, "r" )
                self.lineNumber = 0
                self.line = None
                self.linesRead = []

            def close( self ) :
                "For internal use only."

                self.f.close( )

            def readline( self ) :
                "For internal use only."

                if( len( self.linesRead ) > 0 ) :
                    self.line = self.linesRead[0]
                    del self.linesRead[0]
                else :
                    self.line = self.f.readline( )
                self.lineNumber += 1
                return( self.line )

            def unreadline( self, line ) :

                self.linesRead.append( line )
                self.lineNumber -= 1

            def errMsg( self, routine, str ) :
                "For internal use only."

                return( '\nError in %s: %s at line %d\n%s' % ( routine, str, self.lineNumber, self.line ) )

        if ( not os.path.exists( name ) ) : raise Exception( "\nError in bdfls.read: file does not exist (%s)" % name )
        f = fileInfo( name )
        self.Source = name
        while 1 :                               # Read the group data.
            g = bdfls_group( f )
            if ( g.id == None ) : break
            self.g.append( g )
        while 1 :                               # Read the flux data.
            fl = bdfls_flux( f )
            if ( fl.id == None ) : break
            self.f.append( fl )
        line = f.readline( )                    # Handle new mass format.
        f.unreadline( line )
        s = line.split( "#" )
        self.massFormat = 'old'
        if( len( s ) > 1 ) :
            self.massLabel = s[1][:-1]
            self.massFormat = 'Audi2003'
        while 1 :
            m = bdfls_mass( f, format = self.massFormat )
            if ( m.za == None ) : break
            self.m.append( m )
        while 1 :
            h = bdfls_halflife( f )
            if ( h.za == None ) : break
            self.h.append( h )
        ( self.nConstants, i ) = bdfls_read_int_label( f, 0, "bdfls.read", "constant", "number of constants" )
        i = self.nConstants
        while 1 :
            c = bdfls_constant( f )
            if ( c.value == None ) : break
            self.c.append( c )
            i = i - 1
        if ( i != 0 ) : raise Exception( "\nError in bdfls.read: %d unread constants" % i )
        while 1 :
            t = bdfls_temperature( f )
            if ( t.id == None ) : break
            self.t.append( t )
        ( self.nSubshells, i ) = bdfls_read_int_label( f, 0, "bdfls.read", "subshell", "number of subshell designators" )
        i = self.nSubshells
        while 1 :
            s = bdfls_subshell( f )
            if ( s.line == None ) : break
            self.s.append( s )
            i = i - 1
        if ( i != 0 ) : raise Exception( "\nError in bdfls.read: %d unread subshells" % i )
            
        f.close( )

    def save( self, fileName = None ) :
        """Saves the contents of self to fileName."""

        if ( fileName == None ) : fileName = self.name
        try :
            f = open( fileName, "w" )
        except :
            raise Exception( "\nError in bdfls.save: cannot open file %s" % fileName )
        s = `self.g`
        for g in self.g : f.write( `g` )
        f.write( bdfls_EOD )
        for fl in self.f : f.write( `fl` )
        f.write( bdfls_EOD )
        if( self.massFormat == 'old' ) :
            l = "                 nuclear masses                          *\n"
        elif( self.massFormat == 'Audi2003' ) :
            n = 0
            if( len( self.m ) > 0 ) :
                n = len( `self.m[0]` )
            l = '#' + self.massLabel
            n = 80 - ( len( l ) + n )
            l = n * ' ' + l + '\n'
        else :
            raise Exception( '\nError in bdfls.save: bad format = "%s"' % self.format )
        for m in self.m :
            f.write( `m` + l )
            l = "\n"
        f.write( bdfls_mh_EOD )
        l = "                 half-lives                              *\n"
        for h in self.h :
            f.write( `h` + l )
            l = "\n"
        f.write( bdfls_mh_EOD )
        f.write( "%4d         nuclear constants                                                 *\n" % self.nConstants )
        for c in self.c : f.write( `c` )
        f.write( bdfls_EOD )
        for t in self.t : f.write( `t` )
        f.write( bdfls_EOD )
        f.write( "%4d    subshell designators                                                   *\n" % self.nSubshells )
        for s in self.s : f.write( `s` )
        f.write( bdfls_EOD )
        f.close( )

    def setZAHalflife( self, za, halflife ) :
        """Changes the halflife for ZA. A raise is triggered is za in not in self's halflife list."""

        for h in self.h :
            if ( h.za == za ) :
                h.halflife = halflife
                return
            elif ( h.za > za ) :
                raise Exception( "\nError in bdfls.setZAHalflife: no such ZA = %s" % `za` )

    def setZAMass( self, za, mass ) :
        """Changes the mass for ZA. A raise is triggered is za in not in self's mass list."""

        for m in self.m :
            if ( m.za == za ) :
                m.mass = mass
                return
            elif ( m.za > za ) :
                raise Exception( "\nError in bdfls.setZAMass: no such ZA = %s" % `za` )

    def za( self, za, im = 0, ih = 0 ) :
        """Returns the tuple ( za, mass, halflife ) for za."""

        m = None
        lm = len( self.m )
        while ( im < lm ) :
            if ( self.m[im].za == za ) : m = self.m[im].mass
            if ( self.m[im].za >= za ) : break
            im = im + 1
        h = None
        lh = len( self.h )
        while ( ih < lh ) :
            if ( self.h[ih].za == za ) : h = self.h[ih].halflife
            if ( self.h[ih].za >= za ) : break
            ih = ih + 1
        if ( m == None ) and ( h == None ) : za = None
        return ( za, m, h )

    def zas( self, ZAMin, ZAMax = None ) :
        """Print za, mass and halflife for all za between ZAMin and ZAMax inclusive."""

        if ( ZAMax == None ) : ZAMax = ZAMin
        im = 0
        lm = len( self.m )
        ih = 0
        lh = len( self.h )
        while ( ( im < lm ) and ( self.m[im].za < ZAMin ) ) : im = im + 1
        while ( ( ih < lh ) and ( self.h[ih].za < ZAMin ) ) : ih = ih + 1
        za = ZAMin
        while za <= ZAMax :
            if ( ( im < lm ) and ( za > self.m[im].za ) ) : im = im + 1
            if ( ( ih < lh ) and ( za > self.h[ih].za ) ) : ih = ih + 1
            if ( ( im == lm ) and ( ih == lh ) ) : break
            if ( im == lm ) :
                za = self.h[ih].za
            elif ( ih == lh ) :
                za = self.m[im].za
            else :
                za = min( self.m[im].za, self.h[ih].za )
            if ( za <= ZAMax ) :
                self.printza( za, im, ih )
                za = za + 1

class bdfls_group :
    """A class that contains a single group."""

    def __init__( self, f ) :
        """f must be a python list of number for the group boundaries or an opened python file
        that is a bdfls file positioned at the start of a group."""

        self.id = None
        if( type( f ) == type( [] ) ) :
            self.label  = None
            self.n = len( f )
            self.gb = f
        else :
            self.read( f )

    def __len__( self ) :
        "Return the number of group boundaries."

        return( len( self.gb ) )

    def __getitem__( self, index ) :
        "Returns the (index - 1)^th boundary"

        return( self.gb[index] )

    def __repr__( self ) :
        "Returns a string of the group as it would appear in a bdfls file."

        s = "%4d %-74s*\n" % ( self.id, self.label )
        if ( len( s ) > 81 ) : s = s[:79] + "*\n"
        s = s + "%4d                                                                           *\n" % self.n
        j = 0
        for i in self.gb :
            s = s + "%12.5e" % i
            j = j + 1
            if ( ( j % 6 ) == 0 ) : s = s + "\n"
        if ( ( j % 6 ) != 0 ) : s = s + "\n"
        return s

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None,
        title = None, style = "lines" ) :
        """Plots the group boundaries.

xylog meanings::
    xylog   plot-type
   -----------------------
      0     linear-linear
      1     log-linear
      2     linear-log
      3     log-log"""

        if ( xLabel == None ) : xLabel = 'Energy (MeV)'
        if ( yLabel == None ) : yLabel = 'A. U.'
        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        for p in self.gb: f.write( "%15.7e 0.1\n%15.7e 10.\n\n" % ( p, p ) )
        f.close( )
        p = os.path.join( __file__.split( 'fudge/processing/' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        s = [ 'python', p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def read( self, f ) :
        """Reads in the next group from a bdfls file. Sets self's id to None if not more groups to read in."""

        ( self.id, self.label ) = bdfls_read_int_label( f, 1, "bdfls_group.read", "group", "group id" )
        if( self.id == None ) : return
        self.name = 'LLNL_gid_%d' % self.id
        ( self.n, l ) = bdfls_read_int_label( f, 0, "bdfls_group.read", "group", "number of groups" )
        self.gb = []
        i = 0
        while ( i < self.n ) :
            if ( ( i % 6 ) == 0 ) :
                l = f.readline( )
                if ( l == "" ) : raise Exception( f.errMsg( 'bdfls_group.__init__', 'EOF while reading group data' ) )
                if ( l == bdfls_EOD ) : raise Exception( f.errMsg( 'bdfls_group.__init__', 'EOD detected while reading group data' ) )
                j = 0
            self.gb.append( float( l[j:j+12] ) )
            j = j + 12
            i = i + 1

class bdfls_flux :
    """A class that contains a single flux."""

    def __init__( self, f ) :
        """f must be a python list of number for the with the following structure:
            [
                [ [E0, f0], [E1, f1], ... ], # for l = 0 
                [ [E0, f0], [E1, f1], ... ], # for l = 1 
                ...
            ]
        or an opened python file that is a bdfls file positioned at the start of a flux."""

        self.id = None
        if( type( f ) == type( [] ) ) :
            self.label  = None
            self.lMax = len( f ) - 1
            self.n = [ len(x) for x in f ]
            self.EF_l = f
        else :
            self.id = None
            self.read( f )

    def __len__( self ) :
        """Returns the number of Legendre orders for this flux."""

        return( len( self.EF_l ) )

    def __repr__( self ) :
        """Returns a string of the flux as it would appear in a bdfls file."""

        s = "%4d %-74s*\n" % ( self.id, self.label )
        if ( len( s ) > 81 ) : s = s[:79] + "*\n"
        s = s + "%4d                                                                           *\n" % self.lMax
        r = range( self.lMax + 1 )
        for i in r :
            s = s + "%4d                                                                           *\n" % ( 2 * self.n[i] )
        k = 0
        for i in r :
            for j in self.EF_l[i] :
                s = s + "%12.5e%12.5e" % ( j[0], j[1] )
                k = k + 1
                if ( ( k % 3 ) == 0 ) : s = s + "\n"
        if ( ( k % 3 ) != 0 ) : s = s + "\n"
        return s

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None,
        title = None, style = "lines" ) :
        """Plots the fluxes

xylog meanings::
    xylog   plot-type
   -----------------------
      0     linear-linear
      1     log-linear
      2     linear-log
      3     log-log"""

        if ( xLabel == None ) : xLabel = 'Energy (MeV)'
        if ( yLabel == None ) : yLabel = 'A. U.'
        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        for lGroup in self.EF_l:
            for p in lGroup: f.write( '   '.join( map( str, p ) ) + '\n' )
            f.write( '\n' )
        f.close( )
        p = os.path.join( __file__.split( 'fudge/processing/' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        s = [ 'python', p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def read( self, f ) :
        """Reads in the next flux from a bdfls file. Sets self's id to None if no more fluxes to read in."""

        ( self.id, self.label ) = bdfls_read_int_label( f, 1, "bdfls_flux.read", "flux", "flux id" )
        if( self.id == None ) : return
        self.name = 'LLNL_fid_%d' % self.id
        ( self.lMax, l ) = bdfls_read_int_label( f, 0, "bdfls_flux.read", "flux", "lMax" )
        self.EF_l = []
        self.n = []
        r = range( self.lMax + 1 )
        for i in r :
            ( n, l ) = bdfls_read_int_label( f, 0, "bdfls_flux.read", "flux", "number of flux points" )
            self.n.append( n / 2 )
        for i in r :
            j = 0
            EF = []
            while ( j < self.n[i] ) :
                if ( ( j % 3 ) == 0 ) :
                    l = f.readline( )
                    if ( l == "" ) : raise Exception( f.errMsg( 'bdfls_flux.__init__', 'EOD while reading flux data' ) )
                    if ( l == bdfls_EOD ) : f.errMsg( 'bdfls_flux.__init__', 'EOD detected while reading flux data' )
                    k = 0
                EF.append( [ float( l[k:k+12] ), float( l[k+12:k+24] ) ] )
                k = k + 24
                j = j + 1
            self.EF_l.append( EF )

class bdfls_mass :
    """A class that contains a single za and mass datum."""

    def __init__( self, f_or_za, mass = None, format = 'old' ) :
        """f_or_za must be an integer for the za or an opened python file
        that is a bdfls file positioned at a mass datum. If f_or_za
        is an integer then mass is its mass; otherwise, the za and mass are
        read in from the file."""

        if ( type( f_or_za ) == type( 1 ) ) :
            self.za = f_or_za
            self.mass = mass
        else :
            self.read( f_or_za )
        self.format = format

    def __repr__( self ) :
        """Returns a string of the za and mass as they would appear in a bdfls file."""

        if( self.format == 'old' ) :
            s = "    %6.4d" % self.za
            if ( self.mass == 0 ) :
                s = s + "    0."
            elif ( self.mass < 1. ) :
                s = s + "%14.8f" % self.mass
            elif ( self.mass > ( 1.e4 - 1e-6 ) ) :
                if   ( self.mass <= ( 1.e5  - 1e-6 ) ) :
                    s = s + "%12.5f" % self.mass
                elif ( self.mass <= ( 1.e6  - 1e-5 ) ) :
                    s = s + "%12.4f" % self.mass
                elif ( self.mass <= ( 1.e7  - 1e-4 ) ) :
                    s = s + "%12.3f" % self.mass
                elif ( self.mass <= ( 1.e8  - 1e-3 ) ) :
                    s = s + "%12.2f" % self.mass
                elif ( self.mass <= ( 1.e9  - 1e-2 ) ) :
                    s = s + "%12.1f" % self.mass
                elif ( self.mass <= ( 1.e10 - 1e-1 ) ) :
                    s = s + "%12.0f" % self.mass
                else :
                    s = s + "%12.5e" % self.mass
            else :
                s = s + "%12.6f" % self.mass
                if ( ( ( self.za % 1000 ) == 0 ) or ( self.za == 99120 ) or ( self.za == 99125 ) ) :
                    while ( s[-1] == '0' ) :
                        if ( s[-2] == "." ) : break
                        s = s[:-1]
        elif( self.format == 'Audi2003' ) :
            s = "    %6d" % self.za
            d = ' %17.12f' % self.mass
            try :
                while( d[-1] == '0' ) : d = d[:-1]
            except :
                fudgemisc.printWarning( d )
                raise
            s += d
        else :
            raise Exception( '\nError in bdfls.bdfls_mass.__repr__: bad format = "%s"' % self.format )
        return s

    def read( self, f ) :
        """Reads in the next mass datum from the bdfls file. Returns the tuple (za, mass) or
        (None, None) if at end of mass data."""

        ( self.za, self.mass ) = bdfls_read_int_float( f, "bdfls_mass.read", "mass" )

class bdfls_halflife :
    """A class that contains a single za and halflife datum."""

    def __init__( self, f_or_za, halflife = None ) :
        """f_or_za must be an integer for the za or an opened python file
        that is a bdfls file positioned at a halflife datum. If f_or_za
        is an integer then halflife is its halflife; otherwise, the za and halflife are
        read in from the file."""

        if ( type( f_or_za ) == type( 1 ) ) :
            self.za = f_or_za
            self.halflife = halflife
        else :
            self.read( f_or_za )

    def __repr__( self ) :
        """Returns a string of the za and halflife as they would appear in a bdfls file."""

        if ( self.halflife == None ) :
            s = "    %6d" % self.za
        else :
            s = "    %6d%12.4e" % ( self.za, self.halflife )
        return s

    def read( self, f ) :
        """Reads in the next halflife datum from the bdfls file. Set za and halflife 
        to None if at end of halflife data."""

        ( self.za, self.halflife ) = bdfls_read_int_float( f, "bdfls_halflife.read", "halflife" )

class bdfls_constant :
    """A class that contains a single constant datum."""

    def __init__( self, f ) :
        """f must be an opened python file that is a bdfls file positioned at a constant datum."""

        self.read( f )

    def __repr__( self ) :
        """Returns a string of the constant as it would appear in a bdfls file."""

        return "%15.8e $%s\n" % ( self.value, self.label )

    def read( self, f ) :
        """Reads in the next constand datum from the bdfls file. Sets self's value to None if at end of constant data."""

        self.value = None
        self.label = ""
        l = f.readline( )
        if ( l == bdfls_EOD ) : return
        i = string.find( l, "$" )
        if ( i == -1 ) : raise Exception( f.errMsg( 'bdfls_constant.read', 'missing $' ) )
        try :
            self.value = float( l[:i] )
        except :
            raise Exception( f.errMsg( 'bdfls_constant.read', 'converting constant to float' ) )
        self.label = l[i+1:-1]

class bdfls_temperature :
    """A class that contains a temperaure datum."""

    def __init__( self, f ) :
        """f must be an opened python file that is a bdfls file positioned at a temperature datum."""

        self.read( f )

    def __repr__( self ) :
        """Returns a string of the temperature data as they would appear in a bdfls file."""

        s = "%4d %-74s*\n" % ( self.id, self.label )
        if ( len( s ) > 81 ) : s = s[:79] + "*\n"
        s = s + "%4d                                                                           *\n" % self.n
        for i in self.ts : s = s + "%s\n" % i
        return s

    def read( self, f ) :
        """Reads in the next temperature datum from the bdfls file. Set self's id
        to None if at end of temerature data."""

        ( self.id, self.label ) = bdfls_read_int_label( f, 1, "bdfls_temperature.read", "temperature", "temperature id" )
        if( self.id == None ) : return
        ( self.n, l ) = bdfls_read_int_label( f, 0, "bdfls_temperature.read", "temperature", "number of temperatures" )
        self.ts = []
        i = 0
        while i < self.n :
            l = f.readline( )
            if ( l == "" ) : raise Exception( f.errMsg( 'bdfls_temperature.read', 'EOF while reading temperature data' ) )
            self.ts.append( l[:-1] )
            i = i + string.count( l, "," )
            if ( l[-2] != "," ) : break

class bdfls_subshell :
    """A class that contains the subshell data."""

    def __init__( self, f ) :
        """f must be an opened python file that is a bdfls file positioned at a subshell data."""

        self.read( f )

    def __repr__( self ) :
        """Returns a string of the subshell as it would appear in a bdfls file."""

        return( self.line )

    def read( self, f ) :
        """Reads in the subshell data from the bdfls file. Set self's line
        to None if at end of subshell data."""

        self.line = None
        l = f.readline( )
        if ( l == "" ) : raise Exception( f.errMsg( 'bdfls_subshell.read', 'EOF while reading subshell data' ) )
        if ( l != bdfls_EOD ) : self.line = l

def bdfls_read_int_label( f, allowEOD, callingRoutine, Str1, Str2 ) :

    l = f.readline( )
    if ( l == bdfls_EOD ) :
        if ( allowEOD ) :
            return ( None, None )
        else :
            raise Exception( f.errMsg( callingRoutine, 'EOD detected while reading %s' % Str1 ) )
    if ( l == "" ) : raise Exception( f.errMsg( callingRoutine, 'EOF while reading %s data' % Str1 ) )
    if ( len( l ) != 81 ) or ( l[79] != "*" ) : raise Exception( f.errMsg( callingRoutine, 'invalid %s line' % Str2 ) )
    try :
        n = int( l[:4] )
    except :
        raise Exception( f.errMsg( callingRoutine, 'converting %s to integer' % Str2 ) )
    return( n, l[5:-2] )

def bdfls_read_int_float( f, callingRoutine, Str ) :

    l = f.readline( )
    if ( l == bdfls_mh_EOD ) : return( None, None )
    if ( ( len( l ) == 81 ) and ( l[79] != "*" ) ) : raise Exception( f.errMsg( callingRoutine, 'invalid first za/%s line', Str ) )
    if ( l[:4] != "    " ) : raise Exception( f.errMsg( callingRoutine, 'bad za/%s line' % Str ) )
    s = l.split( )
    try :
        za = int( s[0] )
    except :
        raise Exception( f.errMsg( callingRoutine, 'converting za to integer' ) )
    d = None
    if ( len( s ) > 1 ) :
        try :
            d = float( s[1] )
        except :
            raise Exception( f.errMsg( callingRoutine, 'converting mass to float "%s"' % Str ) )
    return( za, d )

def getDefaultBdfls( name = None, template = None ) :

    global default_bdfls
    if( default_bdfls == None ) :
        if( name == None ) : name = "bdfls.junk"
        default_bdfls = bdfls( name = name, template = template )
    return( default_bdfls )

def getBdflsFile( name = None ) :

    return( bdfls( name = None, template = name ) )
