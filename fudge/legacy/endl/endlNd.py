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
This module contains the class endlNd which contains methods common to data of any (N) dimension.
"""
import os
import copy
import sys
import time

import bdfls
import endlmisc
import endl2

class endlNd :
    """
Useful Members (all of these should be accessed via get/set methods)::

    h               The ENDL header for this data.
    ZA              The ENDL ZA for this data. That is, the target for this file.
    yi              The ENDL yi-value for this data. That is, the projectile for this file.
    yo              The ENDL yo-value for this data. That is, the outgoing particle for this file.
    C               The ENDL C-value for this data.
    I               The ENDL I-value for this data.
    S               The ENDL S-value for this data.
    X1, X2, X3, X4  The X-vales for this data.
    ELevel          The initial excitation level of the target for this data.
    Halflife        The halflife of the target.
    Temperature     The temperature of the target for this data.
"""

    def __init__( self, f, I_, yo, C, I, S, h, points, i0 = 0, i1 = 1, i2 = None, i3 = None, bdflsFile = None ) :
        """For internal use only."""

        if ( I != I_ ) : raise Exception( "\nError in endlI%d.__init__: I = %d != %d" % ( I_, I, I_ ) )
        if( bdflsFile == None ) : bdflsFile = bdfls.getDefaultBdfls( )
        self.bdflsFile = bdflsFile
        self.columns = endlmisc.getNumberOfColumns_( I, "endlNd.__init__" )
        self.yo = yo
        self.C = C
        self.I = I
        self.S = S
        self.format = 12
        if ( f == None ) :                              # Designer data.
            self.h = copy.deepcopy( h )
            data = points
            self.designerData = 1
        elif ( type( f ) == file ) :    # Otherwise, data must be in an endl formatted file.
            self.h = [ f.readline( ), "" ]
            if ( self.h[0] == "" ) : return
            self.h[1] = f.readline( )
            if ( self.h[1] == "" ) : raise Exception( "\nError in endlNd.__init__: end-of-file while reading from %s" % f.name )
            data = readNdEndlData( f, i0, i1, i2, i3 )  # File data.
            self.designerData = 0
        else :
            raise Exception( "\nError from endlNd.__init__: invalid type = %s for f" % type( f ) )
        if( not( ( self.h[0] == "\n" ) and ( C == 1 ) and ( I == 0 ) ) ) :
            if( self.yo != int( self.h[0][9:12] ) ) : raise Exception( "\nError in endlNd.__init__: yo in header differs from requested yo\n" + `self.h` )
            if( self.C  != int( self.h[1][ :2]  ) ) : raise Exception( "\nError in endlNd.__init__: C in header differs from requested C\n" + `self.h` )
            if( self.I  != int( self.h[1][2:5] ) ) : raise Exception( "\nError in endlNd.__init__: I in header differs from requested I\n" + `self.h` )
            if( self.S  != int( self.h[1][5:8] ) ) : raise Exception( "\nError in endlNd.__init__: S in header differs from requested S\n" + `self.h` )
            try :
                self.ZA = int( self.h[0][:6] )
                self.yi = int( self.h[0][6:9] )
                self.Mass = endlmisc.headerString2FunkyDouble( self.h[0],  13, "endlNd.__init__" )
                if( self.h[0][31] == ' ' ) :
                    self.ENDLInterpolation = 0
                else :
                    self.ENDLInterpolation = int( self.h[0][31] )
                self.format = int( self.h[0][32:34] )
                self.ELevel = endlmisc.headerString2FunkyDouble( self.h[0],  35, "endlNd.__init__" )
                self.Halflife = endlmisc.headerString2FunkyDouble( self.h[0],  47, "endlNd.__init__" )
                self.Temperature = endlmisc.headerString2FunkyDouble( self.h[0],  59, "endlNd.__init__" )
                self.Q = 0.
            except :
                endlmisc.printWarning( self.h[0] )
                endlmisc.printWarning( self.h[1] )
                raise
            if ( len( self.h[1] ) >= 20 ) : self.Q = endlmisc.headerString2FunkyDouble( self.h[1],  9, "endlNd.__init__" )
            self.X1 = 0.
            if ( len( self.h[1] ) >= 32 ) : self.X1 = endlmisc.headerString2FunkyDouble( self.h[1], 21, "endlNd.__init__" )
            self.X2 = 0.
            if ( len( self.h[1] ) >= 44 ) : self.X2 = endlmisc.headerString2FunkyDouble( self.h[1], 33, "endlNd.__init__" )
            self.X3 = 0.
            if ( len( self.h[1] ) >= 56 ) : self.X3 = endlmisc.headerString2FunkyDouble( self.h[1], 45, "endlNd.__init__" )
            self.X4 = 0.
            if ( len( self.h[1] ) >= 68 ) : self.X4 = endlmisc.headerString2FunkyDouble( self.h[1], 57, "endlNd.__init__" )
        else :
            if( not ( self.designerData ) ) : endlmisc.printWarning( "Warning: using old style C = 1, I = 0 header file" )
            self.ZA = None
            self.yi = None
            self.yo = 0
            self.C = 1
            self.I = 0
            self.S = 0
            self.Mass = None
            self.ENDLInterpolation = 0
            self.ELevel = None
            self.Halflife= None
            self.Temperature = None
            self.Q = 0.
            self.X1 = 0.
            self.X2 = 0.
            self.X3 = 0.
            self.X4 = 0.
        self.set( data, checkDataType = 0, interpolation = endlmisc.endlToFudgeInterpolation( self.ENDLInterpolation ) )

    def __len__( self ) :
        """Returns the length of data. For 2d data this is the number of data pairs, for 3d data this is the
        number of unique x values and for 4d data this is the number of unique t values."""

        return len( self.data )

    def __repr__( self ) :
        """Returns a string containing yo, C, I, S, X1, X2, X3, X4 and Q. To get a string of the actual data
        use the toString method for that data."""

        return "yo = %2d  C = %3d  I = %3d  S = %3d  X1 = %e  X2 = %e  X3 = %e  X4 = %e  Q = %e" % \
            ( self.yo, self.C, self.I, self.S, self.X1, self.X2, self.X3, self.X4, self.Q )

    def checkHeader( self, index = None ) :
        "For internal use only."

        if( ( index == None ) or ( index == 1 ) ) :
            if( len( self.h[0] ) < 71 ) :
                if( self.h[0][-1] == "\n" ) : self.h[0][:-1]
                self.h[0] = "%-70s\n" % self.h[0]
        if( ( index == None ) or ( index == 2 ) ) :
            if( len( self.h[1] ) < 69 ) :
                if( self.h[1][-1] == "\n" ) : self.h[1][:-1]
                self.h[1] = "%-68s\n" % self.h[1]

    def cmpyoCISX1( self, other, epsX1 ) :
        """Compares yo, C, I, S and X1 values of self to other and returns 1 if they
    are the same and 0 otherwise. The X1 values can differ by as much as epsX1."""

        v1 = 1000 * ( 100 * ( 100 * self.yo + self.C ) + self.I ) + self.S
        v2 = 1000 * ( 100 * ( 100 * other.yo + other.C ) + other.I ) + other.S
        if ( v1 < v2 ) :
            return -1
        elif ( v1 > v2 ) :
            return 1
        else :
            dX1 = self.X1 - other.X1
            sX1 = ( abs( self.X1 ) + abs( other.X1 ) ) * abs( epsX1 )
            if ( abs( dX1 ) <= sX1 ) :
                return 0
            elif ( dX1 < 0 ) :
                return -1
            else :
                return 1

    def compareQXs( self, Q, X1, X2, X3, X4 ) :
        """Return -1 if self's values < (Q, X1, X2, X3, X4), 0 if equal and 1 otherwise."""

        if( Q is not None ) :
            if( self.Q < Q ) : return(  1 )
            if( self.Q > Q ) : return( -1 )
        if( self.X1 < X1 ) : return( -1 )
        if( self.X1 > X1 ) : return(  1 )
        if( self.X2 < X2 ) : return( -1 )
        if( self.X2 > X2 ) : return(  1 )
        if( self.X3 < X3 ) : return( -1 )
        if( self.X3 > X3 ) : return(  1 )
        if( self.X4 < X4 ) : return( -1 )
        if( self.X4 > X4 ) : return(  1 )
        return( 0 )

    def determineQandThreshold( self, Q = None, specialCases = 0 ) :

        if( ( self.C < 8 ) or ( 49 < self.C < 58 ) ) : return( 0., 0. )
        if( Q is None ) :
            Q = endl2.reactionQ( self.yi, self.ZA, self.C, targetELevel = self.getELevel( ), X4 = self.getX4( ), specialCases = specialCases, bdflsFile = self.bdflsFile )
            if( ( self.ZA == 95241 ) and ( self.C == 46 ) and ( self.Q == 5.48 ) ) : Q += 5.48 - 5.4318 # Special case for endl99
        if( Q is None ) : Q = self.getQ( )
        if( self.S in [ 1 ] ) :
            residualELevel = self.getX1( )
        else :
            residualELevel = self.getX4( )
        threshold = endl2.reactionThreshold( self.yi, self.ZA, self.C, targetELevel = self.getELevel( ), residualELevel = residualELevel, Q = Q,
            specialCases = specialCases, S = self.S, bdflsFile = self.bdflsFile )
        return( Q, max( 0., threshold ) )

    def setQandThreshold( self, thresholdCrossSectionIsZero, Q, threshold, fixThresholdMode = None, dThreshold_MeV = 10e-3, 
            EMin = 0., threshold_MeV_shiftWarning = 0.1 ) :

        self.setQ( Q )
        if( hasattr( self, 'fixThreshold' ) ) : self.fixThreshold( thresholdCrossSectionIsZero, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin,
            fixThresholdMode = fixThresholdMode, threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getDate( self ) :
        "Returns the target's date."

        return( int( self.h[0][25:31] ) )

    def getELevel( self ) :
        "Returns the target's excitation level."

        return( self.ELevel )

    def getENDLInterpolation( self ) :
        "Returns the target's ENDL interpolation flag. 0 or 2 is lin-lin, 3 is log-lin, 4 is lin-log and 5 is log-log. Also see getInterpolation."

        return( self.ENDLInterpolation )

    def getInterpolation( self ) :
        "Returns the target's fudge interpolation flag. 0 is lin-lin, 1 is log-lin, 2 is lin-log and 3 is log-log. Also see getENDLInterpolation."

        return( self.interpolation )

    def getFormat( self ) :
        "Returns the target's format specifier."

        return( self.format )

    def getLifetime( self ) :
        "Deprecated method. Use getHalflife instead."

        endlmisc.printWarning( 'getLifetime is a deprecated method. Use getHalflife instead.' )
        return( self.getHalflife( ) )

    def getHalflife( self ) :
        "Returns the target's halflife."

        return( self.Halflife )

    def getMass( self ) :
        "Returns the target's mass."

        return( self.Mass )

    def getQ( self ) :
        "Returns the Q value for self's reaction."

        return( self.Q )

    def getTemperature( self ) :
        "Returns the taget's temperature."

        return( self.Temperature )

    def getX1( self ) :
        "Returns the X1 value for this data."

        return( self.X1 )

    def getX2( self ) :
        "Returns the X2 value for this data."

        return( self.X2 )

    def getX3( self ) :
        "Returns the X3 value for this data."

        return( self.X3 )

    def getX4( self ) :
        "Returns the X4 value for this data."

        return( self.X4 )

    def getYi( self ) :
        "Returns projectile's identifier, yi value."

        return( self.yi )

    def getYo( self ) :
        "Returns the outgoing particle's identifier, yo value."

        return( self.yo )

    def getZA( self ) :
        "Returns the targets ZA value."

        return( self.ZA )

    def info( self ) :
        "Prints information about self's data."

        print "C = %2d  I = %3d  S = %3d  X1 = %e  X2 = %e  X3 = %e  X4 = %e  Q = %e" % \
            ( self.C, self.I, self.S, self.X1, self.X2, self.X3, self.X4, self.Q )
        print self.h[0], self.h[1],
        print "Data has %d columns and %d rows" % ( self.columns, len( self.data ) )

    def save( self, f = None ) :        # ??????
        """Writes header and data to file f where f is a python file type (e.g., sys.stdout).
    If f == None then writes data to sys.stdout."""

        if ( f == None ) : f = sys.stdout
        f.write( self.h[0] )
        f.write( self.h[1] )
        Fmt = None
        if( self.format != 12 ) :
            Fmt = '%%%d.%de' % ( self.format, self.format - 7 )
        s = self.toString( format = Fmt )
        f.write( s )
        f.write( "                                                                       1\n" )

    def setC( self, C ) :
        "Sets self's C value and self's header C value to C."

        self.checkHeader( 2 )
        self.C = C
        self.h[1] = "%2d" % C + self.h[1][2:]

    def setDate( self, date = None ) :
        "Sets self's data and self's header date value to data. If date == None then today's date is used."

        if ( date == None ) :
            date = time.localtime( time.time( ) )
            date = 10000 * ( date[0] - 2000 ) + 100 * date[1] + date[2]
        if ( type( date ) == type( "" ) ) :
            try :
                date = int( date )
            except :
                endlmisc.printWarning( "date cannot be converted to an integer" )
                raise ValueError
        if ( type( date ) == type( 1 ) ) and ( date > 0 ) and ( date < 1000000 ) :
            self.h[0] = "%-25s%6.6d%s" % ( self.h[0][0:25], date, self.h[0][31:] )
        else :
            raise Exception( "\nError in setDate: date must be an integer, or a string convertible to an \ninteger, in the range [1, 999999]" )

    def setELevel( self, ELevel ) :
        "Sets self's ELevel value and self's header ELevel value to ELevel."

        self.checkHeader( 1 )
        self.ELevel = ELevel
        self.h[0] = self.h[0][:35] + endlmisc.headerFunkyDouble2String( ELevel ) + self.h[0][46:]

    def setFormat( self, width ) :
        "Sets self's format specifier and self's header format specifier value to width."

        if( ( width < 7 ) or ( width > 99 ) ) : raise Exception( '\nError in endlNd.setFormat: width = %d not in range [7 to 99].' % width )
        self.checkHeader( 1 )
        self.format = width
        self.h[0] = '%s%2d%s' % ( self.h[0][:32], width, self.h[0][34:] )

    def setI( self, I ) :
        "Sets self's I value and self's header I value to I."

        self.checkHeader( 2 )
        self.I = I
        self.h[1] = self.h[1][:2] + "%3d" % I + self.h[1][5:]

    def setLifetime( self, Halflife) :
        "Deprecated method. Use setHalflife instead."

        endlmisc.printWarning( 'setLifetime is a deprecated method. Use setHalflife instead.' )
        self.setHalflife( Halflife )

    def setHalflife( self, Halflife) :
        "Sets self's Halflife value and self's header Halflife value to Halflife."

        self.checkHeader( 1 )
        self.Halflife = Halflife
        self.h[0] = self.h[0][:47] + endlmisc.headerFunkyDouble2String( Halflife ) + self.h[0][58:]

    def setMass( self, Mass ) :
        "Sets self's Mass value and self's header Mass value to Mass."

        self.checkHeader( 1 )
        self.Mass = Mass
        self.h[0] = self.h[0][:13] + endlmisc.headerFunkyDouble2String( Mass ) + self.h[0][24:]

    def setQ( self, Q ) :
        "Sets self's Q value and self's header Q value to Q."

        self.checkHeader( 2 )
        self.Q = Q
        self.h[1] = self.h[1][:9] + endlmisc.headerFunkyDouble2String( Q ) + self.h[1][20:]

    def setS( self, S ) :
        "Sets self's S value and self's header S value to S."

        self.checkHeader( 2 )
        self.S = S
        self.h[1] = self.h[1][:5] + "%3d" % S + self.h[1][8:]

    def setTemperature( self, T ) :
        "Sets self's Temperature value and self's header Temperature value to T."

        self.checkHeader( 1 )
        self.Temperature = T
        self.h[0] = self.h[0][:59] + endlmisc.headerFunkyDouble2String( T ) + self.h[0][70:]

    def setX1( self, X1 ) :
        "Sets self's X1 value and self's header X1 value to X1."

        self.checkHeader( 2 )
        self.X1 = X1
        self.h[1] = self.h[1][:21] + endlmisc.headerFunkyDouble2String( X1 ) + self.h[1][32:]

    def setX2( self, X2 ) :
        "Sets self's X2 value and self's header X2 value to X2."

        self.checkHeader( 2 )
        self.X2 = X2
        self.h[1] = self.h[1][:33] + endlmisc.headerFunkyDouble2String( X2 ) + self.h[1][44:]

    def setX3( self, X3 ) :
        "Sets self's X3 value and self's header X3 value to X3."

        self.checkHeader( 2 )
        self.X3 = X3
        self.h[1] = self.h[1][:45] + endlmisc.headerFunkyDouble2String( X3 ) + self.h[1][56:]

    def setX4( self, X4 ) :
        "Sets self's X4 value and self's header X4 value to X4."

        self.checkHeader( 2 )
        self.X4 = X4
        self.h[1] = self.h[1][:57] + endlmisc.headerFunkyDouble2String( X4 ) + self.h[1][68:]

    def setYi( self, yi ) :
        """Sets self's yi value and self's header yi value to yi.  yi must be an integer
    or a string convertible to an integer."""

        iyi = endlmisc.incidentParticleTags( yi )[0]
        if( ( iyi < 1 ) or ( iyi > 10 ) ) : raise Exception( "\nError in setYi: yi must be in the range [1, 10]" )
        self.checkHeader( 1 )
        self.yi = iyi
        syi = "%3d" % iyi
        self.h[0] = self.h[0][:6] + syi + self.h[0][9:]

    def setYo( self, yo ) :
        """Sets self's yo value and self's header yo value to yo.  yo must be an integer
    or a string convertible to an integer."""

        iyo = int( yo )
        if( ( iyo < 0 ) or ( iyo > 19 ) ) : raise Exception( "\nError in setYo: yo must be in the range [0, 19]" )
        self.checkHeader( 1 )
        self.yo = iyo
        syo = "%3d" % iyo
        self.h[0] = self.h[0][:9] + syo + self.h[0][12:]

    def setZA( self, ZA ) :
        """Sets self's ZA value and self's header ZA value to ZA.  ZA must be an integer
    or a string convertible to an integer."""

        iZA = int( ZA )
        if( ( iZA < 0 ) or ( iZA > 999999 ) ) : raise Exception( "\nError in setZA: ZA must be in the range [0, 999999]" )
        self.checkHeader( 1 )
        self.ZA = iZA
        sZA = "%6d" % iZA
        self.h[0] = sZA + self.h[0][6:]

    def setYiZAYoCIS( self, yi = None, ZA = None, yo = None, C = None, I = None, S = None ) :
        """For input parameters that are not None, sets self's value to input parameter."""

        if( yi != None ) : self.setYi( yi )
        if( ZA != None ) : self.setZA( ZA )
        if( yo != None ) : self.setYo( yo )
        if( C  != None ) : self.setC( C )
        if( I  != None ) : self.setI( I )
        if( S  != None ) : self.setS( S )

def readNdEndlData( f, i0, i1, i2 = None, i3 = None ) :
    "For internal use only."

    data = []
    while 1 :
        l = f.readline( )
        if ( l == "" ) : raise Exception( "\nError: end-of-file while reading from %s" % f.name )
        if ( l == "                                                                       1\n" ) : break
        data.append( l )
    if ( i3 != None ) : return endlmisc.translate4dStringData( data, i0 = i0, i1 = i1, i2 = i2, i3 = i3 )
    elif ( i2 != None ) : return endlmisc.translate3dStringData( data, i0 = i0, i1 = i1, i2 = i2 )
    else : return endlmisc.translate2dStringData( data, i0 = i0, i1 = i1 )
