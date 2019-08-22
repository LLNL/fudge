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
This module contains the class endlFile. See class endlFile for usage.
"""

import os
import time
from fudge.core.utilities import fudgeExceptions
import bdfls
import endl2
import endlIClasses
import endlmisc

checkQThreshold = True

class endlFile :
    """This class holds data that is contained in an ENDL file. (In general, this class should only be instantiated indirectly via methods in endlZA.)
An ENDL data file can contains one or more related data sets, each specified by a unique set of X1, X2, X3, X4 values, for a specific set of 
yo, C, I and S values.

As example, for the cross section, I = 0, for inelastic scattering 
of a neutron, C = 11, the target may be left in different excited states. In ENDL ths cross sections to specific excited state, S = 1, are stored in a file
named "yo00c11i000s001", and the data for different excited states are distinguished by different excitation levels, X1-value, for the
target. When the read method for this class gets the data from this file it creates a new endlI0 object for each data set and addes it
to the levels member of this class. The number of data sets (i.e., length of levels) is return by the len method for the endlFile object.

No data is read in when the object is instantiated; instead, the read method must be explicitly called.

Useful Members::
    name        The file's name.
    path        The full path, include name, of this file.
    ZA          The ENDL ZA-value for this file. ZA = 1000 * Z + A.
    yi          The ENDL yi-value for this file. That is, the projectile for this file.
    yo          The ENDL yo-value for this file. That is, the outgoing particle for this file.
    C           The ENDL C-value for this file.
    I           The ENDL I-value for this file.
    S           The ENDL S-value for this file.
    levels      A list of endlI* objects for each data sets in file. This is empty until the method read is called.
"""

    def __init__( self, designerFile, name, path, ZA, yi, halflife = None, bdflsFile = None ) :
        """In general, this routine should be used indirectly via endlZA.read( ) or endlZA.addFile( ) and not called directly, 
        unless you know what you are doing."""

        if( bdflsFile is None ) : bdflsFile = bdfls.getDefaultBdfls( )
        self.bdflsFile = bdflsFile
        self.designerFile = designerFile
        self.name = name
        self.path = path
        self.ZA = ZA
        self.yiTags = endlmisc.incidentParticleTags( yi )
        self.yi = self.yiTags[0]
        self.yo = int( name[2:4] )
        self.C = int( name[5:7] )
        self.I = int( name[8:11] )
        self.S = int( name[12:15] )
        self.mode = endlmisc.getmodestring_( self.I, "endlFile.__init__" )
        self.halflife = halflife                            # Added so uses can input for natural ZA because the halflife is not in the bdfls file.
        self.levels = []

    def __len__( self ) :
        """Returns the number of levels for self."""

        return len( self.levels )

    def __getitem__( self, index ) :
        """Returns the (i-1)^th item from self.levels."""

        if( self.levels == [] ) : raise Exception( "File not read in." )
        return( self.levels[index] )

    def __repr__( self ) :
        """Returns the string 'R# filename' where R is 'r' if the data has been read in and is '-' 
    if it has not been, '#' is the number of columns of data and file name is the name of the
    file (e.g. 'r2 yo00c01i000s000')."""

        r = "r"
        if( self.levels == [] ) : r = "-"
        return r + self.mode + " " + self.name

    def check( self, targetParameters, dQ_MeV = 10e-3, dThreshold_MeV = 10e-3, allowZeroE = False, xCloseEps = None, normTolerance = 1e-5,
        maxAbsFloatValue = None ) :
        """Checks self's data if it as been read in. This is, it goes through the list levels, calling the
        appropriate routines to check data."""

        suffix = os.path.basename(self.path.rstrip(self.name).rstrip(os.sep))[8:]
        ErrMsgs = []
        X1 = X2 = X3 = X4 = None
        for l in self.levels :
            if( ( ( self.C == 1 ) and ( self.I == 0 ) ) or ( self.C == 5 ) ) :
                Threshold = 0.
                Q = 0.
            elif( self.C == 6 ) :
                Threshold = 0.
                Q = 0.
            elif( ( endl2.AFromZA( l.ZA ) == 0 ) or ( 99120 <= l.ZA <= 99200 ) or ( self.C == 15 ) or ( self.C >= 50 ) ) :
                Q = l.getQ( )
            else :
                Q = endl2.reactionQ( l.yi, l.ZA, l.C, targetELevel = l.getELevel( ), X4 = l.getX4( ), specialCases = 1, bdflsFile = self.bdflsFile )
                if( l.S in [ 1 ] ) :
                    residualELevel = l.getX1( )
                else :
                    residualELevel = l.getX4( )
                Threshold = endl2.reactionThreshold( l.yi, l.ZA, l.C, targetELevel = l.getELevel( ), residualELevel = residualELevel, Q = l.getQ( ), S = self.S,
                        bdflsFile = self.bdflsFile )
            r = 1. + self.bdflsFile.mass( self.yi ) / self.bdflsFile.mass( self.ZA )
            QThreshold = ( -l.getQ( ) ) * r
            if( QThreshold <= 0. ) : QThreshold = 0.
            if( targetParameters['T'] is None ) :
                targetParameters['T'] = l.Temperature
                targetParameters['ELevel'] = l.ELevel
                targetParameters['halflife'] = l.Halflife
            if( ( targetParameters['T'] != l.Temperature ) or ( targetParameters['ELevel'] != l.ELevel ) or ( targetParameters['halflife'] != l.Halflife ) ) :
                ErrMsgs.append( endlmisc.endlCheckerObject( suffix = suffix, data = l, message = 'target parameters differ' ) )
            elif( X1 is not None ) :
                if( ( X1 == l.X1 ) and ( X2 == l.X2 ) and ( X3 == l.X3 ) and ( X4 == l.X4 ) ) :
                    ErrMsgs.append( endlmisc.endlCheckerObject( suffix = suffix, data = l, message = 'duplicate data' ) )
            X1 = l.X1; X2 = l.X2; X3 = l.X3; X4 = l.X4
            if( abs( l.getQ( ) - Q ) > dQ_MeV ) :
                ErrMsgs.append( endlmisc.endlCheckerObject( suffix = suffix, data = l, message = 'bad Q: should be %.4e' % Q ) )
            if( ( self.C != 6 ) and( hasattr( l, "getThresholdsForChecker" ) ) ) :
                thresholds = l.getThresholdsForChecker( )
                iT = 0
                _dThreshold_MeV = dThreshold_MeV
                for t in thresholds :   # Currently, thresholds is the first min( n, 2 ) Energy points in the data where n is the length of the data.
                    if( checkQThreshold ) :
                        if( t + _dThreshold_MeV < QThreshold ) : ErrMsgs.append( endlmisc.endlCheckerObject( suffix = suffix, data = l, \
                            message = 'bad QThreshold: should be %.8e at index = %d instead of %.8e' % ( QThreshold, iT, t ) ) )
                    else :
                        if( t + _dThreshold_MeV < Threshold ) : ErrMsgs.append( endlmisc.endlCheckerObject( suffix = suffix, data = l, \
                            message = 'bad threshold: should be %.8e at index = %d instead of %.8e' % ( Threshold, iT, t ) ) )
                    if( l.I != 4 ) : _dThreshold_MeV = 0.   # The second point must be higher than the threshold else ndfgen may fail.
                    iT += 1
            ErrMsgs += l.check( printWarning = False, printErrors = False, allowZeroE = allowZeroE, xCloseEps = xCloseEps, normTolerance = normTolerance, \
                maxAbsFloatValue = maxAbsFloatValue )
        return( ErrMsgs )

    def setYiZAYoCISs( self, yi = None, ZA = None, yo = None, C = None, I = None, S = None ) :
        """Fixes the headers in all of file's levels to match yo, C, I and S. For input
        parameters that are None, self's values are used. Also sets self's yi and ZA
        to inputted values if they are other than None."""

        if( yi is None ) : yi = self.yi
        if( ZA is None ) : ZA = self.ZA
        if( yo is None ) : yo = self.yo
        if( C  is None ) : C  = self.C
        if( I  is None ) : I  = self.I
        if( S  is None ) : S  = self.S
        self.yi = yi
        self.ZA = ZA
        self.yo = yo
        self.C = C
        self.I = I
        self.S = S
        self.name = 'yo%.2dc%.2di%.3ds%.3d' % ( yo, C, I, S )
        for l in self.levels : l.setYiZAYoCIS( yi = yi, ZA = ZA, yo = yo, C = C, I = I, S = S )

    def IsyoCIS( self, yo = None, C = None, I = None, S = None ) :
        """Returns 1 if self matches yo, C, I and S, and 0 otherwise.  yo may be a single integer or a list of integers and 
        the same for C, I and S. For example self has yo = 7 and C = 3 then self.IsyoCIS( yo = (5, 6, 7), C=(3,4) ) will return a true."""

        if( yo is not None ) :
            if( ( type( yo ) == type( [] ) ) or ( type( yo ) == type( () ) ) ) :
                yo2 = yo
                yo = []
                for y in yo2 : yo.append( endlmisc.outgoingParticleTags( y )[0] )
            else :
                yo = endlmisc.outgoingParticleTags( yo )[0]
        if ( ( type( yo ) != type( [] ) ) and ( type( yo ) != type( () ) ) ) : yo = ( yo, )
        if ( ( type(  C ) != type( [] ) ) and ( type(  C ) != type( () ) ) ) :  C = (  C, )
        if ( ( type(  I ) != type( [] ) ) and ( type(  I ) != type( () ) ) ) :  I = (  I, )
        if ( ( type(  S ) != type( [] ) ) and ( type(  S ) != type( () ) ) ) :  S = (  S, )
        if ( ( self.yo in yo ) or ( yo == ( None, ) ) ) and \
           ( (  self.C in  C ) or (  C == ( None, ) ) ) and \
           ( (  self.I in  I ) or (  I == ( None, ) ) ) and \
           ( (  self.S in  S ) or (  S == ( None, ) ) ) : return 1
        return 0

    def unreadData( self, X1 = None, X2 = None, X3 = None, X4 = None, Q = None ) :
        """Removes self's reference to data matching X1, X2, X3, X4 and Q. Note, the data may still be in the work directory if
        the file is not saved."""

        i = len( self.levels )
        while( i > 0 ) :
            i -= 1
            l = self.levels[i]
            if ( ( l.X1 == X1 ) or ( X1 is None ) ) and ( ( l.X2 == X2 ) or ( X2 is None ) ) and \
               ( ( l.X3 == X3 ) or ( X3 is None ) ) and ( ( l.X4 == X4 ) or ( X4 is None ) ) and \
                ( ( l.Q == Q ) or ( Q is None ) ) :
                del self.levels[i]

    def addData( self, data = [], date = None, interpolation = 0, eLevel = 0., temperature = 2.586e-8, Q = None, 
        X1 = 0., X2 = 0., X3 = 0., X4 = 0., mass = None, halflife = None ) :
        """Adds data to self with the header information provided."""

        i = 0
        for l in self.levels :
            c = l.compareQXs( Q, X1, X2, X3, X4 )
            if( c == 0 ) :
                raise Exception( "\nError in endlFile.addData: data already exist for X1 = %e, X2 = %e, X3 = %e, X4 = %e, Q = %e" % ( X1, X2, X3, X4, l.Q ) )
            elif ( c > 0 ) :
                break
            i += 1
        if( Q is None ) :
            Q = 0.
            if ( len( self.levels ) > 0 ) : Q = self.levels[0].Q
        if( date is None ) : 
            date = time.localtime( time.time( ) )
            date = 10000 * ( date[0] - 2000 ) + 100 * date[1] + date[2]
        if( type( date ) == type( 1 ) ) : date = "%.6d" % date
        if( mass is None ) : mass = self.bdflsFile.mass( self.ZA )
        if( mass is None ) : raise Exception( "\nError in endlFile.addData: bdfls file does not have mass for ZA = %d" % self.ZA )
        mass = endlmisc.headerFunkyDouble2String( mass )
        if( halflife is None ) : halflife = self.bdflsFile.halflife( self.ZA )
        if( halflife is None ) : halflife = self.halflife
        if( halflife is None ) : raise Exception( "\nError in endlFile.addData: bdfls file does not have halflife for ZA = %d" % self.ZA )
        halflife = endlmisc.headerFunkyDouble2String( halflife )
        h = [ "%6d %2d %2d %s %6s%1d12 %s %s %s\n" % ( self.ZA, self.yi, self.yo, mass, date, interpolation, 
                        endlmisc.headerFunkyDouble2String( eLevel ), halflife, endlmisc.headerFunkyDouble2String( temperature ) ),
              "%2d%3d%3d %s %s %s %s %s\n" % ( self.C, self.I, self.S, endlmisc.headerFunkyDouble2String( Q ), 
                        endlmisc.headerFunkyDouble2String( X1), endlmisc.headerFunkyDouble2String( X2), 
                        endlmisc.headerFunkyDouble2String( X3), endlmisc.headerFunkyDouble2String( X4 ) ) ]
#
#  ???? Need to check that data is proper type of data.
#
        d = endlIClasses.endlAddIObject( None, self.yo, self.C, self.I, self.S, h, data, bdflsFile = self.bdflsFile )
        self.levels.insert( i, d )
        return self.levels[i]

    def addEndlData( self, endlData, overWrite = 0 ) :
        """Adds the data in the endl object endlData to self. enldData must be of an endlI* class that must match self's
        I class. For example, if self is I = 11 data then endlData be must an endlI11 object. Currently, overWrite is not
        implemented."""    # THIS NEEDS WORKD????

        if( ( endlData.I != self.I ) or ( endlData.S != self.S ) or ( endlData.C != self.C ) or ( endlData.yo != self.yo ) ) :
            raise Exception( "( endlData.I != self.I ) or ( endlData.S != self.S ) or ( endlData.C != self.C ) or ( endlData.yo != self.yo )" )
        i = 0
        for l in self.levels :
            c = l.compareQXs( endlData.Q, endlData.X1, endlData.X2, endlData.X3, endlData.X4 )
            if( c == 0 ) :
                raise Exception( "Data already exist for X1 = %e, X2 = %e, X3 = %e, X4 = %e, Q = %e for yo = %d C = %d, I = %d and S =%d" 
                    % ( l.X1, l.X2, l.X3, l.X4, l.Q, self.yo, self.C, self.I, self.S ) )
            elif( c > 0 ) :
                break
            i += 1
        points = endlData.copyData( )
        o = endlIClasses.endlAddIObject( None, endlData.yo, endlData.C, endlData.I, endlData.S, endlData.h, points, bdflsFile = self.bdflsFile )
        self.levels.insert( i, o )
        return o

    def columns( self ) :
        """Returns the number of data columns for self's data. For example, for 3d data a 3 is returned."""

        return endlmisc.getNumberOfColumns_( self.I, "endlFile.columns" )

    def info( self ) :
        """Prints information about each level in self (i.e., call info for each element of self.levels."""

        for i in range( len( self.levels ) ) :
            print "levels[%d]" % i
            self.levels[i].info( )
            print

    def ls( self ) :
        """Prints file name, and level information for self."""

        print  self.name
        s = ""
        for i in range( len( self.levels ) ) : s = s + ( " levels[%d]  %s" % ( i, `self.levels[i]` ) )
        print s

    def read( self ) :
        """Reads in the data from self's file."""

        if ( not os.path.isfile( self.path ) ) : raise Exception( "\nError in endlFile.read: no such file %s" % self.path )
        f = open( self.path )
        level = 1
        while 1 :
            d = endlIClasses.endlAddIObject( f, self.yo, self.C, self.I, self.S, None, [], bdflsFile = self.bdflsFile )
            if ( d.h[0] == "" ) : break
            slevel = None
            if( ( self.C == 11 ) and ( self.S == 1 ) ) : slevel = level
            d.label = endl2.reactionEquations( self.yi, self.ZA, self.C, level = slevel, printQWarning = False )[1]
            self.levels.append( d )
            level += 1
        f.close( )

    def save( self, f ) :
        """Saves data to file f where f is a python file type (e.g., sys.stdout)."""

        for l in self.levels : l.save( f )
