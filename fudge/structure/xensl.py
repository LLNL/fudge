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

try  :
    from xml.etree.cElementTree import ElementTree, tostring
except ImportError  :
    from ElementTree import ElementTree, tostring
import os
from fudge import fudgeDefaults
from fudge.particles import nuclear
import xml2py

__metaclass__ = type

global isomerMinimumHalflife
isomerMinimumHalflife = 1e-6 # in seconds

class XENSL( xml2py.XML2PY ) :
    """Class to manipulate the XENSL database found in fudgeDefaults.XENSL_DATABASE_DIR.
    This location can be overrode with the environment variable XENSLPATH.
    This class inherits from xml2py.XML2PY, overriding the read and write with automatic lookup of file based on ZA.
    """
    def __init__( self, ZA, database = None ) :
        xml2py.XML2PY.__init__( self )
        self.read( ZA, database = database )

    def read( self, ZA, database = None ) :
        if( database is None ) : database = fudgeDefaults.XENSL_DATABASE_DIR
        xmlfile = os.path.join( database, 'za%.6d.xml' % ZA )
        xml2py.XML2PY.read( self, xmlfile )

    def write( self ) :
        xmlfile = os.path.join( fudgeDefaults.XENSL_DATABASE_DIR, 'za%.6d.xml' % ( self.Z * 1000 + self.A ) )
        xml2py.XML2PY.write( self, xmlfile )

class Levels :
    """Class containing the levels for one nucleus.
    The levels can be accessed via their index in the list or via an attribute e0, e1, e2, etc.
    where 0 is the ground state, 1 the first excited state and so on.
    A dict of metastable states is stored in self.metastables with keys m0, m1, m2, etc. which
    map to the associated level e?. m0 is the ground state if this is a metastable state, 
    otherwise the first metastable state is labeled m1. Note stable states are not labeled metastable.
    The metastable states can be defined using the method defineMetastable( seconds ) with the
    argument defining the minimum halflife for a metastable state.
    Default is to use the global variable isomerMinimumHalflife which defaults to 1e-6s.
    If any isomers are listed in XENSL database, then these overide the default cutoff at 1e-6s.
    The metastable states can be also looked up via an attribute m1, m2 etc.
    A call to this class with a numerical argument will return the level with the nearest energy, 
    while a call with a string argument, e.g. 'e1', 'm1' will look up the corresponding level.
    """
    def __init__( self, isotope, database = None ) :
        try  :
            self.ZA = int( isotope )
            self.isotope = nuclear.nucleusNameFromZA( self.ZA )
        except  :
            self.isotope = isotope
            self.ZA = nuclear.getZ_A_suffix_andZAFromName( self.isotope )[3]
        self.xensl = XENSL( self.ZA, database = database )
        self.levels = self.xensl.levels
        self.metastables = {} # example {'m1':'e2'}
        self.levelMetastables = {}      # i.e., { 'e2' : 'm1' }
        try :
            self.getMetastables( )
        except AttributeError :
            self.defineMetastables( isomerMinimumHalflife )

    def __getattr__( self, attr ) :
        try :
            if 'e' in attr or 'm' in attr[0] and attr[1:].isdigit( ) :
                return self.getLevel( attr )
            else : raise TypeError
        except TypeError:
            return object.__getattribute__( self, attr )

    def __getitem__( self, index ) :
        return self.levels[index]

    def __repr__( self ) :
        return 'Metastable levels for ZA=%s\n%4s    Energy   Halflife  Metastable\n' % ( self.ZA, 'id' ) + str( self.toString( True ) )

    def __str__( self ) :
        return 'Levels for ZA=%s\n%4s    Energy   Halflife  Metastable\n' % ( self.ZA, 'id' ) + str( self.toString( ) )

    def __call__( self, arg ) :
        """If arg is a float (i.e. an energy in MeV) then returns getNearestLevelAtEnergy(arg).
        If arg is a str ( i.e. 'e2' or 'm1' ) then returns getLevel( arg )."""
        try :
            return self.getNearestLevelAtEnergy( float( arg ) )
        except ValueError :
            try :
                if isinstance( arg, str ): return self.getLevel( arg )
            except TypeError : pass
        raise TypeError, "Bad call for '%s': must be an energy (in MeV) or level identifier (e.g. 'e1','m2')" % arg

    def toString( self, metastablesOnly=False ) :
        s=[]
        metastableLevels = {}
        for key, val in self.metastables.items( ) :
            assert val not in metastableLevels
            metastableLevels[val] = key
        for l in self.levels :
            if l.label in metastableLevels :
                s.append( '%4s  %s  %s' % ( l.label, self.levelString( l ), metastableLevels[l.label] ) )
            elif not metastablesOnly :
                s.append( '%4s  %s' % ( l.label, self.levelString( l ) ) )
        return '\n'.join( s )

    def levelString( self, level ) :
        halflife = level.xml.get( 'halflife', '' )
        if halflife == '-1.0' : halflife = 'stable'
        return '%10.6f %8s' % ( level.energy, halflife )

    def getLevel( self, i ) :
        """Returns the level (depends on type of argument)
        integer : returns the level with index i
        str     : 'e1' 'm1' returns the level via energy label or metastable label
        """
        if isinstance( i, int ): return self.levels[i]
        if isinstance( i, str ) :
            if i[1:].isdigit( ) :
                if 'e' in i[0]: return self.levels[int( i[1:] )]
                if 'm' in i[0]: return self.getMetastable( i )
            elif i=='m':
                if 'm2' not in self.metastables:
                    return self.getMetastable( 'm1' )
                else:
                    raise ValueError, "More than one isomer exists, please specify which isomer, e.g. m1,m2\n%s" % ( repr( self ) )
        raise TypeError, "%s needs to be either integer of level or 'e1','m1'" % i

    def getMetastable( self, l ) :
        """Returns the metastable level with label l. e.g. 'm1' if in list of known metastables"""
        if l in self.metastables: return self.getLevel( self.metastables[l] )
        else: raise KeyError, '%s not a known metastable state: %s' % ( l, ' '.join( self.metastables.keys( ) ) )

    def getNearestLevelAtEnergy( self, energy, printWarning = True, onMiddleChoose = 0 ) :
        """Returns the level with an energy nearest to the energy given as an argument. If the energy is
        exactly between two levels and onMiddleChoose is 0, a raise is executed. Otherwise, if 
        onMiddleChoose is less ( greater ) than 0 the lower ( upper ) level is returned."""

        if( energy < 0 ) : raise Exception( 'Negative excitation energy = %s not allowed' % energy )
        if( len( self.levels ) == 1 ) : return self.levels[0]
        energy = float( energy )
        l1 = self.levels[0]
        for i in xrange( 1, len( self.levels ) ) :
            l2 = self.levels[i]
            if( l1.energy <= energy <= l2.energy ) :
                printReturnedLevel = False
                average = 0.5 * ( l1.energy + l2.energy )
                if( 0 < ( abs( energy - average ) / ( l2.energy - l1.energy ) ) < 0.2 ) :
                    if( printWarning ) :
                        print 'WARNING: Energy %s is close to midpoint between levels, level returned not unambiguous.' % energy
                        print ' Possible levels: %s  %s' % ( l1.energy, l2.energy )
                        printReturnedLevel = True
                if( energy < average ) :
                    if( printReturnedLevel ) : print '  Returned level: %s' % ( l1.energy )
                    return l1
                elif( energy > average ) :
                    if( printReturnedLevel ) : print '  Returned level: %s' % ( l2.energy )
                    return l2
                else :
                    if( onMiddleChoose == 0 ) :
                        raise RuntimeError, 'Energy %s is exactly at midpoint between levels\n%s\n%s' % ( energy, l1, l2 )
                    elif( onMiddleChoose < 0 ) :
                        return( l1 )
                    else :
                        return( l2 )
            l1 = l2
        if( energy < 1.1 * l1.energy ) : return( l1 )
        raise RuntimeError, 'No nearest level found for %s at E = %s' % ( self.isotope, energy )

    def defineMetastables( self, minHalflife ) :
        """Defines the metastable states, filling the dict mapping the metastable states to the
        actual level, for levels with halflives greater than minHalflife. If ground state is
        metastable then the labeling begins with 'm0', else it starts with 'm1'."""
        self.metastables.clear( )
        metastableCounter = 0
        for i, level in enumerate( self.levels ) :
            if i > 0 and metastableCounter == 0: metastableCounter += 1
            if float( level.xml.get( 'halflife', 0 ) ) > minHalflife :
                self.metastables['m%s' % metastableCounter] = 'e%s' % i
                metastableCounter += 1
        self.defineLevelMetastables( )

    def defineLevelMetastables( self ) :
        self.levelMetastables.clear( )
        for m in self.metastables : self.levelMetastables[self.metastables[m]] = m

    def getMetastables( self ) :
        """Gets the metastable states from XENSL,
        filling the dict mapping the metastable states to the actual level."""
        self.metastables.clear( )
        self.metastables.update( dict( [( key, val ) for key, val in self.xensl.isomers.__dict__.items( ) if key != 'xml'] ) )
        self.defineLevelMetastables( )

    def setIsomerMinimumHalflife( self, halflife ) :
        """Sets the minimum halflife for defining an isomer.
        Recalculates isomer labels for level scheme from this new definition."""
        global isomerMinimumHalflife
        isomerMinimumHalflife = halflife
        self.defineMetastables( isomerMinimumHalflife )

    def getIsomerMinimumHalflife( self ) :
        """Returns the minimum halflife for defining an isomer"""
        global isomerMinimumHalflife
        return isomerMinimumHalflife

if __name__ == '__main__' :
    import fudge
    x=XENSL( 1001 )
    x=Levels( 36079 )
    print "\nCreate the object containing the nucleus structure\n>>>x=Levels(36079)\n>>>x"
    print repr( x )
    print "\nFind level nearest to 0.12 MeV\n>>>x.getNearestLevelAtEnergy(0.12)"
    print x.getNearestLevelAtEnergy( 0.12 )
    print "\nFind level nearest to 0.12 MeV\n>>>x(0.12)"
    print x( 0.12 )
    print "\nReturn metastable level m1\n>>>x.m1"
    print x.m1
    print "\nReturn metastable level m1\n>>>x('m1')"
    print x( 'm1' )
    print "\nReturn 4th level\n>>>x.e4"
    print x.e4
    print "\nLook up level near energy (level returned not unambiguous)\n>>>x(0.135)"
    print x( 0.135 )
    print "\nLook up level exactly at midpoint\n>>>x(0.138415)"
    try: print x( 0.138415 )
    except RuntimeError, err: print 'RuntimeError: %s' % err
    print '\nSome other well known nuclei with isomers'
    print 'Hf178'
    print repr( Levels( 72178 ) )
    print 'Am242'
    print repr( Levels( 95242 ) )
