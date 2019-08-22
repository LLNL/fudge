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
This module contains the particle classes.
"""

import re
import xData.ancestry as ancestryModule
from pqu import PQU
from fudge.particles import nuclear
from fudge.gnd import warning

__metaclass__ = type

families = 0

def addFamily( _class, override = False ) :

    global families

    if( not( hasattr( _class, 'order' ) ) or override ) :
        _class.order = families
        families += 1

class base( ancestryModule.ancestry ) :
    """
    Base class for particles and elements. Every appearing in the gnd particleList must be an instance of base.
    """

    def __init__( self, name, attributes = None ):

        ancestryModule.ancestry.__init__( self )
        self.__name = name
        self.attributes = attributes or {}
        self.configurations = {}

    def __str__( self ) :

        return( self.name )

    def addConfiguration( self, label, configuration ) :

        self.configurations[label] = configuration.__copy__( )

    @property
    def name( self ) :

        return self.__name

    @property
    def id( self ) :

        return self.__name

class particle( base ):
    """
    Base class for particles.
    """

    def __eq__( self, other ) :
        """Returns True if other is same class as self and self._localCompare( other ) == 0 and False otherwise (also see _localCompare)."""

        return( self._localCompare( other, testEqual = True ) == 0 )

    def __ge__( self, other ) :
        """If other is same class as self returns result of self._localCompare( other ) > 0. Else executes a raise (also see _localCompare)."""

        return( self._localCompare( other ) >= 0 )

    def __gt__( self, other ) :
        """If other is same class as self returns result of self._localCompare( other ) => 0. Else executes a raise (also see _localCompare)."""

        return( self._localCompare( other ) > 0 )

    def __le__( self, other ) :
        """If other is same class as self returns result of self._localCompare( other ) <= 0. Else executes a raise (also see _localCompare)."""

        return( self._localCompare( other ) <= 0 )

    def __lt__( self, other ) :
        """If other is same class as self returns result of self._localCompare( other ) < 0. Else executes a raise (also see _localCompare)."""

        return( self._localCompare( other ) < 0 )

    def __ne__( self, other ) :
        """Returns True if other is same class as self and self._localCompare( other ) != 0 and False otherwise (also see _localCompare)."""

        return( self._localCompare( other, testEqual = True ) != 0 )

    def _localCompare( self, other, testEqual = False ) :

        if( isinstance( other, particle ) ) :
            if( self.family == other.family ) : return( self._localCompare2( other ) )
            if( hasattr( self, 'order' ) and hasattr( other, 'order' ) ) :
                return( self.order - other.order )      # Family with lower order is higher on sort ordering.
            elif( hasattr( self, 'order' ) ) :
                return( -1 )
            elif( hasattr( other, 'order' ) ) :
                return(  1 )
            else :
                if( self.family < other.family ) : return( -1 )
                return( 1 )
        elif( testEqual ) :
            return( 1 )
        else :
            raise TypeError( 'other must be an instance of particle: it is type "%s"' % type( other ) )

    def getFamily( self ) :

        return( self.family )

    def getMass( self, unit ) :

        return( self.mass.getValueAs( unit ) )

    def getParity( self ) :

        return( self.parity.__copy__( ) )

    def getSpin( self ) :

        return( self.spin.__copy__( ) )

    def getZ_A_SuffixAndZA( self ) :

        if( not( isinstance( self, ( isotope, nuclearLevel ) ) ) ) : return( 0, 0, '', 0 )
        Z, A, suffix, ZA = nuclear.getZ_A_suffix_andZAFromName( self.name )
        return( Z, A, suffix, 1000 * Z + A )

class isotope( particle ) :
    """This is the particle class for an isotope."""

    moniker = 'isotope'
    family = moniker

    def __init__( self, name, mass = None, attributes = None ) :

        particle.__init__( self, name, attributes = attributes )
        self.setMass( mass )
        self.levels = {}

    def __getitem__( self, key ) :

        return( self.levels[key] )

    def __iter__( self ) :

        particleKeys = sorted( self.levels.keys( ) )
        for particleKey in particleKeys :
            yield self.levels[particleKey]

    def _localCompare2( self, other ) :
        """
        Designed only for internal use. Returns -1 if self is less than other, 0 if it is the same and 1 
        otherwise. Compares self to other based firstly on their ZAs (ZA = 1000 * Z + A).  If that comparison 
        differs, the comparison is returned. Otherwise, returns the comparison of self's and other's suffixes.
        """

        Z, A, suffixSelf,  ZASelf  = self.getZ_A_SuffixAndZA()
        Z, A, suffixOther, ZAOther = other.getZ_A_SuffixAndZA()
        if( ZASelf < ZAOther ) : return( -1 )
        if( ZASelf > ZAOther ) : return(  1 )
        if( suffixSelf < suffixOther ) : return( -1 )
        if( suffixSelf > suffixOther ) : return(  1 )
        return( 0 )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""

        if( self.attributes['mass'] is None ) : return( None )
        return( self.attributes['mass'].getValueAs( unit ) )

    def hasID( self, ID ) :

        for level in self :
            if( level.name == ID ) : return( True )
        return( False )

    def setMass( self, mass ) :
        """Sets mass of the particle. Mass must be either None or a PQU instance."""

        if( mass is not None ) :
            if( not( isinstance( mass, PQU.PQU ) ) ) : raise TypeError( 'mass must be None or PQU instance. Its type is "%s"' % type( mass ) )
        self.attributes['mass'] = mass

    def getSpin( self ):
        """Return particle spin, or None if the spin is unknown"""

        if 'spin' in self.attributes: return self.attributes['spin']
        if self.levels: return self.levels[0].getSpin()
        return None

    def getParity( self ):
        """Return particle parity, or None if parity is unknown"""

        if 'parity' in self.attributes: return self.attributes['parity']
        if self.levels: return self.levels[0].getParity()
        return None
    
    def addLevel( self, level ):

        if isinstance(level, nuclearLevel):
            level.groundState = self
            self.levels[ level.label ] = level
            level.setAncestor( self, 'name' )
        else :
            raise TypeError( 'level not an instance of class nuclearLevel: class = "%s".' % level.__class__ )

    def check( self, info ) :

        warnings = []

        emax = -1
        for level in sorted(self.levels):
            levelWarnings = self.levels[ level ].check( info )
            if levelWarnings:
                warnings.append( warning.context('Level %s' % level, levelWarnings) )

            if type(level) is not int: continue # 'c' and 'u' levels may be out of order
            enow = self.levels[level].getLevelAsFloat('eV')
            if enow <= emax:
                warnings.append( warning.discreteLevelsOutOfOrder( level, self.levels[level] ) )
            else:
                emax = enow

        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        qs = ''                             # order attributes:
        attrs = self.attributes.copy()
        for attr in ( 'mass', 'spin','parity' ) :
            if( attr == 'mass' ) :
                if attr in attrs: qs += ' %s="%s"' % ( attr, attrs.pop( attr ).toString( keepPeriod = False ) )
            else :
                if attr in attrs: qs += ' %s="%s"' % (attr, attrs.pop(attr))
        for q in attrs:
            if type(attrs[q]) is bool:
                qs += ' %s="%s"' % (q, str(attrs[q]).lower())
            else:
                qs += ' %s="%s"' % ( q, attrs[q] )

        xml = [ '%s<%s name="%s"%s' % ( indent, self.family, self.name, qs ) ]
        if self.levels:
            xml[-1] += '>'
            for level in sorted(self.levels.keys()):
                xml.extend( self.levels[level].toXMLList( indent2, **kwargs ) )
            xml[-1] += '</%s>' % self.family
        else :
            xml[-1] += '/>'
        return xml

class nuclearLevel( particle ) :

    moniker = 'level'
    family = moniker

    def __init__( self, name, energy, label, gammas = None, attributes = None, groundState = None ) :

        particle.__init__( self, name, attributes = attributes )
        self.energy = energy
        if( type( label ) != type( 1 ) ) :
            if( label not in [ 'c', 's' ] ) : raise Exception( 'invalid level = "%s" for particle "%s"' % ( label, name ) )
        self.label = label      # could be integer, 'c' for continuum or 's' for sum.
        self.gammas = gammas or []
        self.groundState = groundState  # points to isotope instance

    def _localCompare2( self, other ) :
        """
        Designed only for internal use. Returns -1 if self is less than other, 0 if it is the same and 1 otherwise. 
        Compares self to other based firstly on their ground state comparison (see isotope._localCompare2).  If 
        that comparison differs, the comparison is returned. Otherwise, returns the comparison of self's and other's labels.
        """

        if( self.groundState < other.groundState ) : return( -1 )
        if( self.groundState > other.groundState ) : return(  1 )
        if( self.label  < other.label ) : return( -1 )
        if( self.label == other.label ) : return(  0 )
        return( 1 )
 
    def addGamma( self, gamma ) :

        if not isinstance( gamma, nuclearLevelGamma ):
            raise TypeError( 'gamma must be xParticle.nuclearLevelGamma instance and not class "%s"' % gamma.__class__ )
        self.gammas.append( gamma )
        gamma.setAncestor( self, 'gamma' )

    def getMass( self, unit ):
        """Get the mass (as a float) in the specified unit.
        Result includes ground-state mass + level excitation energy."""

        return self.groundState.getMass( unit ) + self.energy.getValueAs( unit+'*c**2' )
    
    def getLevelIndex( self ) :

        return self.label

    def getLevelAsFloat( self, unit, default = 0. ) :

        try :
            level = float( self.energy.getValueAs( unit ) )
        except :    # for continuum and sum, self.energy = None
            if( default is None ) :
                raise Exception( 'Could not convert level energy = "%s" to float' % self.energy )
            level = default
        return( level )
    
    def getSpin( self ) :

        return self.attributes.get( 'spin' )

    def getParity( self ) :

        return self.attributes.get( 'parity' )

    def check( self, info ):
        '''
        Check parts of nuclearLevel
        
        Options & defaults:
            'branchingRatioSumTolerance'      1e-6
        '''
        warnings = []
        gammaBRSumAbsTol = info.get( 'branchingRatioSumTolerance', 1e-6 )
        if self.gammas and abs( sum([gamma.probability for gamma in self.gammas]) - 1.0 ) > gammaBRSumAbsTol:
            warnings.append( warning.unnormalizedGammas( sum([gamma.probability for gamma in self.gammas]), self ) )
        return warnings

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        qs = ''
        # order attributes:
        attrs = self.attributes.copy()
        for attr in ('energy','spin','parity'):
            if attr in attrs: qs += ' %s="%s"' % (attr, attrs.pop(attr))
        for q in attrs: qs += ' %s="%s"' % (q, attrs[q])
        energy = self.energy.toString( keepPeriod = False )
        xml = [ '%s<level name="%s" label="%s" energy="%s"%s' % (indent, self.name, self.label, energy, qs) ]
        if self.gammas:
            xml[-1] += '>'
            for gamma in self.gammas:
                xml.extend( gamma.toXMLList( indent2, **kwargs ) )
            xml[-1] += '</level>'
        else: xml[-1] += '/>'
        return xml
    
    def getGammaEmission( self ) :
        """
        This method returns a dictionary of all gammas (as nuclearLevelGamma) and their probabilities for the
        current level and all levels in its decay path. Each gamma has a unique (for this target) key in the
        dictionary and the returned dictionary can be combined with those returned from other levels from self's 
        isotope (i.e., parent). The probability for each gamma is the sum from all decay paths. As example, 
        consider the decay from level 5 where the possible decays are (when more than one 'to' level, probability 
        'p' is given as (p)):

                 level 
              from    to
                5     1 (1/2) and 3 (1/2)
                4     0
                3     0 (1/4) and 2 (3/4)
                2     0 (2/3) and 1 (1/3)
                1     0
                0

        In this example, the gamma from level 1 to 0 occurs with probability 1/2 + 1/2 * 3/4 * 1/3 = 5/8.

        Two returned dictionaries name GGE1 and GGE2 can be combined into GGE1 as:

        for key in GGE2 :
            if( key not in GGE1 ) :
                GGE1[key] = copy.copy( GGE2[key] )
            else :
                GGE1[key][1] += GGE2[key][1]
        """

        def _addGamma( levelFrom, gamma, gammas, probability, level ) :
            """Designed for internal use only."""

            levelTo = gamma.finalLevel
            name = levelFrom.name + '__' + levelTo.name
            probability *= gamma.probability
            if( name not in gammas ) : gammas[name] = [ gamma, 0 ]
            gammas[name][1] += probability
            for gamma_ in levelTo.gammas : _addGamma( levelTo, gamma_, gammas, probability, level + 1 )

        gammas = {}
        for gamma in self.gammas : _addGamma( self, gamma, gammas, 1., 0 )
        return( gammas )

class FissionProduct( isotope ) :
    """
    Special case for storing 'representative' fission product information.
    This is partly for backwards-compatibility with ENDL databases, which collect all data
    for prompt and delayed fission products into '99120' and '99125' respectively.
    """

    moniker = 'FissionProduct'
    family = moniker

class thermalNeutronScatteringLawIsotope( isotope ) :
    """
    Paricle family for thermal neutron scattering law particle. These particles are nuclei whose scattering
    cross section depend on the material that they are associated with. For example at protron (i.e., H1)
    in water that is given the name "H1_inH2O_TNSL". All particles must have the name "P_M_TNSL" where
    "P" is the isotope name for the particle (i.e., "H1" in the example above), "M" is the material
    descriptor (i.e., "inH2O" in the example above) and the string "TNSL" is short for
    Thermal Neutron Scattering Law.
    """

    moniker = 'thermalNeutronScatteringLawIsotope'
    family = moniker

    def getZ_A_SuffixAndZA( self ) :

        ZA = { 'H1' : 1901, 'H2' : 1902, 'Be' : 4909, 'C1' : 6912, 'O1' : 8916 }[self.name[:2]]
        Z = ZA // 1000
        A = ZA % 100
        if( ZA == 1901 ) :
            if( 'H2O' in self.name ) : ZA -= 100
        elif( ZA == 4909 ) :
            if( 'Metal' in self.name ) : ZA -= 100
        return( Z, A, '', ZA )

class lepton( particle ) :

    moniker = 'lepton'
    family = moniker
    generations = { 'electronic' : 1, 'muonic' : 2,  'tauonic' : 3 }
    chargeOrder = { -1 : 0, 1 : 1, 0 : 2 }

    def __init__( self, name, generation, mass, charge, anti = False, attributes = None ) :

        particle.__init__( self, name, attributes = attributes )
        self.mass = mass
        self.charge = charge
        if( generation not in self.generations ) : raise Exception( 'Invalid generation "%s" for particle "%s" of family "%s"' % \
                ( generation, name, self.family ) )
        self.generation = generation
        self.anti = anti
        self.spin = spin( 0.5 )
        self.parity = parity( 1 )

    def _localCompare2( self, other ) :
        """For internal use only."""

        if( self.generation != other.generation ) : return( abs( self.generations[other.generation] ) - abs( self.generations[self.generation] ) )
        if( self.charge != other.charge ) : return( self.chargeOrder[other.charge] - self.chargeOrder[self.charge] )
        if( self.anti != other.anti ) :
            if( self.anti ) : return( -1 )
            return( 1 )
        return( 0 )

    def check( self, info ) :

        warnings = []
        return( warnings )

    def toXMLList( self, indent = '', **kwargs ) :

        quantities = ' mass="%s" spin="%s" charge="%s" generation="%s"' % ( self.mass, self.spin, self.charge, self.generation )
        attributes = self.attributes.copy()
        for quantity in attributes :
            value = attributes[quantity]
            if( isinstance( attributes[quantity], bool ) ) : value = str( value ).lower()
            quantities += ' %s="%s"' % ( quantity, value )
        xml = [ '%s<%s name="%s"%s/>' % ( indent, self.family, self.name, quantities ) ]
        return( xml )

class photon( particle ) :

    moniker = 'photon'
    family = moniker

    def __init__( self, name = "gamma", attributes = None ) :

        particle.__init__( self, name, attributes = attributes )
        self.mass = PQU.PQU( 0, "amu" )
        self.spin = spin( 1. )
        self.parity = parity( 1 )

    def _localCompare2( self, other ) :
        """For internal use only."""

        return( 0 )

    def check( self, info ) :

        warnings = []
        return( warnings )

    def toXMLList( self, indent = '', **kwargs ) :

        quantities = ' mass="%s" spin="%s"' % ( self.mass, self.spin )
        attributes = self.attributes.copy()
        for quantity in attributes :
            if( quantity in [ 'mass', 'spin' ] ) : continue
            value = attributes[quantity]
            if( isinstance( attributes[quantity], bool ) ) : value = str( value ).lower()
            quantities += ' %s="%s"' % ( quantity, value )
        return( [ '%s<%s name="%s"%s/>' % ( indent, self.moniker, self.name, quantities ) ] )

class nuclearLevelGamma( ancestryModule.ancestry ) :

    moniker = 'gamma'

    def __init__( self, finalLevel, angularDistribution, probability, nonRadiativeProbability = 0.0 ) :

        ancestryModule.ancestry.__init__( self )
        self.finalLevel = finalLevel
        self.probability = probability
        self.nonRadiativeProbability = nonRadiativeProbability
        self.angularDistribution = angularDistribution

    def __str__( self ) :

        return( 'Gamma energy = %s, final level = %s, gammaEmissionProbability = %s, internalConversionProbability = %s'
                % ( self.getEnergy(), self.finalLevel, self.probability, self.nonRadiativeProbability ) )

    def getEnergy( self, unit = 'eV', correctForRecoil = True ) :
        """Calculate emitted gamma energy (including recoil correction) in specified units."""

        energy = self.getAncestor().getLevelAsFloat( unit ) - self.finalLevel.getLevelAsFloat( unit )
        if( correctForRecoil ) : energy -= 0.5 * energy**2 / self.finalLevel.getMass( unit+'/c**2' )
        return( energy )

    def toXMLList( self, indent = '', **kwargs ) :

        xml =  '%s<%s finalLevel="%s" probability="%s"' % ( indent, self.moniker, self.finalLevel, self.probability )
        if( self.nonRadiativeProbability ) : xml += ' nonRadiativeProbability="%s"' % self.nonRadiativeProbability
        xml += '/>'
        return [xml]

class spin:
    """Class for storing nuclear spin (integer or half-integer).
    Initialize either with j1=spin(2.5) or j1=spin('5/2')
    str(j1) returns '5/2' """

    _re = re.compile(r'^([0-9]+)/?(2?)$')
    def __init__(self, spin):
        if spin is None: self.value = spin
        elif type(spin) == str:
            match = re.search(self._re, spin)
            if not match:
                try:
                    self.value=float(spin)
                except:
                    raise ValueError( "Cannot decipher spin from string: '%s'" % spin )
            else:
                num,den = match.groups()
                #self.spin = float(num)/(2 if den else 1) # this is such a pretty construction that I hate to remove it for Python 2.4
                if den: self.value = float(num)/(2)
                else:   self.value = float(num)/(1)
        else:
            if not ( spin * 2. ).is_integer() :
                raise ValueError( "Spin '%g' is not supported (must be integer or half-integer)" % spin )
            self.value = spin
    
    def __str__(self):
        """ return J as string: '4' or '3/2' for example """

        if type( self.value ) == type( 1 ):
            rets = '%i' % self.value
        elif self.value is None: rets = '?'
        elif self.value.is_integer():
            rets = '%i' % self.value
        else: rets = '%i/2' % (2*self.value)
        return rets

    def __float__(self): return self.value

    def __copy__( self ) :

        return( spin( self.value ) )

    __deepcopy__ = __copy__

class parity:
    """Store the parity for a nuclear level. Allowed values for the parity are '+','-' and '?' (for unknown)."""

    def __init__(self, parity):
        if type(parity) == str:
            try: self.value = {'+':1, '-':-1, '?':None}[parity]
            except KeyError:
                raise ValueError( "Can't understand parity: %s" % parity )
        else:
            self.value = parity

    def __str__(self):

        return {1:'+', -1:'-', None:'?'}[ self.value ]

    def __copy__( self ) :

        return( parity( self.value ) )

    __deepcopy__ = __copy__

class undefinedLevel :

    def __init__( self, level ) :

        self.level = level

    def __str__( self ) :

        return( 'u:%s' % str( self.level.toString( keepPeriod = False ) ) )

    def getValueAs( self, unit ) :

        return( self.level.getValueAs( unit ) )

    def toString( self, keepPeriod = False ) :

        return( str( self ) )

class element( base ) :

    moniker = 'element'

    def getZ_A_SuffixAndZA( self ) :

        return( nuclear.elementZFromSymbol( self.name.split('{')[0] ), 0, '', 0 )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xmlString = [ '%s<%s id="%s">' % ( indent, self.moniker, self.name ) ]
        ACs = sorted( self.configurations.keys( ) )
        if( len( ACs ) > 0 ) :
            xmlString.append( '%s<configurations>' % indent2 )
            for AC in ACs : xmlString += self.configurations[AC].toXMLList( indent3, **kwargs )
            xmlString[-1] += '</configurations>'
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class atomicConfiguration :

    moniker = 'atomicConfiguration'

    def __init__( self, subshell, bindingEnergy, electronNumber ) :

        self.subshell = subshell
        self.bindingEnergy = bindingEnergy
        self.electronNumber = electronNumber
        self.decays = []

    def __copy__( self ) :

        AC = atomicConfiguration( self.subshell, self.bindingEnergy, self.electronNumber )
        for decay in self.decays : AC.addDecay( decay )
        return( AC )

    def addDecay( self, decay ) :

        self.decays.append( decay.__copy__( ) )
        
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        ender = ''
        if( len( self.decays ) == 0 ) : ender = '/'
        xmlString = [ '%s<%s subshell="%s" bindingEnergy="%s" electronNumber="%s"%s>' % 
            ( indent, self.moniker, self.subshell, self.bindingEnergy, self.electronNumber, ender ) ]
        for decay in self.decays : xmlString += decay.toXMLList( indent2, **kwargs )
        if( len( ender ) == 0 ) : xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class atomicDecay :

    moniker = 'decay'

    def __init__( self, probability, daughters = None ) :

        self.probability = probability
        self.daughters = daughters or []

    def __copy__( self ) :

        decay = atomicDecay( self.probability, self.daughters )
        return( decay )

    __deepcopy__ = __copy__

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        ender = ''
        if( len( self.daughters ) == 0 ) : ender = '/'
        xmlString = [ '%s<%s probability="%s"%s>' % ( indent, self.moniker, self.probability, ender ) ]
        for daughter in self.daughters : xmlString.append( '%s<daughter particle="%s"/>' % ( indent2, daughter ) )
        if( len( ender ) == 0 ) : xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

addFamily( photon )
addFamily( lepton )
addFamily( isotope )
addFamily( nuclearLevel )
addFamily( FissionProduct, override = True )
addFamily( thermalNeutronScatteringLawIsotope, override = True )
