# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the xParticle class.

Some stuff done in the product class will be put here.
"""

import re
from fudge.core.ancestry import ancestry
from pqu import physicalQuantityWithUncertainty
from fudge.core.utilities import fudgeZA, reactionStrings
from fudge.legacy.converting import endfFormats
from fudge.gnd import warning

__metaclass__ = type

particleType_Unknown = 'unknown'
particleType_Nucleon = 'nucleon'
particleType_Nuclear = 'nucleus'
particleType_Photon = 'photon'
particleType_Lepton = 'lepton'

monikerParticle = 'particle'
monikerLevel = 'level'

class baseParticle( ancestry ):
    """ abstract base class for particles and excited levels.
    Every particle appearing in the gnd particleList should be an instance of baseParticle """
    def __init__( self, token, attribute=None ):
        ancestry.__init__( self, token, None, attribute=attribute )

class xParticle(baseParticle) :
    """This class contains information about a particle (e.g., mass, charge, spin)."""

    def __init__( self, name, genre, mass = None, attributes = {} ) :

        baseParticle.__init__( self, monikerParticle, attribute = 'name' )
        self.name = name
        self.genre = genre
        self.attributes = attributes or {}
        self.setMass( mass )
        self.levels = {}

    def __getitem__( self, key ) :

        return( self.levels[key] )

    def __eq__( self, other ) :

        return( self.localCompare( other ) == 0 )

    def __ge__( self, other ) :

        return( self.localCompare( other ) >= 0 )

    def __gt__( self, other ) :

        return( self.localCompare( other ) > 0 )

    def __le__( self, other ) :

        return( self.localCompare( other ) <= 0 )

    def __lt__( self, other ) :

        return( self.localCompare( other ) < 0 )

    def __ne__( self, other ) :

        return( self.localCompare( other ) != 0 )

    def localCompare( self, other ) :

        Z, A, suffixSelf,  ZAMSelf  = self.getZ_A_SuffixAndZA()
        Z, A, suffixOther, ZAMOther = other.getZ_A_SuffixAndZA()
        if( ZAMSelf < ZAMOther ) : return( -1 )
        if( ZAMSelf > ZAMOther ) : return(  1 )
        if( suffixSelf < suffixOther ) : return( -1 )
        if( suffixSelf > suffixOther ) : return(  1 )
        return( 0 )

    def __str__( self ) :
        return( self.name )

    def getName( self ) :
        return( self.name )

    def getToken( self ) :
        return( self.name )

    def getGenre( self ) :
        return( self.genre )

    def getMass( self, unit ) :
        """Returns the mass of the particle if possible, otherwise None is returned."""
        if( self.attributes['mass'] is None ) : return( None )
        return( self.attributes['mass'].getValueAs( unit ) )

    def setMass( self, mass ) :
        """Sets mass of the particle. Mass must be a physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty instance or None."""
        if( not( mass is None ) ) :
            if( not isinstance( mass, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ):
                raise Exception( 'Mass must be a sub-class of PhysicalQuantityWithUncertainty' )
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
    
    # three functions for compatibility with excited states:
    def getLevelIndex( self ) : return 0

    def getLevelAsFloat( self, unit, default = 0. ): return 0

    def getMetaStableIndex( self ): return 0

    def addLevel( self, level ):

        if isinstance(level, nuclearLevel):
            level.groundState = self
            self.levels[ level.label ] = level
            level.setParent( self )
        else: raise TypeError( "Energy levels should be of class nuclearLevel!" )

    def getZ_A_SuffixAndZA( self ) :

        Z, A, suffix, ZA = fudgeZA.gndNameToZ_A_Suffix( self.name )
        return( Z, A, suffix, 1000 * Z + A )

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
                warnings.append( warning.discreteLevelsOutOfOrder( level ) )
            else:
                emax = enow

        return warnings
    
    def toXMLList( self, indent = '' ) :
        qs = ''
        # order attributes:
        attrs = self.attributes.copy()
        for attr in ('mass','spin','parity'):
            if attr in attrs: qs += ' %s="%s"' % (attr, attrs.pop(attr))
        for q in attrs:
            if type(attrs[q]) is bool:
                qs += ' %s="%s"' % (q, str(attrs[q]).lower())
            else:
                qs += ' %s="%s"' % ( q, attrs[q] )
        xml = [ '%s<particle name="%s" genre="%s"%s' % ( indent, self.name, self.genre, qs ) ]
        if self.levels:
            xml[-1] += '>'
            for level in sorted(self.levels.keys()):
                xml.extend( self.levels[level].toXMLList(indent+'  ') )
            xml[-1] += '</particle>'
        else: xml[-1] += '/>'
        return xml

class nuclearLevel(baseParticle):

    def __init__( self, name, energy, label, gammas=[], attributes={}, groundState=None):

        baseParticle.__init__( self, monikerLevel, attribute = 'name' )
        self.name = name
        self.energy = energy
        if( type( label ) != type( 1 ) ) :
            if( label not in [ 'c', 's' ] ) : raise Exception( 'invalid level = "%s" for particle "%s"' % ( label, name ) )
        self.label = label      # could be integer, 'c' for continuum or 's' for sum.
        self.gammas = gammas or []
        self.attributes = attributes or {}
        self.groundState = groundState  # points to xParticle instance

    def __lt__(self, other):
        if isinstance(other, nuclearLevel):
            if self.groundState == other.groundState: return self.label < other.label
            return self.groundState < other.groundState
        elif isinstance(other, xParticle):
            return self.groundState < other
        else:
            raise Exception("Can't compare nuclear level to %s instance" % other.__class__)

    def __str__( self ) :

        return( self.name )

    def addGamma( self, gamma ) :

        if not isinstance( gamma, nuclearLevelGamma ):
            raise TypeError( 'gamma must be xParticle.nuclearLevelGamma instance' )
        self.gammas.append( gamma )
        gamma.setParent( self )

    def getName( self ): return self.name

    def getToken( self ): return self.name

    def getMass( self, unit ):
        """Get the mass (as a float) in the specified unit.
        Result includes ground-state mass + level excitation energy."""

        return self.groundState.getMass( unit ) + self.energy.getValueAs( unit+'*c**2' )
    
    def getLevelIndex( self ) :
        return self.label

    def getMetaStableIndex( self ) :
        if( 'metaStableIndex' in self.attributes ) : return( int( self.attributes['metaStableIndex'] ) )
        return( 0 )

    def getLevelAsFloat( self, unit, default = 0. ) :
        try :
            level = float( self.energy.getValueAs( unit ) )
        except :    # for continuum and sum, self.energy = None
            if( default is None ) :
                raise Exception( 'Could not convert level energy = "%s" to float' % self.energy )
            level = default
        return( level )
    
    def getZ_A_SuffixAndZA( self ):
        Z, A, suffix, ZA = fudgeZA.gndNameToZ_A_Suffix( self.name )
        return (Z, A, suffix, 1000 * Z + A )

    def getSpin( self ): return self.attributes.get('spin')

    def getParity( self ): return self.attributes.get('parity')

    def check( self, info ):
        '''
        Check parts of nuclearLevel
        
        Options & defaults:
            'gammaBRSumAbsTol'      1e-6
        '''
        warnings = []
        gammaBRSumAbsTol = info.get( 'gammaBRSumAbsTol', 1e-6 )
        if self.gammas and abs( sum([gamma.probability for gamma in self.gammas]) - 1.0 ) > gammaBRSumAbsTol:
            warnings.append( warning.unnormalizedGammas( sum([gamma.probability for gamma in self.gammas]), self ) )
        return warnings

    def toXMLList( self, indent = '' ):
        qs = ''
        # order attributes:
        attrs = self.attributes.copy()
        for attr in ('energy','spin','parity'):
            if attr in attrs: qs += ' %s="%s"' % (attr, attrs.pop(attr))
        for q in attrs: qs += ' %s="%s"' % (q, attrs[q])
        xml = [ '%s<level name="%s" label="%s" energy="%s"%s' % (indent, self.name, self.label, self.energy, qs) ]
        if self.gammas:
            xml[-1] += '>'
            for gamma in self.gammas:
                xml.extend( gamma.toXMLList(indent+'  ') )
            xml[-1] += '</level>'
        else: xml[-1] += '/>'
        return xml
    
    def toENDF6( self, baseMT, endfMFList, flags, tempInfo ) :

        MF, MT, LP = 12, baseMT + int( self.label ), 0
        nGammas, gammaData, levelEnergy_eV = len( self.gammas ), [], self.energy.getValueAs( 'eV' )
        for gamma in self.gammas : gammaData.append( gamma.toENDF6List( ) )
        LGp = len( gammaData[0] )
        for gamma in gammaData : gamma[0] = levelEnergy_eV - gamma[0]
        endfMFList[MF][MT] = [ endfFormats.endfHeadLine( tempInfo['ZA'], tempInfo['mass'], 2, LGp - 1, MT - baseMT, 0 ),
            endfFormats.endfHeadLine( levelEnergy_eV, 0., LP, 0, LGp * nGammas, nGammas ) ]
        gammaData.sort( reverse = True )
        endfMFList[MF][MT] += endfFormats.endfNdDataList( gammaData )
        endfMFList[MF][MT].append( endfFormats.endfSENDLineNumber( ) )

        # Currently, assume all distributions are isotropic
        endfMFList[14][MT] = [ endfFormats.endfHeadLine( tempInfo['ZA'], tempInfo['mass'], 1, 0, nGammas, 0 ) ]
        endfMFList[14][MT].append( endfFormats.endfSENDLineNumber( ) )

class nuclearLevelGamma( ancestry ) :

    def __init__( self, finalLevel, angularDistribution, probability, nonRadiativeProbability = 0.0 ) :

        ancestry.__init__( self, 'gamma', parent=None, attribute='finalLevel' )
        self.finalLevel = finalLevel
        self.probability = probability
        self.nonRadiativeProbability = nonRadiativeProbability
        self.angularDistribution = angularDistribution

    def __str__( self ) :

        return( 'Gamma energy = %s, final level = %s, gammaEmissionProbability = %s, internalConversionProbability = %s' % \
               ( self.getEnergy(), self.finalLevel, self.probability, self.nonRadiativeProbability ) )

    def getEnergy( self, unit='eV' ) :
        """Calculate emitted gamma energy (including recoil correction) in specified units."""

        deltaE = self.getParent().getLevelAsFloat( unit ) - self.finalLevel.getLevelAsFloat( unit )
        correction = 0.5 * deltaE**2 / self.finalLevel.getMass( unit+'/c**2' )
        return deltaE - correction  # would be '+' for absorption

    def toENDF6List( self ) :

        # ENDF ignores recoil correction:
        energy = self.getParent().getLevelAsFloat('eV') - self.finalLevel.getLevelAsFloat('eV')
        return( [ energy, self.probability ] )
    
    def toXMLList( self, indent = '' ) :
        xml =  '%s<gamma finalLevel="%s" probability="%s"' % (indent, self.finalLevel, self.probability)
        if self.nonRadiativeProbability: xml += ' nonRadiativeProbability="%s"' % self.nonRadiativeProbability
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
                raise ValueError, "Can't decipher spin from string: %s" % spin
            num,den = match.groups()
            #self.spin = float(num)/(2 if den else 1) # this is such a pretty construction that I hate to remove it for Python 2.4
            if den: self.value = float(num)/(2)
            else:   self.value = float(num)/(1)
        else:
            if not (spin*2).is_integer():
                raise ValueError, "Spin '%g' is not supported (must be integer or half-integer)" % spin
            self.value = spin
    
    def __str__(self):
        """ return J as string: '4' or '3/2' for example """
        if type( self.value ) == type( 1 ):
            rets = '%i' % self.value
        elif self.value == None: rets = '?'
        elif self.value.is_integer():
            rets = '%i' % self.value
        else: rets = '%i/2' % (2*self.value)
        return rets
    
class parity:
    """Store the parity for a nuclear level. Allowed values for the parity are '+','-' and '?' (for unknown)."""

    def __init__(self, parity):
        if type(parity) == str:
            try: self.value = {'+':1, '-':-1, '?':None}[parity]
            except KeyError:
                raise ValueError, "Can't understand parity: %s" % parity
        else:
            self.value = parity

    def __str__(self):

        return {1:'+', -1:'-', None:'?'}[ self.value ]
