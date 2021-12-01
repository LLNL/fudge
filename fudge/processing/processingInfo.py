# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

conserveParticle = 'conserveParticle'
conserveEnergy = 'conserveEnergy'
conserveParticleAndEnergy = 'conserveParticleAndEnergy'

from PoPs import IDs as IDsPoPsModule

from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

class tempInfo :

    def __init__( self ) :

        self.dict = {}

    def __contains__( self, key ) :

        return( key in self.dict )

    def __delitem__( self, key ) :

        del self.dict[key]

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) : 

        self.dict[key] = value

    def get( self, key, default=None ):

        return self.dict.get( key, default )

    def keys( self ) :

        return( list( self.dict.keys( ) ) )

class processInfo2 :

    def __init__( self, style, reactionSuite, flux = None, logFile = None, verbosity = 0, verbosityIndentStep = '  ',
            energyAccuracy = 1e-6, momentumAccuracy = 1e-3 ) :

        if( not( isinstance( style, stylesModule.style ) ) ) : raise TypeError( 'invalid style' )
        self.__style = style

        if( not( isinstance( reactionSuite, reactionSuiteModule.reactionSuite ) ) ) : raise TypeError( 'invalid reactionSuite' )
        self.__reactionSuite = reactionSuite

        self.__flux = flux
        self.__logFile = logFile
        self.__verbosity = int( verbosity )

        if( not( isinstance( verbosityIndentStep, str ) ) ) : raise TypeError( 'invalid verbosityIndentStep' )
        if( len( verbosityIndentStep ) != verbosityIndentStep.count( ' ' ) ) :
            raise ValueError( 'verbosityIndentStep must only contain the space character' )
        self.__verbosityIndentStep = verbosityIndentStep

        self.__energyAccuracy = energyAccuracy
        self.__momentumAccuracy = momentumAccuracy

    @property
    def energyAccuracy( self ) :

        return( self.__energyAccuracy )

    @property
    def momentumAccuracy( self ) :

        return( self.__momentumAccuracy )

    @property
    def reactionSuite( self ) :

        return( self.__reactionSuite )

    @property
    def style( self ) :

        return( self.__style )

    @property
    def verbosityIndentStep( self ) :

        return( self.__verbosityIndentStep )

    @property
    def flux( self ) :

        return( self.__flux )

    @property
    def logFile( self ) :

        return( self.__logFile )

    @property
    def verbosity( self ) :

        return( self.__verbosity )

    def checkStyle( self, _class ) :

        if( not( issubclass( _class, stylesModule.style ) ) ) : raise TypeError( 'class not a sub-class of style' )
        if( isinstance( self.style, _class ) ) : return
        raise TypeError( 'style %s not of style %s' % ( _class.moniker, self.style.moniker ) )

class processInfo :

    def __init__( self, target, particles = None, flux = None, logFile = None, verbosity = 0 ) :

        self.target = target
        self.dict = {}
        self['particles'] = {}
        if( not( particles is None ) ) :
            for name in particles : self.addParticle( particles[name] )
        if( not flux is None ) : self['flux'] = flux
        self['verbosityIndent'] = ''
        self['verbosity'] = verbosity
        self['logFile'] = logFile

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) :

        self.dict[key] = value

    def addParticle( self, particle ) :

        self['particles'][particle.id]  = particle

    def getParticleGroups( self, name ) :

        return( self['particles'][name].groups )

    def getParticle_lMax( self, name ) :

        return( self['particles'][name].lMax )

    def getProjectileName( self ) :

        return( self.target.projectile.pid )

    def getTargetName( self ) :

        return( self.target.target.pid )

    def isProcessParticle( self, name ) :

        return( name in self['particles'] )

    def process( self, verbosity = None, verbosityIndent = None ) :

        verbositySave = self['verbosity']
        if( not ( verbosity is None ) ) : self['verbosity'] = verbosity
        if( verbosityIndent is None ) : verbosityIndent = self['verbosityIndent']
        try :
            self.target.process( self, verbosityIndent )
        finally:
            self['verbosity'] = verbositySave

class processInfoParticle :

    def __init__( self, name, groups, lMax, conservationFlag = conserveParticle ) :

        self.pid = name
        self.groups = groups
        self.lMax = lMax
        self.conservationFlag = conservationFlag

    def __repr__( self ) :

        s = '\n%s %s %s\n' % ( self.pid, self.lMax, self.conservationFlag )
        s += str( self.groups )
        return( s )

class processInfoLLNL( processInfo ) :

    def __init__( self, target, groups = None, flux = None, LLNL_Pn = True, lMax = 3, logFile = None, verbosity = 0 ) :

        particles = {}
        for particle in groups :
            lMax_ = 0
            if( particle == target.projectile.pid ) : lMax_ = lMax
            conservationFlag = conserveParticleAndEnergy
            if( particle == 'n' ) : conservationFlag = conserveParticle
            if( particle == IDsPoPsModule.photon ) : conservationFlag = conserveEnergy
            particles[particle] = processInfoParticle( particle, groups[particle], lMax_, conservationFlag )
        processInfo.__init__( self, target, particles, flux = flux, logFile = logFile, verbosity = verbosity )
        self['workDir'] = 'xndfgen.work'
        styles = []
        if( LLNL_Pn ) : styles.append( 'LLNL_Pn' )
        self['styles'] = styles
