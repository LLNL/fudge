# <<BEGIN-copyright>>
# <<END-copyright>>

__metaclass__ = type

conserveParticle = 'conserveParticle'
conserveEnergy = 'conserveEnergy'
conserveParticleAndEnergy = 'conserveParticleAndEnergy'

class tempInfo :

    def __init__( self ) :

        self.dict = {}

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) : 

        self.dict[key] = value

    def get( self, key, default=None ):

        return self.dict.get( key, default )

    def keys( self ) :

        return( self.dict.keys( ) )

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

        self['particles'][particle.name]  = particle

    def getParticleGroups( self, name ) :

        return( self['particles'][name].groups )

    def getParticle_lMax( self, name ) :

        return( self['particles'][name].lMax )

    def getProjectileName( self ) :

        return( self.target.projectile.getToken( ) )

    def getTargetName( self ) :

        return( self.target.target.getToken( ) )

    def isProcessParticle( self, name ) :

        return( name in self['particles'] )

    def process( self, verbosity = None, verbosityIndent = None ) :

        verbositySave = self['verbosity']
        if( not ( verbosity is None ) ) : self['verbosity'] = verbosity
        if( verbosityIndent is None ) : verbosityIndent = self['verbosityIndent']
        doRaise = False
        try :
            self.target.process( self, verbosityIndent )
        except :
            doRaise = True
        self['verbosity'] = verbositySave
        if( doRaise ) : raise

class processInfoParticle :

    def __init__( self, name, groups, lMax, conservationFlag = conserveParticle ) :

        self.name = name
        self.groups = groups
        self.lMax = lMax
        self.conservationFlag = conservationFlag

    def __repr__( self ) :

        s = '\n%s %s %s\n' % ( self.name, self.lMax, self.conservationFlag )
        s += `self.groups`
        return( s )

class processInfoLLNL( processInfo ) :

    def __init__( self, target, groups = None, flux = None, LLNL_MC = True, LLNL_Pn = True, lMax = 3, logFile = None, verbosity = 0 ) :

        particles = {}
        for particle in groups :
            lMax_ = 0
            if( particle == target.projectile.getToken( ) ) : lMax_ = lMax
            conservationFlag = conserveParticleAndEnergy
            if( particle == 'n' ) : conservationFlag = conserveParticle
            if( particle == 'gamma' ) : conservationFlag = conserveEnergy
            particles[particle] = processInfoParticle( particle, groups[particle], lMax_, conservationFlag )
        processInfo.__init__( self, target, particles, flux = flux, logFile = logFile, verbosity = verbosity )
        self['workDir'] = 'xndfgen.work'
        styles = []
        if( LLNL_MC ) : styles.append( 'LLNL_MC' )
        if( LLNL_Pn ) : styles.append( 'LLNL_Pn' )
        self['styles'] = styles
