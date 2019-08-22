"""Modules to interact with the TALYS structure database"""

import os

def find_talys_distribution_dir(  ):
    which_result = " ".join( os.popen("which talys").readlines( ) )
    #print which_result
    if which_result.find("Command not found") == -1 and which_result.find("talys") != -1:
        ## get base and path names
        talys_dir = os.path.dirname(which_result)
        #print talys_dir
        #print os.path.dirname(talys_dir)
        talys_dist_dir = talys_dir+'/structure'
        #print talys_dist_dir
        
        
        #talys_dist_dir = " ".join( os.popen( 'strings `which talys` | grep /structure' ) )
        #talys_dist_dir = talys_dist_dir[ : talys_dist_dir.find( "structure" )+10 ]
        #talys_dist_dir = talys_dist_dir.split()[-1]
        #talys_dist_dir = talys_dist_dir.strip("?") ## don't know when or why this started appearing but seems to be coming right out of the talys grep.
        return talys_dist_dir
    else:
        raise ImportError( "Cannot find talys structure files, cannot import talysStructurePath. check for talys in your path." )

talysStructurePath = 'talys'
#print 'talys structure: ',talysStructurePath

def talysZFile( path, Z ) :
    return os.path.join( path, 'z'+str(Z).zfill(3) )

class fileLine :
    def __repr__( self ) :
        return self.line

class levelsMassLine( fileLine ) :
    def __init__( self, line ) :
        self.line = line
        if line[18:74].strip() : return
        self.Z = int( line[ 0: 4] )
        self.A = int( line[ 4: 8] )
        self.nLines = int( line[ 8:13] )
        self.nLevels = int( line[13:18] )
        self.symbol = line[74:80]
        self.levels = []

class levelsLine( fileLine ) :
    def __init__( self, line ) :
        self.line = line
        self.level     = int  ( line[ 0: 4] )
        self.energy    = float( line[ 4:15] )
        self.spin      = float( line[15:21] )
        self.parity    = float( line[24:26] )
        self.nBranches = int  ( line[26:29] )
        self.branches  = {}

class levelsBranchingLine( fileLine ) :
    def __init__( self, line ) :
        self.line = line
        self.level = int  ( line[29:32] )
        self.BR    = float( line[32:42] )

    def __repr__( self ) :
        return str( self.BR )

class talysLevelsFile( list ) :
    def __init__( self, path, Z ) :
        #print path
        talysLevelsFile = open( talysZFile( path, Z ) )
        list.__init__( self, talysLevelsFile )

    def __repr__( self ) :
        return ''.join( map( str, self ) )

class levels :
    """Module to interact with the levels folder in the TALYS structure database"""
    
    talysLevels = os.path.join( talysStructurePath, 'levels' )
    
    def __init__( self, talysPath=talysLevels ) :
        self.talysPath = talysPath
        
    def levelScheme( self, Z, A ) :
        levelsFile = talysLevelsFile( self.talysPath, Z )
        while levelsFile :
            try : ZA = levelsMassLine( levelsFile.pop( 0 ) )
            except : continue
            if ZA :
                try :
                    if ZA.Z == Z and ZA.A == A : break
                except : continue
            else : continue
        for iLevel in range( ZA.nLevels + 1 ) :
            level = levelsLine( levelsFile.pop(0) )
            if iLevel != level.level : raise Exception( 'Level number wrong', iLevel ,'!=',level.level )
            ZA.levels.append( level )
            for iBranch in range( level.nBranches ) :
                branch = levelsBranchingLine( levelsFile.pop( 0 ) )
                #branch.Egamma = level.energy - ZA[branch.level].energy
                level.branches[branch.level] = branch
        return ZA

class resonances( dict ) :
    """Module to interact with the resonance folder in the TALYS structure database"""
    
    talysResonances = os.path.join( talysStructurePath, 'resonances' )
    
    def __init__( self, talysPath=talysResonances ) :
        dict.__init__( self )
        self.talysPath = talysPath
    
    def loadZ( self, Z ) :
        resonancesFile = talysResonancesFile( self.talysPath, Z )
        if resonancesFile :
            self[Z] = ( talysResonancesFile( self.talysPath, Z ) )

    def getGamgam( self, Z, A ) :
        if Z not in self : self.loadZ( Z )
        try : return self[Z].getGamgam( A )
        except : return None

class talysResonancesLine( fileLine ) :
    def __init__( self, line ) :
        self.line = line
        self.Z = int( line[ 0: 4] )
        self.A = int( line[ 4: 8] )
        try :
            self.D0          = float( line[ 8:17] )
            self.deltaD0     = float( line[17:26] )
        except : pass
        try:
            self.S0          = float( line[26:31] )
            self.deltaS0     = float( line[31:36] )
        except : pass
        try:
            self.Gamgam      = float( line[36:45] )
            self.deltaGamgam = float( line[45:54] )
        except : pass
        self.restOfLine = line[54:80]

class talysResonancesFile( list ) :
    def __init__( self, path, Z ) :
        list.__init__( self )
        try : talysResonancesFile = open( talysZFile( path, Z ) )
        except : return None
        self += map( talysResonancesLine, talysResonancesFile )

    def __repr__( self ) :
        return ''.join( map( str, self ) )

    def getLineForA( self, A ) :
        for line in self :
            if line.A == A : return line

    def getGamgam( self, A ) :
        try :    return self.getLineForA( A ).Gamgam
        except : return None

class talysAbundance( dict ) :
    """Module to interact with the abundance folder in the TALYS structure database"""
    
    talysAbundances = os.path.join( talysStructurePath, 'abundance' )
    
    def __init__( self, talysPath=talysAbundances ) :
        dict.__init__( self )
        self.talysPath = talysPath
    
    def loadZ( self, Z ) :
        try: abundanceFile = open( talysZFile( self.talysPath, Z ) )
        except: return None
        for line in abundanceFile :
            data = line.split()
            A = int( data[1] )
            ZA = Z*1000+A
            self[ZA] = float( data[2] )/100

### Transplanted the mass utilities to the structure module
##
## Usage:
##    pathToTalysMasses = where is your talys installed? look in there in structure/masses/
##    could just get it from talysStructurePath = find_talys_distribution_dir() if your talys is installed according to the instructions
##    a = talysMasses(Z = 36, path = pathToTalysMasses)
##    yourMass = a.getMass( Z, A )
##

def talysZFile( path, Z ) :
    return os.path.join( path, 'z'+str(Z).zfill(3) )

class talysMasses :
    def __init__( self, Z=None, path = None ) :
        self.massFiles = []
        if not Z: 
            self.zrange = range( 111 ) 
        else: 
            self.zrange = [Z]
            
        if path: 
            if not os.path.exists( path ) : raise Exception( '%s does not exist' % path )
        else:
            raise Exception( 'must specify a path to find talys structure info')
            
        for z in self.zrange:
            self.massFiles.append( talysMassFile( path, z ) )

    def __getitem__( self, Z ) :
        return self.massFiles[Z]

    def setMass( self, Z, A, mass ) :
        if mass != None : self[Z].setMass( A, mass )

    def getMass( self, Z, A ) :
        return self[Z].getMass( A )

    def toFile( self ) :
        os.mkdir( newMassesPath )
        for Z in self.zrange :
            open( talysZFile( newMassesPath, Z ), mode='w' ).write( str( self[Z] ) )

class talysMassFile :
    def __init__( self, path, Z ) :
        talysMassFile = open( talysZFile( path, Z ) )
        self.FileLines = []
        for line in talysMassFile :
            self.FileLines.append( talysMassLine( line ) )
    
    def __getitem__( self, i ) :
        return self.FileLines[i]
    
    def __repr__( self ) :
        return ''.join( map( str, self ) )

    def getLineForA( self, A ) :
        for line in self :
            if line.A == A : return line

    def setMass( self, A, mass ) :
        if mass == None or self.getLineForA( A ) == None : return
        self.getLineForA( A ).setMass( mass )

    def getMass( self, A ) :
        if self.getLineForA( A ) == None : return None
        if self.getLineForA( A ).hasExperimentalMass == True :
            return self.getLineForA( A ).mass
        else :
            return None

class talysMassLine :
    def __init__( self, line ) :
        self.Z      = int(   line[ 0: 4] )
        self.A      = int(   line[ 4: 8] )
        self.Moller =        line[20:32]
        self.restOfLine =    line[44:80]
        excess = line[32:44] # look at excess instead of mass because of special case for 12C as it has a mass but is not experimental
        if excess != '            ' :
            self.hasExperimentalMass = True
            self.mass   = float( line[ 8:20] )
            self.excess = float( line[32:44] )
        else :
            self.hasExperimentalMass = False
            self.line = line[0:80]

    def __repr__( self ) :
        if self.hasExperimentalMass :
            line = "%4i%4i%12.6f%s%12.6f%s\n" % ( self.Z, self.A, self.mass, self.Moller, self.excess, self.restOfLine )
        else :
            line = "%s\n" % ( self.line )
        return line

    def setMass( self, mass ) :
        
        amu = 931.49386

        if mass == None : return
        if hasattr( self, 'mass' ) :
            self.excess = mass - self.mass + self.excess
            self.mass = mass
        else :
            self.mass = mass
            self.excess = ( mass - self.A ) * amu
            self.hasExperimentalMass = True
