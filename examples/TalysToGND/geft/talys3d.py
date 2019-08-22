import math
import talysFiles, endlLikeDeltaFunc, fudge, talysStructure
from talysMisc import *
#from fudge.legacy.endl.endl2dmathClasses import *
#from fudge.legacy.endl.endl3dmathClasses import *
#import fudge.legacy.endl.endl2dmathClasses 
#import fudge.legacy.endl.endl3dmathClasses 

class Spectra :
    """
Bass class for continuum spectra creation from TALYS output files.
Provides some functions to process the raw data from TALYS.
"""
    def fixEQbug( self, SPFile ) :
        """
Trim function on 2d data leaves a 1 trailing zero after the data.
This point falls outside the energy range allowable from E' <= Ein + Q.
Drop this point if this is the case.
"""
        energy = SPFile.energy
        for label, data in SPFile.items() :
            if data.xMax() > energy + self.Q :
                if debug :
                    print '%s:%s Chopping Ep at %s for E=%s' % ( SPFile.filename, label, (energy + self.Q), SPFile.energy )
                    print 'The value of the point dropped was %s' % ( data.getValue( data.xMax() ) )
            if 'Probability' in label :
                data.set( data.slice( xMax=( energy + self.Q ) ).normalize() )
            elif 'CrossSection' in label :
                data.set( data.slice( xMax=( energy + self.Q ) ) )
            if len( data ) == 1 or data.max() == 0 :
                del SPFile[label]

    def insertThreshold( self, energy=None ) :
        """
Inserts a (fake) delta function at threshold for the energy spectra,
or if energy keyword is specified, at this energy, which is useful for filling in energies
to match I=0 energy ranges, rather that use fixQThreshold in fudge which copies across higher
energy probability which then violates the E'<=E+Q condition.
"""
        if hasattr( self, 'threshold' ) and not energy : energy = self.threshold
        if energy :
            size = math.pow(10,math.floor(math.log10(float("%.1e"%energy))-3))
            for label, data in self.items() :
                if debug : print 'Setting threshold for %s at %s' % ( label, energy )
                data.setAtX( energy, endl2dmathmisc.ramp2d( size, 2*size, down=True ).normalize() )

    def process( self ) :
        """Default list of functions for processing particle spectra."""
        self.removeEmptyData()

class ContinuumEnergySpectra( talysFiles.OnePerArgList, Spectra, talysFiles.OutputFile ) :
    """
Class for continuum particle spectra. Combines the list of 2d particle spectra found
in a different file for each incident energy, into endl3dmath objects (ie L=0 distributions)
For particles, combines the Probability distributions and endl3dmath objects
are found in particleName+EnergySpectrum.
For Gammas, the raw CrossSection is returned in GammaEnergySpectrum since the discrete
Gammas need to be added before normalizing.
"""
    def __init__( self, talysZA, particleList ) :
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            fmin, fmax = 9, 16
            self.filelist = getFileList( talysZA, 'sp'+particleList+'E*.tot' )

            try    : CrossSectionFile = talysFiles.pe_CON( talysZA, talysYTags( talysZA.yi ) + eFromParticleList( particleList ) )
            except : CrossSectionFile = talysFiles.XSnpdtha_TOT( talysZA, particleList )
            try:
                self.threshold = CrossSectionFile.threshold
                self.Q = CrossSectionFile.Q
            except AttributeError : pass

            particleNameList = map( particleName, [1,2,3,4,5,6,7] )
            for particle in particleNameList :
                label3d = particle+'EnergySpectrum'
                self[label3d] = fudge.legacy.endl.endl3dmath( data = [] )
            self['GammaCrossSection'] = fudge.legacy.endl.endl3dmath( data = [] )

            for filename in self.filelist :
                energy = filename[fmin:fmax]
                SPFile = talysFiles.SPnpdthaEenergy_TOT( talysZA, particleList, energy )
                self.fixEQbug( SPFile )
                for particle in particleNameList :
                    label3d = particle+'EnergySpectrum'
                    label2d = particle+'Probability'
                    if label2d in SPFile :
                        self[label3d].setAtX( SPFile.energy, SPFile[label2d] )
                label2d = 'GammaCrossSection'
                label3d = label2d
                if label2d in SPFile :
                    self[label3d].setAtX( SPFile.energy, SPFile[label2d] )
            self.process()
    
    def process( self ) :
        self.removeEmptyData()
        self.insertThreshold()
        
class DiscreteEnergySpectra( talysFiles.OnePerArgList, Spectra, talysFiles.OutputFile ) :
    """
Class for discrete gamma distributions.
Processes the discrete gammas found in xs*.gam into endl3dmath objects.
"""
    def __init__( self, talysZA, particleList ) :
        #print 'new DiscreteEnergySpectra'
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            gammaList = talysFiles.XSnpdtha_GAM( talysZA, particleList ).gammaList
            label = 'DiscreteGammaCrossSection'
            self[label] = fudge.legacy.endl.endl3dmath( data = [] )
            for gammas in gammaList :
                E = gammas.energy
                gammas.sort()
                EgammaP = fudge.legacy.endl.endl2dmath( data = [] )
                for gamma in gammas :
                    EgammaP += endlLikeDeltaFunc.deltafn( gamma.energy, h=gamma.CrossSection )
                if EgammaP :
                    self[label].setAtX( E, EgammaP )
            self.process()
    
class DiscreteLevelEnergySpectra( talysFiles.OnePerArgList, Spectra, talysFiles.OutputFile ) :
    """
Class for discrete gamma distributions per level.
Processes the discrete gammas found in xs*.gam into endl3dmath objects for each level.
"""
    def __init__( self, talysZA, particleList, level ) :
        #print 'new DiscreteLevelEnergySpectra'
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, particleList, level )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            pe = talysYTags( talysZA.yi ) + eFromParticleList( particleList )
            gammaList = talysFiles.XSnpdtha_GAM( talysZA, particleList ).gammaList
            #print 'gammaList from XSnpdtha_GAM: ',type(gammaList)
            levelScheme = talysStructure.levels().levelScheme( talysZA.Z, talysZA.A )
            label = 'DiscreteGammaEnergySpectrum'
            discreteFile = BinaryDiscreteList( talysZA, pe ).getDiscreteForLevel( level )
            self.X1 = discreteFile.X1
            self.threshold = discreteFile.threshold
            self[label] = fudge.legacy.endl.endl3dmath( data = [] )
            if len(gammaList)==0:
                gammas = []
                return
            else:
                gammas = gammaList[-1]
                E = gammas.energy
                EgammaP = fudge.legacy.endl.endl2dmath( data = [] )
                EgammaP = self.getGammaCascade( gammas, level, levelScheme, EgammaP )
                if EgammaP :
                    self[label].setAtX( self.threshold, EgammaP.normalize() )
                    self[label].setAtX( E, EgammaP.normalize() )
            self.process()
    
    def getGammaCascade( self, gammas, level_i, levelScheme, EgammaP, BRi=1.0 ) :
        for gamma in filter( lambda gamma: gamma.level_i == level_i, gammas ) :
            #print levelScheme.levels[gamma.level_i].branches
                    
            try:
                BR = BRi*levelScheme.levels[gamma.level_i].branches[gamma.level_f].BR
                #print gamma.level_i,gamma.level_f
            #print BR, levelScheme.levels[gamma.level_i].branches[gamma.level_f].BR
            except:
                #print 'OUT OF RANGE: ',gamma.level_i,gamma.level_f
                #print levelScheme.levels[gamma.level_i].branches
                ##print levelScheme.levels[gamma.level_i].branches[gamma.level_f]
                #k = levelScheme.levels[gamma.level_i].branches.keys()
                #print k
                #print levelScheme.levels[gamma.level_i].branches[k[0]]
                #print levelScheme.levels[gamma.level_i].branches[k[0]].BR
                BR = 0.
            EgammaP += endlLikeDeltaFunc.deltafn( gamma.energy, h=BR )
            EgammaP = self.getGammaCascade( gammas, gamma.level_f, levelScheme, EgammaP, BRi=BR )
        return EgammaP
    
    def process( self ) :
        self.removeEmptyData()

class S0GammaEnergySpectra( talysFiles.OnePerArgList, Spectra, talysFiles.OutputFile ) :
    """
Calculates the S=0 continuum gamma energy spectra, which includes discrete gammas which come from cascases originating in the continuum .
The S=1 discrete gammas have to be subtracted from the total discrete gammas given in the TALYS .gam output files.
"""
    def __init__( self, talys, particleList ) :
        print 'new S0GammaEnergySpectra'
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            ContinuumGammaCrossSection = ContinuumEnergySpectra( talysZA, particleList ).get( 'GammaCrossSection' )
            gammaList = talysFiles.XSnpdtha_GAM( talysZA, particleList ).gammaList
            DiscreteGammaCrossSection = fudge.legacy.endl.endl3dmath( data = [] )
            levelScheme = talysStructure.levels().levelScheme( talysZA.Z, talysZA.A )
            for gammas in gammaList :
                E = gammas.energy
                gammas.sort()
                EgammaP = fudge.legacy.endl.endl2dmath( data = [] )
                for gamma in gammas :
                    S1GammaCrossSection = self.getS1GammaCrossSection( levelScheme, gamma.level_i, gamma.level_f )
                    EgammaP += endlLikeDeltaFunc.deltafn( gamma.energy, h=gamma.CrossSection-S1GammaCrossSection )
                if EgammaP :
                    DiscreteGammaCrossSection.setAtX( E, EgammaP )
            
            self['GammaEnergySpectrum'] = add3d( ContinuumGammaCrossSection, DiscreteGammaCrossSection )
            self.process()

class GammaEnergySpectra( talysFiles.OnePerArgList, Spectra, talysFiles.OutputFile ) :
    """
Calculates the total Gamma energy spectra by adding the continuum gamma and discrete gamma
energy spectra to produce an endl3dmath object of the sum.
"""
    def __init__( self, talysZA, particleList ) :
        def add3d( a, b ) :
            if not a and not b : return
            elif not a : a = fudge.legacy.endl.endl3dmath( data = [] ) #still need to run through to normalize data at end
            elif not b : b = fudge.legacy.endl.endl3dmath( data = [] ) #still need to run through to normalize data at end
            def union1d( a, b ) :
                d = a.copyData()
                for x in b.data :
                    if x not in d.data : d.data.append( x )
                d.data.sort()
                return d
            d = fudge.legacy.endl.endl3dmath( data = [] )
            x_ = union1d( a.xArray(), b.xArray() )
            for x in x_ :
                yz_a = None
                yz_b = None
                for xyz1 in a.data :
                    if x == xyz1[0] :
                        yz_a = fudge.legacy.endl.endl2dmath( xyz1[1] )
                        break
                for xyz2 in b.data :
                    if x == xyz2[0] :
                        yz_b = fudge.legacy.endl.endl2dmath( xyz2[1] )
                        break
                if not yz_b :
                    d.setAtX( x, normalize( yz_a ) )
                else :
                    if not yz_a :
                        d.setAtX( x, normalize( yz_b ) )
                    else :
                        yz_a += yz_b
                        d.setAtX( x, normalize( yz_a ) )
            return d
        def normalize(xy_):
            xy = fudge.legacy.endl.endl2dmath()
            for x,y in xy_.data:
                x = float(fudge.legacy.endl.endl4dmathmisc.endl4d_repr_yFormat%x)
                y = float(fudge.legacy.endl.endl4dmathmisc.endl4d_repr_zFormat%x)
                xy.setValue(x,y)
            return xy.normalize()
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            self['GammaEnergySpectrum'] = \
                add3d( ContinuumEnergySpectra( talysZA, particleList ).get( 'GammaCrossSection' ), \
                        DiscreteEnergySpectra( talysZA, particleList ).get( 'DiscreteGammaCrossSection' ) )
            self.process()
    
class PreEquilibriumRatio( talysFiles.OnePerArgList, talysFiles.OutputFile ) :
    """
Class for pre-equilibrium ratios for particle spectra used in Kalbach systematics.
The 2d data for each file from eSPECenergy_TOT classes is processed for each energy
and put into an endl3dmath object obtained by the key 'PreEquilibriumRatio'.
"""
    def __init__( self, talysZA, ejectile ) :
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, ejectile )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            fmin, fmax = 5, 12
            self.ejectile = ejectile
            self.filelist = getFileList( talysZA, ejectile+'spec*.tot' )
            label = 'PreEquilibriumRatio'
            self[label] = fudge.legacy.endl.endl3dmath( data = [] )
            for filename in self.filelist :
                energy = filename[fmin:fmax]
                SPECFile = talysFiles.eSPECenergy_TOT( talysZA, ejectile, energy )
                if label in SPECFile: self[label].setAtX( SPECFile.energy, SPECFile[label] )
            self.process()

    def process( self ) :
        """Default list of functions for processing particle spectra."""
        self.removeEmptyData()

class BinaryAngles( talysFiles.OnePerArgList, talysFiles.OutputFile ) :
    """
Class for angular distributions for binary reaction where only one particle is emitted.
The 2d data from pe_energyANG_LXX classes is processed for each energy and put into
an endl3dmath object obtained by the key 'AngularSpectrum', this is the I=1 P(E|mu) distribution.
"""
    def __init__( self, talysZA, inelasticPrefix, level ) :
        #print 'new BinaryAngles'
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix, level )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            self.pe = inelasticPrefix
            self.level = int( level )
            fmin,fmax = 2, 9
            self.filelist = getFileList( talysZA, self.pe+'*ang.L'+str(level).zfill(2) )
            label3d = 'AngularSpectrum'
            self[label3d] = fudge.legacy.endl.endl3dmath( data = [] )
            try :
                CrossSectionFile = talysFiles.pe_LXX( talysZA, self.pe, level )
                self.threshold = CrossSectionFile.threshold
                self.Q = CrossSectionFile.Q
                discreteList = BinaryDiscreteList( talysZA, self.pe )
                self.X1 = discreteList.getDiscreteForLevel( level ).X1
            except :
                self.threshold = talysFiles.TOTAL_TOT( talysZA )['Elastic'].xMin()

            for filename in self.filelist :
                energy = filename[fmin:fmax]
                ANGFile = talysFiles.pe_energyANG_LXX( talysZA, self.pe, energy, level )
                label2d = 'Probability'
                if label2d in ANGFile :
                    self[label3d].setAtX( ANGFile.energy, ANGFile[label2d] )

            self.process()

    def process( self ) :
        self.removeEmptyData()
        self.insertThreshold()
        
    def insertThreshold( self ) :
        """Inserts a flat angular distribution at threshold or at Emin."""
        if hasattr( self, 'threshold' ) : 
            if self.threshold :
                for label, data in self.items() :
                    if debug : print 'Setting threshold for %s at %s' % ( label, self.threshold )
                    data.setAtX( self.threshold, [[-1.,0.5],[1.,0.5]] )
            

class BinaryLegendreCoeffs( talysFiles.OnePerArgList, talysFiles.OutputFile ) :
    """
Class for angular distributions for binary reaction where only one particle is emitted.
The 2d data from pe_energyANG_LXX classes is processed for each energy and put into
an endl3dmath object obtained by the key 'AngularSpectrum', this is the I=1 P(E|mu) distribution.
"""
    def __init__( self, talysZA, inelasticPrefix, level ) :
        talysFiles.OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix, level )
        if self.newInstance :
            talysFiles.OutputFile.__init__( self, talysZA )
            self.pe = inelasticPrefix
            self.level = int( level )
            fmin,fmax = 2, 9
            self.filelist = getFileList( talysZA, self.pe+'*leg.L'+str(level).zfill(2) )
            label3d = 'LegendreSpectrum'
            self[label3d] = fudge.legacy.endl.endl3dmath( data = [] )
            ### look for the threshold
            try :
                CrossSectionFile = talysFiles.pe_LXX( talysZA, self.pe, level )
                self.threshold = CrossSectionFile.threshold
                self.Q = CrossSectionFile.Q
                discreteList = BinaryDiscreteList( talysZA, self.pe )
                self.X1 = discreteList.getDiscreteForLevel( level ).X1
            except :
                self.threshold = talysFiles.TOTAL_TOT( talysZA )['Elastic'].xMin()
            ### now load and process each file
            for filename in self.filelist :
                energy = filename[fmin:fmax]
                LEGFile = talysFiles.pe_energyLEG_LXX( talysZA, self.pe, energy, level )
                label2d = 'Probability'
                if label2d in LEGFile :
                    self[label3d].setAtX( LEGFile.energy, LEGFile[label2d] )

            self.process()

    def process( self ) :
        self.removeEmptyData()
        self.insertThreshold()
        
    def insertThreshold( self ) :
        """Inserts a flat angular distribution at threshold or at Emin."""
        if hasattr( self, 'threshold' ) : 
            if self.threshold :
                for label, data in self.items() :
                    if debug : print 'Setting threshold for %s at %s' % ( label, self.threshold )
                    data.setAtX( self.threshold, [[-1.,0.5],[1.,0.5]] )
            
class IsomerCrossSectionList( list ) :
    """
Class contains a list of cross section to ground and isomeric states.
Provides easy access and grouping for a set of talysFiles.XSnpdtha_LXX classes.
"""
    def __init__( self, talysZA, particleList ) :
        list.__init__( self )
        self.npdtha = particleList
        self.filelist = getFileList( talysZA, 'xs'+self.npdtha+'.L*' )
        for filename in self.filelist :
            level = int( filename.split('.')[-1][1:3] )
            self.append( talysFiles.XSnpdtha_LXX( talysZA, self.npdtha, level ) )

    def __repr__( self ) :
        return '\n'.join( self.filelist )

class BinaryDiscreteList( list ) :
    """
Class contains a list of discrete level cross sections for a particular outgoing particle.
Provides easy access and grouping for a set of talysFiles.pe_LXX classes.
"""
    def __init__( self, talysZA, inelasticPrefix ) :
        #print 'new BinaryDiscreteList'
        list.__init__( self )
        self.pe = inelasticPrefix
        self.filelist = getFileList( talysZA, self.pe+'.L*' )
        #print 'pe',self.pe,'filelist',self.filelist
        for filename in self.filelist :
            level = int( filename.split('.')[-1][1:] )
            #print 'filename',filename,'level:',level
            
            self.append( talysFiles.pe_LXX( talysZA, self.pe, level ) )
        self.process()
    
    def __repr__( self ) :
        return '\n'.join( self.filelist )

    def process( self ) :
        self.fixQIdenticalLevels()
    
    def fixQIdenticalLevels( self ) :
        """
Nudges the levels by the last dp that is printed out in endlformat
so that degenerate levels can be entered into fudge.
"""
        for l, level in enumerate( self ) :
            #if level.level == 0 : continue # level=0 may be an inelastic level
            levelp = self[l-1]
            if levelp < 0 : continue
            if level.Q == levelp.Q :
                if debug : print 'Identical levels found for %s at X1 = %s Q = %s' % ( level, level.Q, level.X1 )
                nudge = endlLikeDeltaFunc.nudgex( "%12.5e", level.Q )
                level.Q = level.Q + nudge
                level.X1 = level.X1 - nudge
                if debug : print 'Nudging %s to X1 = %s Q = %s' % ( level, level.Q, level.X1 )
                levelp.Q = levelp.Q - nudge
                levelp.X1 = levelp.X1 + nudge
                if debug : print 'Nudging %s to X1 = %s Q = %s' % ( levelp, levelp.Q, levelp.X1 )

    def getDiscreteForLevel( self, level ) :
        """Returns the pe_LXX instance for a given level."""
        for discreteLevel in self :
            if discreteLevel.level == level :
                return discreteLevel

class BinaryDiscreteAngleList( list ) :    
    """
Class contains a list of discrete level angular distributions for a particular outgoing particle.
Provides easy access and grouping for a set of BinaryAngles classes.
"""
    def __init__( self, talysZA, inelasticPrefix ) :
        list.__init__( self )
        self.pe = inelasticPrefix
        discreteList = BinaryDiscreteList( talysZA, self.pe )
        for level in discreteList :
            self.append( BinaryAngles( talysZA, self.pe, level.level ) )

    def __repr__( self ) :
        def fileString( i ) :
            return i.pe+'*ang.L'+str(i.level).zfill(2)
        return '\n'.join( map( fileString, self ) )

    def getDiscreteAngForLevel( self, level ) :
        """Returns the pe_LXX instance for a given level."""
        for BinaryAngles in self :
            if BinaryAngles.level == level :
                return BinaryAngles
