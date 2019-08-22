"""
Module for interacting with TALYS to create evaluations.
Can run TALYS creating inputs and extracting the data from the output files and exporting evaluation into ENDL format.

"""

import fudge, os, subprocess, shutil, tarfile
import talysGndData
import talysInput,talysStructure
from talysMisc import *
from talysGndMisc import *
from fudge.structure import xensl
from fudge.particles.nuclear import *

from fudge.structure import masses
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as pqu
from fudge.core.utilities import fudgeZA

from fudge.evaluationBuilder import *

class talysZA :
    """Proxy for endlZA with a few things added to run TALYS"""
    def __init__( self, ZA, yi, path='.', talysDir=None, eLevel=None ) :
        
        self.endlZA = endlZA(ZA,yi)
        self.CList = CList
        self.reactionList = reactionList
        self.binaryCList = getBinaryCList( self.yi )
        self.cwd = os.getcwd()
        self.workDir = path
        self.runGrp = False
        if talysDir :
            self.talysDir = talysDir
        else :
            self.talysDir = os.path.join( 'talys', endlmisc.incidentParticleTags( self.yi )[1], self.sZA )
        self.calculationSuccessful = False
        print 'cwd = %s\nworkDir = %s\ntalysDir = %s' % ( self.cwd, self.workDir, self.talysDir )
        self.endlZA.workDir = self.workDir
        self.endlZA.workDir = os.path.join( self.cwd, 'endl_output')
        print 'endlZA.workDir = %s' % ( self.endlZA.workDir)
        if not os.path.exists( self.endlZA.workDir ) :
            os.makedirs( self.endlZA.workDir )
        else:
            print 'Using existing output directory: %s' % self.endlZA.workDir 

        print self.endlZA.suffix
        if self.endlZA.suffix :
            if eLevel : self.targetELevel = eLevel
            else : self.targetELevel = xensl.Levels( endlZA.ZA )( endlZA.suffix ).energy
        else: self.targetELevel = 0
        
        # initialize the reactionSuite:
        print 'incident name: ',incidentnames[yi]
        self.pg = particleLibrary( database=None )
        self.rsb = reactionSuiteBuilder( self.pg, projectile=GNDincidentnames[yi], target=nucleusNameFromZA(ZA), library="talys evaluation", documentation="talys output documentation")
        
        ### add reactant masses
        self.pg.newProduct( GNDincidentnames[yi], mass=getMassFromTalys(GNDincidentnames[yi]) )   
        self.pg.newProduct( nucleusNameFromZA(ZA), mass=getMassFromTalys(nucleusNameFromZA(ZA)) )

    def __getattr__( self, attr ) :
        try : return getattr( self.endlZA, attr )
        except AttributeError : raise AttributeError ("'%s' object has no attribute '%s'" % (self.__class__,attr))

    def mkdir( self ) :
        """Makes directories for running TALYS in. Does not return errors if already there."""
        if not os.path.exists( self.workDir ) :
            os.makedirs( self.workDir )
        talysDir = os.path.join( self.workDir, self.talysDir )
        if not os.path.exists( talysDir ) :
            print 'Creating new TALYS work directory: %s' % talysDir
            os.makedirs( talysDir )
        else:
            print 'Using existing TALYS work directory: %s' % talysDir

    def setup( self, inputs=None, endf=True, restart=False, **kw ) :
        """
Creates a default TALYS input, creates the directory to do the run,
and writes the input files to the directory.
Does not return errors if directory of input files already exist.
Can pass custom input parameters to be added to the
default header and output flags via the input keyword.
Any arguments to the energyList input can be passed as keywords.
"""
        TALYS_INPUTS = talysInput.InputHeader( self )
        if inputs : TALYS_INPUTS += inputs
        outkw = {}
        for i in [ 'spectra', 'gammas', 'discrete', 'channels', 'elastic' ] :
            val = kw.pop(i,None)
            if val : outkw[i] = val
        if endf : outkw['endf'] = endf
        TALYS_INPUTS += talysInput.OutputFlags( **outkw )
        self.mkdir()
        if not restart:
            try :
                talysDir = os.path.join( self.workDir, self.talysDir )
                TALYS_INPUTS.toFile( talysDir )
                talysInput.energyList( **kw ).toFile( talysDir )
                print 'Writing TALYS input files to %s' % talysDir
            except IOError :
                print 'Existing input files found in %s, new input file NOT written' % talysDir
        else:
            import talysFiles
            talysDir = os.path.join( self.workDir, self.talysDir )
            if not os.path.exists( os.path.join( talysDir, 'energies' ) ):
                print 'energies file not found'
                return
            eFile = open( os.path.join( talysDir, 'energies' ), 'r' )
            maxEnergy = talysFiles.TOTAL_TOT( self )['Total'].xMax()
            energies = [float(e.rstrip('\n').strip()) for e in eFile.readlines()]
            energies = [e for e in energies if e > maxEnergy ]
            os.remove( os.path.join( talysDir, 'energies' ) )
            talysInput.energyList( energyList=energies, endf=False ).toFile( talysDir )

    def run( self ) :
        """
Runs TALYS unless talysDir already contains total.tot file.
If directories or input are not there, creates directories with default inputs,
then runs TALYS.
"""
        talysDir = os.path.join( self.workDir, self.talysDir )
        contents = None
        try :
            contents = os.listdir( talysDir )
            if 'total.tot' in contents :
                print 'Using existing TALYS calculation in %s' % talysDir
                return
        except OSError :
            self.setup()
        self.exe()

    def rerun( self ) :
        """
Reruns TALYS after premature ending of the run. Looks to see which energy 
the calculation ended at, and picks up where it left off.
"""
        talysDir = os.path.join( self.workDir, self.talysDir )
        try :
            contents = os.listdir( talysDir )
            if 'total.tot' in contents :
                self.setup( restart=True )
        except OSError :
            self.setup()
        from glob import glob
        fileList = glob(talysDir+'/total.tot')
        fileList += glob(talysDir+'/fission.tot')
        fileList += glob(talysDir+'/xs*.tot')
        fileList += glob(talysDir+'/xs*.L??')
        fileList += glob(talysDir+'/xs*.gam')
        fileList += glob(talysDir+'/??.tot')
        fileList += glob(talysDir+'/??.con')
        fileList += glob(talysDir+'/??.L??')
        if not fileList: raise RuntimeError, 'Bad file list'
        for f in fileList : os.rename( f, f+'.old' )
        self.exe()
        for f in fileList : os.rename( f, f+'.new' )
        for f in fileList :
            ff = open( f+'.old', 'r' )
            fold = ff.readlines()
            ff.close()
            ff = open( f+'.new', 'r' )
            fnew = ff.readlines()
            ff.close()
            for line in fnew:
                if line.startswith('#'):
                    if line.startswith('# # energies'):
                        numen_extra = int(line.split('# # energies =')[1])
                else:
                    fold.append(line)
            for i,line in enumerate(fold):
                if line.startswith('# # energies'):
                    numen_old = int(line.split('# # energies =')[1])
                    fold[i] = line.replace( str(numen_old), str(numen_old+numen_extra) )
                    break
            fcombined = open( f, 'w' )
            fcombined.write( ''.join(fold) )
            fcombined.close()
        for f in fileList : os.remove( f, f+'.old' )
        for f in fileList : os.rename( f, f+'.new' )

    def exe( self ) :
        """Executes the TALYS code"""
        talysDir = os.path.join( self.workDir, self.talysDir )
        try :
            stdin = os.path.join( talysDir, 'input' )
            stdout = os.path.join( talysDir, 'output' )
            print 'Running TALYS calculation in %s' % talysDir
            geftDir = os.path.split( __file__ )[0]
            geftExe = os.path.join( geftDir, 'talys', 'bin', 'talys' )
            if os.path.isfile( geftExe ) :
                print ' using %s' % geftExe
                cmd = geftExe
            else: cmd = 'talys'
            if self.runGrp: cmd = ['sg',self.runGrp,cmd]
            returnCode = subprocess.call( cmd, stdin=open(stdin,'r'), stdout=open(stdout,'w'), \
                                                   stderr=subprocess.STDOUT, cwd=talysDir )
            if returnCode < 0 :
                print 'TALYS terminated by signal',-returnCode
            else :
                self.calculationSuccessful = True
        except OSError, e :
            print 'TALYS execution failed for %s:%s' % ( self.sZA, e )

    def cleanup( self ) :
        """Removes talysDir"""
        os.chdir( self.workDir )
        shutil.rmtree( self.talysDir )
        os.chdir( self.cwd )
        print 'Removed TALYS work directory %s' % self.talysDir

    def shrink( self ) :
        """Creates tarball of talysDir contents."""
        self.tarball = os.path.join( self.workDir, 'talys', fudge.endlmisc.incidentParticleTags( self.yi )[1], self.sZA+'.tar' )
        tarFile = tarfile.open( self.tarball, 'w:gz' )
        tarFile.add( self.talysDir )
        tarFile.close()
        print 'Created tarball of TALYS calculation: %s' % self.tarball

    def expand( self ) :
        """Extracts talys output files from previously shrunk data."""
        tarFile = tarfile.open( self.tarball, 'r:gz' )
        tarFile.extractall()
        tarFile.close()
        print 'Extracted tarball of TALYS calculation: %s' % self.tarball

    def process( self ) :
        """Runs fixQandThresholds on ENDL data."""
        print 'Running fixQandThresholds:'
        self.endlZA.fixQandThresholds( fixThresholdMode = fudge.legacy.endl.endlIClasses.fixThresholdMode_All )
        # Sometimes puts in new threshold when one at same value already exists
        # since thresholds appear different beyond the accuracy given in the ENDL file,
        #so run removeClosePoints to remove it
        for d in self.endlZA.findDatas(I=0):
            d.removeClosePoints()

    def check( self ) :
        """Runs fudge checker on the ENDL data."""
        print 'Running fudge checker:'
        print '\n'.join( map( str, self.endlZA.check() ) )
    
    def getGndReactionSuite(self):
        a = talysGndData.CompleteGndEvaluation(talysZA=self, summed=False, isomers=False)
        
        return self.rsb.finalize()
    
    def write( self, filename ) :
        """ writes GND Reaction Suite to file """
        reactionSuite = self.rsb.finalize()
        reactionSuite.saveToFile( filename+".xml", flags={'verbosity':51} )
