import fudge,talysFiles,talys3d,copy,os
from talysMisc import *
from talysGndMisc import *
from fudge.structure import xensl
import geftMathClasses 
import endlLikeDeltaFunc
from fudge.legacy.endl.endl2dmathClasses import *
from fudge.legacy.endl.endl2 import *
#from fudge.gnd import *
from fudge.particles.nuclear import *
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as pqu
from fudge.evaluationBuilder import *

class CrossSectionType :
    """Base class for types of data that can be put into ENDL from TALYS. Can act as a proxy for the underlying TALYS output file."""
    def __init__( self, talysZA, talysFile, label, **kw ) :
        self.type = type(self)
        self.data = talysFile.get( label )
        #print 'label',label,'talysFile',talysFile 
        self.endlZA = talysZA.endlZA
        self.targetELevel = talysZA.targetELevel
        self.yi = talysZA.endlZA.yi
        self.ZA = talysZA.endlZA.ZA
        self.suffix = talysZA.endlZA.suffix
        self.talysFile = talysFile
        #self.rs = talysZA.gndRS
        for key in kw : setattr( self, key, kw[key] )  ## any other keywords that you pass in the subclass instantiations
        self.MT = -1 ## until caleb fixes the MT lookup/requirement issue
        self.genre = getGenreName(self.C) ## get the genre name from somewhere
        X1 = getattr( self, 'X1', 0. ) # n,n' out energy level
        X4 = getattr( self, 'X4', 0. ) # isomer energy level (S=0)
        
        
        ### NEED: a modern (gnd) function to get the Q value
        self.Q = fudge.legacy.endl.endl2.reactionQ( self.yi, self.ZA, self.C, specialCases = 1, printQWarning = True ) 
        try : self.Q -= self.elevel # isomer energy level
        except AttributeError : pass
        if self.Q is None and self.C == 'fission': 
            print "WARNING: have no fission Q value, setting to canonical 180 MeV"
            self.Q = 180.0
        if self.Q is None:
            self.Q = 0.0    
        
        ### check for excited state
        if X1 : 
            level = X1
        if hasattr(self,'level') : 
            level = 'm%s'%self.level 
            
        #print 'instance data',self.data
            
        
    def __getattr__( self, attr ) :
        try : return getattr( self.talysFile, attr )
        except AttributeError : raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__,attr))

    def __repr__( self ) :
        try :
            reaction = '(' + str( talysYTags( self.endlZA.yi ) ) + ',' + reactionnames[self.C] + ')'
            data = [ reaction, propertynames[self.I] ]
            try : data.append( 'L=%s E=%s Q=%s' % ( self.level, self.elevel, self.Q ) )
            except AttributeError : pass
            if self.yo :
                data.append( incidentnames[self.yo] )
            if self.S :
                data.append( modifiernames[self.S] )
                data.append( self.X1 )
        except : data = [self.type]
        string = ' '.join( map( str, data ) )
        return string

    def toEndl( self ) :
        """
Puts data in self into ENDL. Returns a list of endlIClasses.
Data can then be manipulated from within fudge.endlZA object or by talysZA by proxy.
"""
        
        if self.data == None : return
        if isinstance( self, NonEndlSummedCrossSection ) :
            print "Ignoring %s : No place in ENDL" % self.type
            return
        yi = self.yi
        ZA = self.ZA
        eLevel = self.targetELevel
        yo = self.yo
        C = self.C
        I = self.I
        S = self.S  # qualifier for extra data (S=1 for inelastic)
        X1 = getattr( self, 'X1', 0. ) # n,n' out energy level
        X4 = getattr( self, 'X4', 0. ) # isomer energy level (S=0)
        Q = fudge.legacy.endl.endl2.reactionQ( yi, ZA, C, specialCases = 1, printQWarning = True )
        
        
        try : Q -= self.elevel # isomer energy level
        except AttributeError : pass
        if Q is None and C == 15: 
            print "WARNING: have no fission Q value, setting to canonical 180 MeV"
            Q = 180.0
        try : endlFile = self.endlZA.addFile( yo, C, I, S )
        except: endlFile = self.endlZA.findFile( yo, C, I, S )
        if self.I == 4 :
            l = endlFile.addData( [], eLevel=eLevel, X1 = X1, X4 = X4, Q = Q, temperature = 0 )
            l.setlData( 0, self.data )
        else :
            l = endlFile.addData( self.data, eLevel=eLevel, X1 = X1, X4 = X4, Q = Q, temperature = 0 )
        level = None
        if X1 : level = X1
        if hasattr(self,'level') : level = 'm%s'%self.level
        l.label = fudge.legacy.endl.endl2.reactionEquations( yi, ZA, C, level = level, printQWarning = False )[1]
        return l


        


    def syncLevelEnergies( self, checkOnly = True, isomerMinimumHalflife = None ) :
        """Synchronizes the level energies of the target, and discrete states of the residual
        with a database containing the accepted level energies for all nuclei. The residual
        level energies are synchronized for both the inelastic S=1 data (X1), and isomers (X4).
        Checks that the isomer level labels m? are consistent with the levels given in the database.
        If checkOnly = True then only check against database and print differences,
        else if checkOnly = False then make the changes to the data.
        If isomerMinimumHalflife is given then sets this variable in the Levels class before performing sync."""
        
        if isomerMinimumHalflife : xensl.isomerMinimumHalflife = isomerMinimumHalflife
        targetLevels = xensl.Levels( self.ZA )
        if self.suffix :
            targetEnergy = targetLevels( self.suffix ).energy
            # check isomer level is actually labeled as a metastable state
            if 'm' in self.suffix :
                for m in targetLevels.metastables.keys( ) :
                    if targetEnergy == targetLevels( m ).energy : break
                else:
                    print 'WARNING: Isomer level energy %s, not a known isomer for target ZA=%s' % ( targetEnergy, self.sZA )
        else :
            targetEnergy = 0.

        # check target energy levels
        if targetLevels( self.targetELevel ).energy != targetEnergy :
            print 'WARNING: Isomer label %s with level energy %s for ZA=%s, does not match level energy %s given in database' % ( self.suffix, self.ELevel, self.ZA, targetElevel )
        if self.targetELevel != targetEnergy :
            print 'Changing level energy from %s to %s for target ZA=%s\n%s' % ( self.ELevel, targetEnergy, self.sZA, self )
            if not checkOnly : self.targetELevel = targetEnergy
        residualZA = endl2.residualZA_yos_Q( self.yi, self.ZA, self.C )[0]
        if not residualZA : return
        residualLevels = xensl.Levels( residualZA )
        # check S=1 levels
        if hasattr( self, 'X1' ) :
            if self.S == 1 and self.X1 != residualLevels( self.X1 ).energy :
                print 'Changing discrete state energy from %s to %s for residual ZA=%s\n%s' % ( self.X1, residualLevels( self.X1 ).energy, residualZA, self )
                if not checkOnly : self.X1 = residualLevels( self.X1 ).energy
        # check residual isomer levels
        if hasattr( self, 'X4' ) :
            if self.X4 != residualLevels( self.X4 ).energy :
                print 'Changing isomer level energy from %s to %s for residual ZA=%s\n%s' % ( self.X4, residualLevels( self.X4 ).energy, residualZA, self )
                if not checkOnly : self.X4 = residualLevels( self.X4 ).energy
            # check isomer level is actually labeled as a metastable state
            if self.X4 :
                for m in residualLevels.metastables.keys( ) :
                    if self.X4 == residualLevels( m ).energy : break
                else:
                    print 'WARNING: Isomer level energy %s, not a known isomer for residual ZA=%s\n%s' % ( self.X4, residualZA, self )

    def plot( self ) :
        self.data.plot()

class NonEndlSummedCrossSection( CrossSectionType ) :
    def __init__( self, talysZA, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, type=label, C=None, yo=None, I=None, S=None )   

class ExcitationFunction( CrossSectionType ) :
    def __init__( self, talysZA, C, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=0, I=0, S=0 )
        
class ElasticAngularDistribution( CrossSectionType ) :
    def __init__( self, talysZA, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=10, yo=talysZA.yi, I=1, S=0 )
 
class ElasticLegendreDistribution( CrossSectionType ) :
    def __init__( self, talysZA, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=10, yo=talysZA.yi, I=4, S=0 )
   
class KalbachMannDistribution( CrossSectionType ) : 
    def __init__( self, talysZA, talysFile, label ) :
        pass
        #CrossSectionType.__init__( self, talysZA, talysFile, label, C=10, yo=talysZA.yi, I=1, S=0 )
    
class EnergySpectra( CrossSectionType ) :
    def __init__( self, talysZA, C, yo, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=yo, I=4, S=0 )
    
class GammaMultiplicity( CrossSectionType ) :
    def __init__( self, talysZA, C, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=7, I=9, S=0 )

class DiscreteLevelExcitation( CrossSectionType ) :
    def __init__( self, talysZA, C, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=0, I=0, S=1, X1=talysFile.X1 )
    
class DiscreteLevelAngularDistribution( CrossSectionType ) :
    def __init__( self, talysZA, C, yo, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=yo, I=1, S=1, X1=talysFile.X1 )
    
class DiscreteLevelGammaSpectra( CrossSectionType ) :
    def __init__( self, talysZA, C, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=7, I=4, S=1, X1=talysFile.X1 )
    
class DiscreteLevelGammaMultiplicity( CrossSectionType ) :
    def __init__( self, talysZA, C, talysFile, label ) :
        CrossSectionType.__init__( self, talysZA, talysFile, label, C=C, yo=7, I=9, S=1, X1=talysFile.X1 )
    

def findFunction( x, **kw ) :
    for key in kw :
        try :
            if getattr( x, key ) not in kw[key] : return False
        except TypeError :
            try :
                if getattr( x, key ) != kw[key] : return False
            except : raise
        except AttributeError :
            return False
    return True

class endlList( list ) :
    """Base class for lists of data to be put into ENDL."""
    def __init__( self, dataList=[] ) :
        list.__init__( self, dataList )

    def __add__( self, other ) :
        return endlList( list.__add__( self, other ) )

    def __repr__( self ) :
        return '\n'.join( map( str, enumerate( self ) ) )

    def toEndl( self ) :
        """Converts endlList of CrossSectionType's to fudge.endlIClasses, returning a list of correspsonding endlIClasses"""
        return [ data.toEndl() for data in self ]
    
    def toGnd( self ) :
        """Converts endlList of CrossSectionType's to fudge.endlIClasses, returning a list of correspsonding endlIClasses"""
        return [ data.toGnd() for data in self ]
    
    def append( self, other ) :
        """Modifies built in append function to only add other to list if other actually contains data"""
        if other.data : self += [other]

    def filter( self, function ) :
        """
Variation on the filter function, calling this function as a method of the endlList, returns and endlList object
(if you use the filter function on the endlList as filter(function,endlList) you just get a normal list returned).
"""
        return endlList( filter( function, self ) )

    def findDatas( self, **kw ) :
        """Returns endlList of data matching the arguments given. Arguements can have unique values or be a list of tuple of values."""
        return self.filter( lambda x: findFunction( x, **kw ) )

    def findData( self, **kw ) :
        """Returns data of endlList item that matches the arguments given. Raises exception if None or more than 1 data found."""
        datalist = self.findDatas( **kw )
        if not datalist : raise RuntimeError, "Error in endlList.findData : No data found for %s" % kw
        if len(datalist)==1 : return datalist[0].data
        else: raise RuntimeError, "Error in endlList.findData : More than one data found for %s\n%s" % (kw,datalist)

    def getIndex( self, **kw ) :
        """Gets the index in endlList of the data matching the arguments given. Raises exception if None or more than 1 data found."""
        self.findData( **kw ) # call this first to catch errors in None or more than one data found
        for i,d in enumerate( self ) :
            if findFunction( d, **kw ) : return i

    def fixXArrays( self ):
        """Fixes xArrays for a group of data so that they all match and fudge checker doesn't complain."""
        I0List = endlList( [ d for d in self if isinstance( d, ExcitationFunction ) ] )
        I4List = endlList( [ d for d in self if isinstance( d, EnergySpectra ) ] )
        energyList = []
        for I0 in I0List :
            energyList = I0.data.xArray().data
        energyList.reverse()
        for I4 in I4List :
            for energy in energyList :
                if energy < I4.data.xMin() :
                    if debug : print 'Inserting (fake) delta function for %s at %s' % ( I4, energy )
                    I4.insertThreshold( energy=energy )

    def removeSmallChannels( self, smallXS=1e-6, excludeCList=[], reaction=None, fraction=0.01 ) :
        """
Removes all data for a channel if the integrated cross section is below smallXS.
Or if reaction != None, then the integrated cross section is less than fraction
of the reaction cross section, at all energies.
"""
        smallCList = []
        for d in filter( lambda x: isinstance( x, ExcitationFunction ), self ) :
            if not reaction:
                if d.data.max() < smallXS :
                    if d.C not in excludeCList : smallCList.append( d.C )
            else:
                fudge.endl2dmathClasses.doSafeDivide = True
                if ( d.data / reaction ).max() < fraction :
                    if d.C not in excludeCList : smallCList.append( d.C )
        if smallCList == []: return self
        print "Removing all small channels:",smallCList
        return self.filter( lambda x: x.C not in smallCList )

    def syncLevelEnergies( self, checkOnly=True, isomerMinimumHalflife=None ) :
        """"""
        for data in self :
            data.syncLevelEnergies( checkOnly=checkOnly, isomerMinimumHalflife=isomerMinimumHalflife )
    

class Summed( endlList ) :
    """endlList subclass for the total cross section."""
    def __init__( self, talysZA, C, *args ) :
        endlList.__init__( self )
        for label in args :
            if label == 'Total' :
                print 'channel:',C
                xs =  ExcitationFunction( talysZA, C, talysFiles.TOTAL_TOT( talysZA ), label ) 
                self.append( xs )
               
                ### build GND version
                #rb = summedReactionBuilder( name= 'total' )
                #rb.addCrossSection( crossSectionBuilder( data = xs.data.data, form = "Pointwise" ) )  ### add cross section from endl2dmath
                ### need finalize method for summed
                
                #x = nucleusNameFromZA(ZAforYTags(talysZA.yi))
                #part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x) , multiplicity=1)
                #if not CrossSectionOnly: 
                #   part.addDistribution( 'AngularPointwise', data = angDist.data.data ) 
                #   #part.addDistribution( 'AngularLegendre', data = legDist.data.data ) ### Do we need Legendre dists?
                #rb.addProduct( part )

                #rb.addResidual( )  
                #talysZA.rsb.addSummedReaction( rb )


            elif label == 'Elastic' :
               print "Ignoring Elastic: elastic cross section accessed via ElasticChannel class"
            elif label in talysFiles.TOTAL_TOT( talysZA ) :
               self.append( NonEndlSummedCrossSection( talysZA, talysFiles.TOTAL_TOT( talysZA ), label ) )

class FissionChannel( endlList ) :
    """endlList subclass for the fission channel."""
    def __init__( self, talysZA, C ) :
        endlList.__init__( self )
        #C = 15
        #C = 'fission'
        print 'fission channel:',C
        xs =  ExcitationFunction( talysZA, C, talysFiles.FISSION_TOT( talysZA ), 'CrossSection' ) 
        self.append( xs )
        
        ### now add to GND
        #rb = fissionReactionBuilder()
        #rb.addCrossSection( crossSectionBuilder( data = xs.data.data, form = "Pointwise" ) )  ### add cross section from endl2dmath
        ### need fission products
        
        #x = nucleusNameFromZA(ZAforYTags(talysZA.yi))
        #x = 'n'
        #m = multiplicityFromFission
        #part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x), multiplicity=m )
        #part.addDistribution( 'EnergySpectrum', data = angDist.data.data ) 
        #rb.addProduct( part )
        
        
        #rb.addResidual( )  
        #talysZA.rsb.addReaction( rb )


class ElasticChannel( endlList ) :
    """endlList subclass for the elastic channel."""
    def __init__( self, talysZA, C, CrossSectionOnly=False ) :
        endlList.__init__( self )
        #C = 10
        #C = 'elastic'
        print 'elastic channel:',C
        #npdtha = npdthagFromName('n')
        pe = talysYTags( talysZA.yi ) * 2
        if talysZA.suffix :
            if talysZA.suffix == 'm' : level = 1
            else : level = int( talysZA.suffix[1:] )
        else : level = 0
        xs = ExcitationFunction( talysZA, C, talysFiles.TOTAL_TOT( talysZA ), 'Elastic' )
        self.append( xs )
           
        ### GND version
        rb = reactionBuilder( )
        rb.addCrossSection( crossSectionBuilder( data = xs.data.data, form = "Pointwise" ) )  ### add cross section from endl2dmath
        
        x = nucleusNameFromZA(ZAforYTags(talysZA.yi))
        part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x) , multiplicity=1)
        if not CrossSectionOnly: 
            angDist = ElasticAngularDistribution( talysZA, talys3d.BinaryAngles( talysZA, pe, level ), 'AngularSpectrum' ) 
            self.append( angDist )
            part.addDistribution( 'AngularPointwise', data = angDist.data.data ) # endl3dmath object passes OK
            ### Do we need Legendre dists?
            #legDist = ElasticLegendreDistribution( talysZA, talys3d.BinaryLegendreCoeffs( talysZA, pe, level ), 'LegendreSpectrum' ) 
            #self.append( legDist )
            #print "legDist ",type(legDist), legDist.data
            #part.addDistribution( 'AngularLegendre', data = legDist.data.data ) # endl3dmath object passes OK
        rb.addProduct( part )
        
        rb.addResidual( )  # by default the residual particle distribution is 'recoil'
        talysZA.rsb.addReaction( rb )


class PlainChannel( endlList ) :
    """endlList subclass for a plain channel."""
    def __init__( self, talysZA, C, CrossSectionOnly=False, isomers=False ) :
        endlList.__init__( self )
        rb = reactionBuilder( )
        if sum([int(x) for x in npdthaFromName(C)]) <= 1 : return
        #print 'plain channel:',C
        npdtha = particleList( C )
        
        npdtha_val = [int(x) for x in npdthaFromName(C)] ## npdtha number, as a string ( 2nd = ['2','0','1','0','0','0'] )
        ZA_adjust = sum( [ ZAforYTags(i+1)*npdtha_val[i]  for i in range(6) ] )
        prod = nucleusNameFromZA(talysZA.ZA + ZAforYTags(talysZA.yi) - ZA_adjust)
        talysZA.pg.newProduct( prod, mass=getMassFromTalys(prod)) 


        #print npdtha
        isomerList = talys3d.IsomerCrossSectionList( talysZA, C )
        if isomers and isomerList :
            for isomer in isomerList :
                print 'C:',C
                xsIso = ExcitationFunction( talysZA, C, isomer, 'CrossSection' ) 
                self.append(xsIso)
                rb.addCrossSection( crossSectionBuilder( data = xsIso.data.data, form = "Pointwise" ) )   
                talysZA.pg.setExcitationEnergy( prod+'_e'+str(isomer.level), isomer.X4, 'MeV' )

                
            #if not isomerList:
            #    xs = ExcitationFunction( talysZA, C, talysFiles.XSnpdtha_TOT( talysZA, npdtha ), 'CrossSection' )     
            #    self.append(xs)
        else :        
            #print 'entering here'
            xs = ExcitationFunction( talysZA, C, talysFiles.XSnpdtha_TOT( talysZA, npdtha ), 'CrossSection' ) 
            self.append(xs)
        
        #print "What is xs type: ",type(xs.data)    
            
        if not xs.data:
           print 'no channel:',C
           return
        else:   
           print 'plain channel:',C
           rb.addCrossSection( crossSectionBuilder( data = xs.data.data, form = "Pointwise" ) )   # endl2dmath object passes OK
           
           if not CrossSectionOnly: 
               gammaMult =  GammaMultiplicity( talysZA, C, talysFiles.XSnpdtha_TOT( talysZA, npdtha ), 'GammaMultiplicity' ) 
               self.append( gammaMult )
               
               npdtha_val = [int(x) for x in npdthaFromName(C)] ## npdtha number, as a string ( 2nd = ['2','0','1','0','0','0'] )
               for yo in yoList( C ) :
                   ### add particle distribution
                   yoSpec = EnergySpectra( talysZA, C, yo, talys3d.ContinuumEnergySpectra( talysZA, npdtha ), particleName( yo )+'EnergySpectrum' ) 
                   self.append( yoSpec )
                   
                   x = nucleusNameFromZA(ZAforYTags(yo))
                   m = npdtha_val[yo-1] 
                   part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x) , multiplicity=m  )
                   if yoSpec: 
                      ### Need to change to energy-energy distribution and add isotropic angular as uncorrelated
                      part.addDistribution( 'Uncorrelated', data = {'Energy': ('EnergyPointwise',yoSpec.data.data), 'Angular':('Isotropic',None) }) 
                      
                      #Edist = distributionBuilder( 'EnergyPointwise', data = yoSpec.data.data )
                      #Adist = distributionBuilder( 'Angular', data = None )
                      #part.addDistribution( 'Uncorrelated', data = {'Energy': Edist, 'Angular': Adist})  ### Nice to read
                   rb.addProduct( part )

               ### now get gammas 
               yo=7
               gammaSpec = EnergySpectra( talysZA, C, yo, talys3d.GammaEnergySpectra( talysZA, npdtha ), particleName( yo )+'EnergySpectrum' ) 
               part = talysZA.pg.newProduct( 'gamma', mass='0 amu' , multiplicity=m )
               part.addDistribution( 'EnergyPointwise', data = gammaSpec.data.data )
               #part.addDistribution( 'Isotropic' ) ### add that they are uncorrelated
               rb.addProduct( part )
       
               self.append( gammaSpec )
               self.fixXArrays()
        
        
        
        rb.addResidual( )  # by default the residual particle distribution is 'recoil'
        talysZA.rsb.addReaction( rb )
        
class BinaryChannel( endlList ) :
    """endlList subclass for a binary channel with discrete level cross sections."""
    def __init__( self, talysZA, C, **kw ) :
        endlList.__init__( self )
        npdtha = particleList( C )
        npdtha_val = [int(x) for x in npdthaFromName(C)]                            ### npdtha number, as a string ( 2nd = ['2','0','1','0','0','0'] )
        ZA_adjust = sum( [ ZAforYTags(i+1)*npdtha_val[i]  for i in range(6) ] )     ### get the effective ZA of the outgoing particle
        prod = nucleusNameFromZA(talysZA.ZA + ZAforYTags(talysZA.yi) - ZA_adjust)   ### get the residual nucleus name
        talysZA.pg.newProduct( prod, mass=getMassFromTalys(prod))                   ### 
        z,a = fudge.particles.nuclear.getZandAFromName(prod)                        ### for looking up items in talys structure files
        x = nucleusNameFromZA(ZA_adjust)                                            ### name of outgoing particle (besides residual nucleus)
        pe = inelasticPrefix( talysZA.yi, C )                                       ### alphabetic tag for binary channel (np,na,nt,etc.)
        yo = yoList( C )[0]                                                         ### numerical tag of outgoing part
        #print 'yo',yo
        
        gammas = True
        if 'gammas' in kw: gammas = kw['gammas']
        
        print 'binary channel:',C,'products:',prod,x
        
        ### New reaction for total binary
        rb = reactionBuilder( )
         
        #### getting the total XS. Is this the right function for total binary
        #### **** Need something better for the total cross section
        xs = ExcitationFunction( talysZA, C, talysFiles.pe_TOT( talysZA, pe ), 'CrossSection' ) 
        self.append( xs )
        if xs.data:
            ### check if the data is 0.0
            rb.addCrossSection( crossSectionBuilder( data = xs.data.copyDataToXYs(), form = "Pointwise" ) )   
        
        #### get particle spectrum
        partSpec = EnergySpectra( talysZA, C, yo, talys3d.ContinuumEnergySpectra( talysZA, npdtha ), particleName( yo )+'EnergySpectrum' ) 
        self.append( partSpec )
        part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x) , multiplicity=1  )
        if partSpec.data : 
            part.addDistribution( 'EnergyPointwise', data = partSpec.data.data ) 
        rb.addProduct( part )
        
        
        #### Get the gamma spectrum
        if gammas:
            gamSpec = EnergySpectra( talysZA, C, 7, talys3d.GammaEnergySpectra( talysZA, npdtha ), 'GammaEnergySpectrum' )
            self.append( gamSpec )
            part = talysZA.pg.newProduct( 'gamma', mass='0 amu' )
            if gamSpec.data : 
                part.addDistribution( 'EnergyPointwise', data = gamSpec.data.data )
            rb.addProduct( part )
        
        
        rb.addResidual( )  # by default the residual particle distribution is 'recoil'
        talysZA.rsb.addReaction( rb )
        
    
        ### each discrete level transition
        for discreteLevel in talys3d.BinaryDiscreteList( talysZA, pe ) :
            print '    discrete level:', discreteLevel.level

            rb = reactionBuilder( ) ## new reaction for this 

            ### add excited state and excitation energy
            ExSt = prod+"_e"+str(discreteLevel.level)
            level = getLevelFromTalys(prod,discreteLevel.level) 
            talysZA.pg.newProduct( ExSt, mass=getMassFromTalys(prod), levelIndex = level.level, levelEnergy=pqu(level.energy,'MeV') ) 

            levelXS = DiscreteLevelExcitation( talysZA, C, discreteLevel, 'CrossSection' ) 
            self.append( levelXS )
            if levelXS.data : 
                rb.addCrossSection( crossSectionBuilder( data = levelXS.data.copyDataToXYs(), form = "Pointwise" ) )   

            ### get ang Distribution of outgoing particle
            part = talysZA.pg.newProduct( x, mass=getMassFromTalys(x) , multiplicity=1  )
            for discreteAngLevel in talys3d.BinaryDiscreteAngleList( talysZA, pe ) :
                if discreteAngLevel.level == discreteLevel.level :
                    angDist = DiscreteLevelAngularDistribution( talysZA, C, yo, discreteAngLevel, 'AngularSpectrum' )
                    #part.addDistribution( 'AngularPointwise', data = angDist.data.data )
                    try:
                        part.addDistribution( 'AngularPointwise', data = angDist.data.data )
                    except:
                        print '   **** Still have a problem with this angular distribution: level=', discreteAngLevel.level
                        print '   **** is your instance of talys data complete? '
                        pass
            rb.addProduct( part )

            rb.addResidual( excitationLevel = discreteLevel.level)  # by default the residual particle distribution is 'recoil'
            talysZA.rsb.addReaction( rb )


class CompleteGndEvaluation( endlList ) :
    """endlList subclass for a complete evaluation."""
    def __init__( self, talysZA, summed=False, isomers=False, **kw ) :
        endlList.__init__( self )
        talysDir = os.path.join( talysZA.workDir, talysZA.talysDir )
        contents = None
        try :
            contents = os.listdir( talysDir )
            if 'total.tot' not in contents :
                err = "No evaluation found in %s" % talysDir
                raise RuntimeError, err
        except OSError :
            err = "talys directory (%s) does not exist" % talysDir
            raise RuntimeError, err
        for reaction in talysZA.reactionList : #talysZA.CList :
            if reaction == 'total' :
                self += Summed( talysZA, reaction, 'Total' )
            elif reaction == 'elastic' :
               self += ElasticChannel( talysZA, reaction )
            elif reaction == 'fission' :
                self += FissionChannel( talysZA, reaction )
            #elif reaction == 'capture' :               # capture channel?
            #    pass #self += PlainChannel( talysZA, reaction, isomers=isomers )
            elif sum([int(x) for x in npdthaFromName(reaction)]) == 1:  # binary channel         #reaction in talysZA.binaryCList :
                self += BinaryChannel( talysZA, reaction, **kw )
            elif sum([int(x) for x in npdthaFromName(reaction)]) > 1 : # other regular channel    #11<=C<=14 or ( 16<=C<=49 ) :
                self += PlainChannel( talysZA, reaction, isomers=isomers )
            else :
                print "Warning: Channel %s type not accounted for" % reaction 
        if summed :
            self += Summed( talysZA, 'NonElastic' )
            self += Summed( talysZA, 'ShapeElastic' )
            self += Summed( talysZA, 'CompElastic' )
            self += Summed( talysZA, 'Reaction' )

class CompleteEvaluation( endlList ) :
    """endlList subclass for a complete evaluation."""
    def __init__( self, talysZA, summed=False, isomers=False, **kw ) :
        endlList.__init__( self )
        talysDir = os.path.join( talysZA.workDir, talysZA.talysDir )
        contents = None
        try :
            contents = os.listdir( talysDir )
            if 'total.tot' not in contents :
                err = "No evaluation found in %s" % talysDir
                raise RuntimeError, err
        except OSError :
            err = "talys directory (%s) does not exist" % talysDir
            raise RuntimeError, err
        for C in talysZA.CList :
            if C == 1 :
                self += Summed( talysZA, 'Total' )
            elif C == 10 :
                self += ElasticChannel( talysZA )
            elif C == 15 :
                self += FissionChannel( talysZA )
            elif C in talysZA.binaryCList :
                self += BinaryChannel( talysZA, C, **kw )
            elif 11<=C<=14 or ( 16<=C<=49 ) :
                self += PlainChannel( talysZA, C, isomers=isomers )
            else :
                raise Exception( "Channel %s type not accounted for" % C )
        if summed :
            self += Summed( talysZA, 'NonElastic' )
            self += Summed( talysZA, 'ShapeElastic' )
            self += Summed( talysZA, 'CompElastic' )
            self += Summed( talysZA, 'Reaction' )

class CrossSections( endlList ):
    """endlList subclass for a Cross Sections evaluation."""
    def __init__( self, talysZA, CList=None, isomers=False ) :
        endlList.__init__( self )
        talysDir = os.path.join( talysZA.workDir, talysZA.talysDir )
        contents = None
        try :
            contents = os.listdir( talysDir )
            if 'total.tot' not in contents :
                err = "No evaluation found in %s" % talysDir
                raise RuntimeError, err
        except OSError :
            err = "talys directory (%s) does not exist" % talysDir
            raise RuntimeError, err
        if not CList: CList = talysZA.CList
        for C in CList :
            if C == 1 :
                self += Summed( talysZA, 'Total' )
            elif C == 10 :
                self += ElasticChannel( talysZA, CrossSectionOnly=True )
            elif C == 15 :
                self += FissionChannel( talysZA )
            elif 11<=C<=14 or ( 16<=C<=49 ) :
                self += PlainChannel( talysZA, C, CrossSectionOnly=True, isomers=isomers )
            else :
                raise Exception( "Channel %s type not accounted for" % C )
