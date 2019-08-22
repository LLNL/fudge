import os
import math
import copy
import weakref
import UserDict

import talysStructure
import talysMisc
#from talysMisc import *
from fudge.legacy.endl.endl2dmathClasses import endl2dmath

class datatype :
    """Bass class for datatypes read in from TALYS output files and put into self.data."""
    def getValue( self ) :
        """Returns data read in from TALYS."""
        return self.data

class Energy( datatype ) :
    """Energy data is stored as is."""
    def __init__( self, datapoint ) :
        self.data = float( datapoint )

class CrossSection( datatype ) :
    """
Cross Section data is converted from mb in TALYS to b for ENDL.
Special case 1.00000E-07 is changed to zero.
"""
    def __init__( self, datapoint ) :
        self.data = float( datapoint ) / 1000 # mb -> b
        # Sometimes TALYS 1e-7 when it could well be 0, bug ? replace with zero is a workaround
        if datapoint == "1.00000E-07" : self.data = float( 0 )

class Angle( datatype ) :
    """Angle data is converted from degrees in TALYS to mu=cos(angle)."""
    def __init__( self, datapoint ) :
        angleDegrees = float( datapoint )
        mu = math.cos( math.radians( angleDegrees ) )
        self.data = mu

class Coeff( datatype ) :
    """Legendre Coefficients returned as is for now. (May have to add some normalization later.) """
    def __init__( self, datapoint ) :
        #angleDegrees = float( datapoint )
        #mu = math.cos( math.radians( angleDegrees ) )
        self.data = datapoint

class Order( datatype ) :
    """Legendre Coefficient order returned as is."""
    def __init__( self, datapoint ) :
        self.data = datapoint

class Ratio( datatype ) :
    """Ratio data is stored as is."""
    def __init__( self, datapoint ) :
        self.data = float( datapoint )

class Counter( datatype ) :
    """Counter data is stored as is."""
    def __init__( self, datapoint ) :
        self.data = int( datapoint )

class OutputFile( UserDict.IterableUserDict ) :
    """
Class containing general methods for reading TALYS output files,
first stripping off the header lines starting with #,
while storing any auxilary info contained in the header,
possible info stored in header are:
Q, spin, parity, elevel (Energy of level), threshold, numen(number of energies).
The data is stored in self.datalines as a list, with each item in the
list being a list corresponding to the list of data points on 1 line.
All data read in from the talys output files are endl2dmath objects.
"""
    normalize = False
    
    def __init__( self, talysZA ) :
        UserDict.IterableUserDict.__init__( self )
        self.datalines = []
        #print ' **** New File: ',self.filename
        try : 
            OutputFile = open( os.path.join( talysZA.workDir, talysZA.talysDir, self.filename ), 'r')
        except IOError: 
            #print 'WARNING: failed to open %s' % self.filename
            return
        except (AttributeError,IOError): 
            return
        #except IOError: 
        #    print 'WARNING: failed to open %s' % self.filename
        #    return
        for dataline in OutputFile :
            if dataline[0] == '#' :
                if dataline[2:9] == 'Q-value' :
                    self.Q = float( dataline[14:26] ) # Note Q value not used in the end, but useful for calculating X1 for inelastic levels
                    try :
                        self.spin = float( dataline[32:37] )
                        self.parity = str( dataline[46:47] )
                    except :
                        try : self.elevel = float( dataline[34:45] )
                        except : pass
                if dataline[2:13] == 'E-threshold' :
                    self.threshold = float( dataline[14:26] )
                if dataline[2:12] == '# energies' :
                    self.numen = int( dataline.split( '=' )[1] )
                continue
            elif dataline[0] == '' :
                break
            #print dataline    
            self.datalines.append( dataline.split() )
        OutputFile.close()

    def __repr__( self ) :
        output = ''
        for label, data in self.items() :
            output += '\n'.join( map( str, [ label, data ] ) )
        return output
        
    def set2dData( self, xdatatype, columnLabels ) :
        """
Reads in the data from datalines list in self and puts data into endl2dmath object.
The data is then processed to tidy up data for ENDL formats.
"""
        for column, label in enumerate( columnLabels ) :
            try:
                data = []
                for datapoints in self.datalines :
                    if column+1 >= len(datapoints) : continue
                    if xdatatype == 'energy' :
                        if 'Branching' in label or 'Ratio' in label :
                            data.append( [ Energy( datapoints[0] ).getValue(), Ratio( datapoints[column+1] ).getValue() ] )
                        else :
                            data.append( [ Energy( datapoints[0] ).getValue(), CrossSection( datapoints[column+1] ).getValue() ] )
                    elif xdatatype == 'angle' :
                        data.insert( 0, [ Angle( datapoints[0] ).getValue(), CrossSection( datapoints[column+1] ).getValue() ] )
                    elif xdatatype == 'coeff' :
                        data.append( [Order( datapoints[0] ).getValue(), Coeff( datapoints[column+1] ).getValue() ] )    
                if data :
                    self[label] = endl2dmath( data, label=label )
            except (ValueError,IndexError): print 'set2dData failed for column %s of file %s'%(label,self.filename)
        
        self.process()

    def process( self ) :
        """List of function calls to process data ready for ENDL format."""
        self.removeZeros()
        self.removeEmptyData()
        self.removeBelowEmin()
        self.fixSingleEnergyPoint()
        self.insertThreshold()
        self.calculateMultiplicity()
        self.normalizeData()

    def fixThresholdBuginTalys( self, talysZA ) :
        """
Fix bug in talys: if XSnpdtha.TOT file corresponds to a inelastic channel,
then the threshold is given as the threshold of the first inelastic level,
rather than the continuum threshold. The threshold given in the pe.CON file is used instead.
"""
        if isinstance( self, XSnpdtha_TOT ) :
            if talysMisc.eFromParticleList( self.npdtha ) :
                CONFile = pe_CON( talysZA, talysMisc.talysYTags( talysZA.yi ) + talysMisc.eFromParticleList( self.npdtha ) )
            else : return
            try:                                                                                                                                       
                if talysMisc.debug : print "Changing threshold for %s(%s) from %s to %s" % ( self.__class__, self.npdtha, self.threshold, CONFile.threshold )    
                self.threshold = CONFile.threshold                                                                                                     
            except AttributeError: pass                                                                                                                

    def fixSingleEnergyPoint( self ) :
        """
For some channels e.g. na, there is sometimes only a single energy point in the outgoing spectrum,
change this to a (fake) delta function.
"""
        for label, data in self.items() :
            if len( data ) == 1 :
                if talysMisc.debug : print 'Found single Eout for %s, replacing with (fake) delta function' % label
                x, y = data[0]
                deltafn = endlLikeDeltaFunction.deltafn( x, h=y )
                data.set( deltafn )

    def removeBelowEmin( self ) :
        """
The minimum energy that is possible for outgoing spectra and angular distributions is fixed due to the
filename formatting containing this energy. So slice data below this energy.
The minimum energy is set as a variable in talysMisc as smallestEnergyForSpectra.
"""
        for label,data in self.items() :
            if hasattr( self, 'energy' ) :
                if self.energy < talysMisc.smallestEnergyForSpectra :
                    if talysMisc.debug : print 'Removing Spectra/Angle data %s at %s because below min energy %s' % ( label, self.energy, smallestEnergyForSpectra )
                    del self[label]

    def removeZeros( self ) :
        """Trim extra zeros from data."""
        for data in self.values() :
            data.trim()
    
    def insertThreshold( self ) :
        """Insert threshold value if found in self."""
        if hasattr( self, 'threshold' ) : 
            if self.threshold :
                for label, data in self.items() :
                    thresholdValue = 0
                    if isinstance( self, pe_LXX ) : thresholdValue = data.getValue( data.xMin() )
                    if talysMisc.debug : print 'Setting threshold for %s.%s at %s' % ( self.filename, label, self.threshold )
                    data.setValue( self.threshold, thresholdValue )
                    if data.xMin() < self.threshold :
                        data.trim()

    def removeEmptyData( self ) :
        """Delete data containing all zeros."""
        for label,data in self.items() :
            if not data :
                if talysMisc.debug : print 'Removing empty data %s' % label
                del self[label]

    def normalizeData( self ) :
        """If self.normalize is set in subclass, renormalize data into a probability distribution."""
        for label,data in self.items() :
            if self.normalize :
                newLabel = label.split('CrossSection')[0]+'Probability'
                if talysMisc.debug : print '%s = Normalized data' % newLabel
                self[newLabel] = data.normalize()
                self[newLabel].label = newLabel

    def calculateMultiplicity( self ) :
        """If find Gamma data in self, then add to self the corresponding multiplicity by dividing by the cross section."""
        for label,data in self.items() :
            if 'Gamma' in label :
                labelSplit = label.split('Gamma')
                denominatorLabel = labelSplit[0]+labelSplit[1]
                newLabel = label.split('CrossSection')[0]+'Multiplicity'
                if denominatorLabel in self.data :
                    if talysMisc.debug : print 'Setting %s = %s / %s' % ( newLabel, label, denominatorLabel )
                    self[newLabel] = talysMisc.getMultiplicity( self[label], self[denominatorLabel] )
                    self[newLabel].label = newLabel

class OnePerArgList :
    """
Bass class for all TALYS output files.
All output files are singletons, so all output file subclasses should inherit from this class.
Their dicts are stored in a weakref dict which gets garbage collected once talysZA is dereferenced.
There is one unique instance for each set of talysZA, subclass, argument list.
"""
    __dict = weakref.WeakKeyDictionary()
    def __init__( self, talysZA, *args ) :
        newInstance = False
        if talysZA not in self.__dict :
            self.__dict[talysZA] = {}
        if args not in self.__dict[talysZA] :
            self.__dict[talysZA][args] = {}
            newInstance = True
        self.__dict__ = self.__dict[talysZA][args]
        self.newInstance = newInstance

class TOTAL_TOT( OnePerArgList, OutputFile ) :
    """
Class for 'total.tot'.
Contains various summed cross sections.
KEYS:
    NonElastic
    Elastic
    Total
    CompElastic
    ShapeElastic
    Reaction
    CompNonElastic
    Direct
    PreEquilibrium
"""
    def __init__( self, talysZA ) :
        OnePerArgList.__init__( self, talysZA, self.__class__ )
        if self.newInstance :
            self.columnLabels = ( "NonElastic", \
                                  "Elastic", \
                                  "Total", \
                                  "CompElastic", \
                                  "ShapeElastic", \
                                  "Reaction", \
                                  "CompNonElastic", \
                                  "Direct", \
                                  "PreEquilibrium" )
            self.filename = "total.tot"
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )

class FISSION_TOT( OnePerArgList, OutputFile ) :
    """
Class for 'fission.tot'.
Contains the fission cross section.
KEYS:
    CrossSection
"""
    def __init__( self, talysZA ) :
        OnePerArgList.__init__( self, talysZA, self.__class__ )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", )
            self.filename = "fission.tot"
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )

class XSnpdtha_TOT( OnePerArgList, OutputFile ) :
    """
Class for 'xs??????.tot'.
Contains the total exclusive cross section for each channel
and the total gamma cross section for the exclusive channel.
INPUTS:
    particleList = particle list of the form npdtha where the number
                   describes the number of particles of each type.
KEYS:
    CrossSection
    GammaCrossSection
""" 
    def __init__( self, talysZA, particleList ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "GammaCrossSection" )
            self.npdtha = particleList
            self.filename = 'xs' + particleList + '.tot'
            OutputFile.__init__( self, talysZA )
            self.fixThresholdBuginTalys( talysZA )
            self.set2dData( 'energy', self.columnLabels )



class XSnpdtha_LXX( OnePerArgList, OutputFile ) :
    """
Class for 'xs??????.L??'.
Contains the cross section to each isomer labeled by the level L??,
which sum up to the total cross section in 'xs??????.tot',
and the branching ratio to that isomer.
Should be called from talys3d.IsomerCrossSectionList( talysZA, particleList )
which will then contain a list of these classes for all levels found.
INPUTS:
    particleList = particle list of the form npdtha where the number
                   describes the number of particles of each type.
    level        = level number of isomer.
KEYS:
    CrossSection
    Branching
""" 
    def __init__( self, talysZA, particleList, level ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, particleList, level )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "Branching" )
            self.npdtha = particleList
            self.level = int(level)
            self.filename = 'xs' + particleList + '.L' + str( level ).zfill( 2 )
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )
            self.X4 = self.elevel

class eSPECenergy_TOT( OnePerArgList, OutputFile ) :
    """
Class for '?spec???.???.tot'.
Contains the particle spectra cross sections and ratios as function of outgoing energy
for different processes contributing to the particle spectra.
Should be called from talys3d.PreEquilibriumRatio( talysZA, ejectile )
which will then contain a list of these classes for all incident energies found.
INPUTS:
    ejectile     = e, string containing the ejectile type (e=npdthag)
    energy       = incident energy
KEYS:
    Total
    Direct
    PreEquilibrium
    MultiplePreEquilibrium
    Compound
    PreEquilibriumRatio
""" 
    def __init__( self, talysZA, ejectile, energy ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, ejectile, energy )
        if self.newInstance :
            self.columnLabels = ( "Total", \
                                  "Direct", \
                                  "PreEquilibrium", \
                                  "MultiplePreEquilibrium", \
                                  "Compound", \
                                  "PreEquilibriumRatio" )
            self.energy = Energy( energy ).getValue()
            ENERGY = ("%.3f" % self.energy ).zfill(7)
            self.filename = ejectile + 'spec' + ENERGY + '.tot'
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )
            
class SPnpdthaEenergy_TOT( OnePerArgList, OutputFile ) :
    """
Class for 'sp??????E???.???.tot'.
Contains the paricle cross section as function of outgoing energy for all outgoing particles,
for a single incident energy defined in the filename. 
Should be called from talys3d.ContinuumEnergySpectra( talysZA, particleList )
which will then contain a list of these classes for all incident energies found.
INPUTS:
    particleList = particle list of the form npdtha where the number
                   describes the number of particles of each type.
    energy       = incident energy
KEYS:
    GammaCrossSection, GammaProbability
    NeutronCrossSection, NeutronProbability
    ProtonCrossSection, ProtonProbability
    DeuteronCrossSection, DeuteronProbability
    TritonCrossSection, TritonProbability
    Helium3CrossSection, Helium3Probability
    AlphaCrossSection, AlphaProbability
""" 
    normalize = True

    def __init__( self, talysZA, particleList, energy ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, particleList, energy )
        if self.newInstance :
            self.columnLabels = ( "GammaCrossSection", \
                                  "NeutronCrossSection", \
                                  "ProtonCrossSection", \
                                  "DeuteronCrossSection", \
                                  "TritonCrossSection", \
                                  "Helium3CrossSection", \
                                  "AlphaCrossSection" )
            self.energy = Energy( energy ).getValue()
            ENERGY = ("%.3f" % self.energy ).zfill(7)
            self.filename = 'sp' + particleList + 'E' + ENERGY + '.tot'
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )
            
class pe_CON( OnePerArgList, OutputFile ) :
    """
Class for 'pe.con'.
Contains the total continuum cross section for binary channels,
the particle continuum cross section, and continuum gamma cross section.
INPUTS:
    inelasticPrefix = pe, string containing the projectile (p) and ejectile (e)
                      using the TALYS particle notation npdthag.
KEYS:
    CrossSection
    ContinuumCrossSection
    ContinuumGammaCrossSection
""" 
    def __init__( self, talysZA, inelasticPrefix ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "ContinuumCrossSection", \
                                  "ContinuumGammaCrossSection" )
            self.pe = inelasticPrefix
            self.filename = inelasticPrefix + '.con'
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )

class pe_TOT( OnePerArgList, OutputFile ) :
    """
Class for 'pe.tot'.
Contains the total cross section for binary channels,
the particle discrete cross section, the particle continuum cross section,
and the total gamma cross section.
INPUTS:
    inelasticPrefix = pe, string containing the projectile (p) and ejectile (e)
                      using the TALYS particle notation npdthag.
KEYS:
    CrossSection
    DiscreteCrossSection
    ContinuumCrossSection
    GammaCrossSection
""" 
    def __init__( self, talysZA, inelasticPrefix ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "DiscreteCrossSection", \
                                  "ContinuumCrossSection", \
                                  "GammaCrossSection" )
            self.pe = inelasticPrefix
            self.filename = inelasticPrefix + '.tot'
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'energy', self.columnLabels )

class pe_LXX( OnePerArgList, OutputFile ) :
    """
Class for 'pe.L??'.
Contains the total cross section for binary channels to discrete levels,
the direct cross section and the compound cross section.
The gamma multiplicity is also calculated from the level scheme for each level.
INPUTS:
    inelasticPrefix = pe, string containing the projectile (p) and ejectile (e)
                      using the TALYS particle notation npdthag.
    level           = level number of discrete state.
KEYS:
    CrossSection
    DirectCrossSection
    CompoundCrossSection
    GammaMultiplicity
""" 
    def __init__( self, talysZA, inelasticPrefix, level ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix, level )
        
        if self.newInstance :
        
            self.columnLabels = ( "CrossSection", \
                                  "DirectCrossSection", \
                                  "CompoundCrossSection" )
            self.pe = inelasticPrefix
            self.filename = inelasticPrefix + '.L' + str( level ).zfill( 2 )
            #print 'new pe_LXX:',self.filename
            self.level = int( level )
            OutputFile.__init__( self, talysZA )
            
            #print 'new pe_LXX:',level
            gl = XSnpdtha_GAM( talysZA, talysMisc.particleListFromE( self.pe[1] ) ).gammaList
            #print 'pe =',self.pe,'type(gammaList) =',type(gl),'len(gammalist) =',len(gl)
            #print 'gammalist =',gl
            #try:
            #print type(gl),type(gl[-1])
            
            
            #print 'gl type and size:',type(gl), len(gl)
                
            try:
                #print 'gl type and size:',type(gl[-1]), len(gl[-1])
                self.gammaList = copy.copy( gl[-1] )
            except:
                #print 'gl type and size:',gl,type(gl), len(gl)
                self.gammaList = copy.copy( gl )
                #self.gammaList = [0.0,0.0]
                
                
            ##except:
            #    self.gammaList = []
            #    print "no gammalist"    
            
            self.levelScheme = talysStructure.levels().levelScheme( talysZA.Z, talysZA.A )
            self.set2dData( 'energy', self.columnLabels )
            self.X1 = pe_CON( talysZA, self.pe ).Q - self.Q

    def calculateMultiplicity( self ) :
        """
Redefines the calculateMultiplicity function in process to calculate the gamma multiplicity
from the function Mi.
"""
        label = 'GammaMultiplicity'
        self[label] = endl2dmath( [], label=label )
        gammas = self.gammaList
        if self.level in [ gamma.level_i for gamma in gammas ] :
            M = self.Mi( self.level )
            self[label].setValue( self.threshold, M )
            self[label].setValue( gammas.energy, M )

    def Mi( self, i ) :
        """
Recursive function to calculate the gamma multiplicity for this level
directly from the level scheme and branching ratios.
"""
        if i == 0 : return 0
        else :
            M = 1
            for f in self.levelScheme.levels[i].branches.keys() :
                M += self.levelScheme.levels[i].branches[f].BR * self.Mi( f )
            return M

class decay_e( OnePerArgList, OutputFile ) :
    """
Class for 'decay.e'.
Contains the discrete gamma multiplicities for the binary channels (where e donated the channel).
INPUTS:
    ejectile = e, string containing the ejectile (e=npdthag)
KEYS:
    GammaMultiplicity
"""
    def __init__( self, talysZA, ejectile ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, ejectile )
        if self.newInstance :
            self.filename = 'decay.' + ejectile
            OutputFile.__init__( self, talysZA )
            
class pe_energyANG_LXX( OnePerArgList, OutputFile ) :
    """
Class for 'pe???.???ANG.L??'.
Contains the angular distributions corresponding the pe.L?? cross sections.
INPUTS:
    inelasticPrefix = pe, string containing the projectile (p) and ejectile (e)
                      using the TALYS particle notation npdthag.
    energy          = incident energy
    level           = level number of discrete state.
KEYS:
    CrossSection
    DirectCrossSection
    CompoundCrossSection
""" 
    normalize = True

    def __init__( self, talysZA, inelasticPrefix, energy, level ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix, energy, level )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "DirectCrossSection", \
                                  "CompoundCrossSection" )
            self.energy = float( energy )
            self.level = int( level )
            ENERGY = ("%.3f" % self.energy ).zfill(7)
            self.filename = inelasticPrefix + ENERGY + 'ang.L' + str( level ).zfill( 2 )
            self.pe = inelasticPrefix
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'angle', self.columnLabels )

            
class pe_energyLEG_LXX( OnePerArgList, OutputFile ) :
    """
Class for 'pe???.???LEG.L??'.
Contains the legendre coefficients corresponding the pe.L?? cross sections.
INPUTS:
    inelasticPrefix = pe, string containing the projectile (p) and ejectile (e)
                      using the TALYS particle notation npdthag.
    energy          = incident energy
    level           = level number of discrete state.
KEYS:
    CrossSection
    DirectCrossSection
    CompoundCrossSection
""" 
    normalize = False
    
    def __init__( self, talysZA, inelasticPrefix, energy, level ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, inelasticPrefix, energy, level )
        if self.newInstance :
            self.columnLabels = ( "CrossSection", \
                                  "DirectCrossSection", \
                                  "CompoundCrossSection", \
                                  "Normalized" )
            self.energy = float( energy )
            self.level = int( level )
            ENERGY = ("%.3f" % self.energy ).zfill(7)
            self.filename = inelasticPrefix + ENERGY + 'leg.L' + str( level ).zfill( 2 )
            self.pe = inelasticPrefix
            OutputFile.__init__( self, talysZA )
            self.set2dData( 'coeff', self.columnLabels )
            #print "leg data", 

class XSnpdtha_GAM( OnePerArgList, OutputFile ) :
    """
Class for 'xs??????.gam'.
Contains the discrete gammas for exclusive channels.
Should be called from talys3d.DiscreteEnergySpectra( talysZA, particleList ).
INPUTS:
    particleList = particle list of the form npdtha where the number
                   describes the number of particles of each type.
KEYS:
    DiscreteCrossSection
DATA:
    gammaList = list of __GamDisList objects, one for each incident energy.
""" 
    
    class __GamDisList( list ) :
        """
Class contains a list of discrete gamma rays, all for single incident energy, that come listed after
Einc ngam in xs*.gam file
"""
        class __GamDis :
            """Class contains a single discrete gamma ray, info obtained from single dataline in xs*.gam file"""
            def __init__( self, data ) :
                self.level_i      = Counter( data[0] ).getValue()
                self.level_f      = Counter( data[1] ).getValue()
                self.CrossSection = CrossSection( data[2] ).getValue()
                self.X1_i         = Energy( data[3] ).getValue()
                self.X1_f         = Energy( data[4] ).getValue()
                self.energy       = float( str( self.X1_i - self.X1_f ) ) # float(str()) forces degenerate gammas to be equal, even if X1_i and X1_f were different

            def __repr__( self ) :
                return str( self.energy ) + " " + str( self.CrossSection )

        def __init__( self, datalines, energy ) :
            list.__init__( self )
            self.energy = energy
            self += map( self.__GamDis, datalines )
            self.CrossSection = sum( [ gamma.CrossSection for gamma in self ] )

        def __repr__( self ) :
            string = str( self.energy ) + " " + str( len( self ) )
            if self : string += '\n'
            string += '\n'.join( map( str, self ) )
            return string

    def __init__( self, talysZA, particleList ) :
        OnePerArgList.__init__( self, talysZA, self.__class__, particleList )
        if self.newInstance :
            self.filename = 'xs' + particleList + '.gam'
            #print 'filename',self.filename
            OutputFile.__init__( self, talysZA )
            self.gammaList = []
            #self.gammaList.append( self.__GamDisList( [0.0], 0.0 ) )
            #print 'numen =',self.numen
            if not hasattr( self, 'numen' ) : return
            #print 'numen',self.numen
            for i in range( self.numen ) :
                #print i
                #if 'xs001000.gam' in self.filename : print self.datalines[0]
                try:
                    datapoints = self.datalines.pop( 0 )
                except:
                    print 'Ran out of file content',self.filename
                    return
                    
                energy = Energy( datapoints[0] ).getValue()
                ngam = Counter( datapoints[1] ).getValue()
                gammalines = []
                for j in range( ngam ) :
                    gammalines.append( self.datalines.pop(0) )
                #print 'gammalines',gammalines    
                self.gammaList.append( self.__GamDisList( gammalines, energy ) )
            #print len(self.gammaList)
            self.process()
            if len(self.gammaList)==0: 
                #print 'no gammaList'
                self.gammaList = []
                return
                #self.gammaList = [self.__GamDisList( [0.0], 0.0 )]
            
            #print self.gammaList[0]
            
    def __repr__( self ) :
        return '\n'.join( map( str, self.gammaList ) )
