import fudge, os, glob
from operator import add
from fudge.legacy.endl import *
import talysGndMisc 


debug = False
smallestEnergyForSpectra = 0.001

#grabbed the following 4 from endlNameMap in the XEndl package
#these are the yi (yo) names
incidentnames = ["none", "neutron", "proton", "deuteron", "triton", "helium-3", "alpha", "photon", "positron", "electron","electronCapture",
                          "neutronasRecoil","protonasRecoil","deuteronasRecoil", "tritonasRecoil","helium-3asRecoil","alphaasRecoil","photonasRecoil","positronasRecoil","electronasRecoil"]
#these are the C numbers (hmm 36 n3a was "" ?)
reactionnames = ["", "total", "", "", "", "production", "destruction", "", "lacs", "nuclearInterface", "elastic", "n", "2n", "3n",
        "4n", "fission", "3np", "n2p", "2p", "pd", "np", "pn", "nd", "nda", "nt", "nHe3", "na", "n2a", "nta", "2np", "gna",
        "2npa", "2nd", "2na", "npa", "dn", "n3a", "2a", "He3a", "tp", "p", "d", "t", "ta", "He3", "a", "g", "da", "pa", "2pa", "xp",
        "xd", "xt", "xHe3", "xa", "modify_xg", "xn", "xe", "", "", "", "", "", "", "", "act", "yield", "", "", "", "totp", "coh", "incoh",
        "photo", "pair", "triplet", "", "", "ic", "", "", "ion", "brem", "exit", "coll", "", "", "", "", "", "", "shell",
        "trans", "whole"]
#these are the I numbers
propertynames = {0:"integrated", 1:"angular", 3:"energyAngle", 4:"energyAngleLegendre", 7:"nubar", 8:"energyHistogram", 9:"multiplicity",
        10:"averageEnergyToLightProduct", 11:"averageEnergyToHeavyResidual", 12:"energyAvailable", 20:"energyLoss", 21:"energyDistribution",
        22:"angularDistribution", 30:"stragglingCrossSection", 80:"maxwellianReactionRates", 81:"inFlightCrossSection", 
        84:"energyDistForThermalEquil", 89:"particleMultiplicityForThermalEquil", 90:"maxwellianEnergyToLightProduct", 
        91:"maxwellianEnergyToHeavyResidual", 92:"maxwellianTotalEnergy", 99:"supplementalParameter", 912:"electronPerSubshell", 
        913:"subshellBindingEnergy", 914:"subshellKeneticEnergy", 915:"subshellAverageRadius", 921:"subshellRadiativeLevelWidth", 
        922:"subshellNonradiativeLevelWidth", 931:"subshellRadiativeTransitionProb", 932:"subshellNonradiativeTransitionProb", 
        933:"particlesPerVacancy", 934:"particleEnergyPerVacancy", 935:"energyToResidualAtom", 941:"formFactor",
        942:"scatteringFunction", 943:"imaginaryAnomalousScatteringFactor", 944:"realAnomalousScatteringFactor"}
#these are the S numbers?
modifiernames = {0:"unSpecified", 1:"level", 2:"preEquilibrium", 3:"discreteGamma", 5:"activation", 
        7:"delayed", 8:"cluster", 10:"wideLevelExcitation", 11:"secondNeutron", 15:"densityTemperatureDependent", 
        91:"subshell"}
        
def eFromParticleList( npdtha ) :
    """Return the ejectile from the particleList of numbers npdtha. If more than one particle emitted returns None."""
    sum = reduce( add, map( float, npdtha ), 0 )
    if sum > 1 : return None
    if sum == 0 : return 'g'
    for position, count in enumerate( npdtha ) :
        if count == '1' : return talysYTags( position+1 )

def particleListFromE( e ) :
    """Return the particleList of numbers npdtha from the ejectile label."""
    npdtha = ['0']*6
    if e == 'g' : return ''.join( npdtha )
    for p in range( 6 ) :
        if talysYTags( p+1 ) == e :
            npdtha[p] = '1'
            return ''.join( npdtha )

def isNaN( x ) : # both <= and >= return False when comparing to nan
    """Return True if x is NaN and False if its a real number."""
    if not (x <= 0.) and not (x >= 0.) : return True

def getMultiplicity( GammaCrossSection, CrossSection ) :
    for d in [ GammaCrossSection, CrossSection ] : d.trim()
    fudge.legacy.endl.endl2dmathClasses.doSafeDivide = True
    GammaMultiplicity = GammaCrossSection / CrossSection
    removeNaN( GammaMultiplicity )
    return GammaMultiplicity

def removeNaN( data ) :
    for x,y in data.data :
        if isNaN( y ) : data.set( data.slice( xMin=x+1e-11 ) )

def getFileList( talysZA, globber ) :
    """Returns a filelist from the talysDir found in talysZA for a given glob string."""
    try : os.chdir( os.path.join( talysZA.workDir, talysZA.talysDir ) )
    except : return []
    filelist = glob.glob( globber )
    os.chdir( talysZA.cwd )
    return filelist

def particleName( y ) :
    """Returns the name of the particle used in geft for a given endl_y value."""
    if y == 1 : return 'Neutron'
    if y == 2 : return 'Proton'
    if y == 3 : return 'Deuteron'
    if y == 4 : return 'Triton'
    if y == 5 : return 'Helium3'
    if y == 6 : return 'Alpha'
    if y == 7 : return 'Gamma'

def talysYTags( y ) :
    """Return y tags in talys language, 1234567 = npdthag"""
    if y == 1 : return 'n'
    if y == 2 : return 'p'
    if y == 3 : return 'd'
    if y == 4 : return 't'
    if y == 5 : return 'h'
    if y == 6 : return 'a'
    if y == 7 : return 'g'
    return None

def ZAforYTags( y ):
    """Return ZA for y tags to calculate the products """
   
   ### had leading zeros which made it give octal answers
    
    if y == 1 : return 1
    if y == 2 : return 1001
    if y == 3 : return 1002
    if y == 4 : return 1003
    if y == 5 : return 2003
    if y == 6 : return 2004
    if y == 7 : return 0      # no net nucleons
    if y == 8 : return -1000+1     # positron means you lost a proton to a neutron
    if y == 9 : return +1000-1     # electron out means neutron decay
    return 0


def inelasticPrefix( yi, C ) :
    """Return a string projectile+ejectile for use in TALYS binary reactions."""
    if C in [talysYTags( x+1 ) for x in range(7)]:
        prefix = talysYTags( yi ) + C
    elif C in CList:
        prefix = talysYTags( yi ) + talysYTags( binaryYDict[C] )
    else:
        print 'WARNING unknown channel designator: ',C     
    return prefix

def particleList( C ) :
    """
Returns a string corresponding to the particles emitted for each reaction C
in talys format of 6 numbers 'npdtha' i.e. (x,2nd) 201000
"""
    if C in CList:
        yoList = []
        yos = endl_C.endl_C_yoInfo( C )
        for yo in yos :
            if yo in endl_C.yoDict : yoList.append( endl_C.yoDict[yo] )
        particleList = ''
        for p in range( 6 ) :
            particleList += str( yoList.count( p+1 ) )
        return particleList
    elif C in talysGndMisc.reactionSymbolList:
        return talysGndMisc.npdthaFromName(C)   
    else:
        print 'WARNING unknown channel designator: ',C   
        return ''   

def yoList( C ) :
    """Returns list of yo's (appearing only once) for each reaction C"""
    
    if C in CList:
        yoList = []
        yos = endl_C.endl_C_yoInfo( C )
        for yo in yos :
            if yo != 7 : yop = endl_C.yoDict[yo] # leave out gammas as they are always there and are treated differently
            if yo in endl_C.yoDict and yop not in yoList : yoList.append( yop )
        return yoList
    elif C in talysGndMisc.reactionSymbolList:
        return [ [talysYTags( y+1 ) for y in range(7) ].index(x)+1 for x in talysGndMisc.npdthaUniqueName(C) ]
    else:
       print 'WARNING unknown channel designator: ',C   
       return ''   
  

CList = [1]+range(10,50)
excludeCList = [21,35,30] # these channels are duplicates but different ordering which are summed together in TALYS, or C=30 is C=26 with gammas
CList = filter( lambda C: C not in excludeCList, CList )
binaryCList = [11,40,41,42,44,45,46]
binaryYDict = { 11:1, 40:2, 41:3, 42:4, 44:5, 45:6, 46:7 }

binaryPartList = [talysYTags( y ) for y in range(7)]


def getBinaryCList( yi, S1Data=False, onlyInelasticDiscrete=True ) : #note at the moment don't change from default, until levelScheme can get all level schemes in discreteLevelEnergySpectra
    if onlyInelasticDiscrete == False : return binaryCList
    for C in binaryCList :
        if binaryYDict[C] == yi :
            for Cp in binaryCList :
                if reactionnames[Cp][-1] == "'" : reactionnames[Cp] = reactionnames[Cp][0:-1] 
            reactionnames[C] += "'"
            if not S1Data : return []
            return [C]


