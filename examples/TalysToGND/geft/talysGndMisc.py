import fudge, os, glob
from operator import add
from fudge.legacy.endl import *
import talysStructure
#from talysMisc import *
import talysMisc 

# incident names for GND inputs        
GNDincidentnames = {0:'', 1:'n', 2:'H1', 3:'H2', 4:'H3', 5:'He3', 6:'He4', 7:'gamma', 8:'eplus', 9:'eminus', 10:'ec'}
        

reactionSymbolList = ['n', '2n','3n','4n', 'p', 'd', 't', 'h', 'a', 'g',
                      'pd', 'np', 'pn', 'ha', 'ta', 'tp', 'da', 'pa', 'dn', 'nt', 'nh', 'na', 'nd', '2a', '2p',
                      'nta', 'nda', 'npa', '2np', '2npa', '2nd', '2na', 'n2a', '2pa', 'n2p', 'n3a', '3np', 
                       'h2a', '3a', 't2a', 'd2a', 'pt2a', 'p2a', '2p2a', 'np2a', '2n2a', '3na'
                        ]

reactionList = ['total','elastic','fission']+reactionSymbolList
print reactionList


def npdthaFromName(name):
    npdtha = [0]*6
    if name in reactionSymbolList:
        for i,x in enumerate(name):
            if x in 'npdtha':
                if name[i-1] in '123456789':
                    npdtha['npdtha'.find(x)] =  int(name[i-1])
                else:
                    npdtha['npdtha'.find(x)] =  1
    return ''.join([ str(x) for x in npdtha ])        

def npdthaUniqueName(name):
    npdtha = []
    if name in reactionSymbolList:
        for i,x in enumerate(name):
            if x in 'npdtha':
                npdtha.append(x)
    return npdtha

def getGenreName(C):
    #import talysMisc 
    if C in reactionSymbolList:
        npdtha = talysMisc.particleList(C)
        if sum([int(x) for x in npdtha.split()]) == 1:
            return "twoBody"
        else:     
            return "NBody"
    else:
        return "unknown"
    
    
def getMTNumber(C):
    import fudge.legacy.converting.endf_endl as FL

    return 1 # filler for now
            
    
def getMassFromTalys(name):

    Z,A = fudge.particles.nuclear.getZandAFromName(name)
    #print 'name,Z,A:',name,Z,A
    #talysPath =  '/g/g16/jurgenso/cnp_src/talys/trunk/structure/masses/'
    talysPath = os.path.join(talysStructure.talysStructurePath,'masses')
    m = talysStructure.talysMasses( path = talysPath)
    return  str(m.getMass( Z, A ))+" amu"

def getLevelFromTalys(name,level):
    Z,A = fudge.particles.nuclear.getZandAFromName(name)
    L = talysStructure.levels().levelScheme(Z,A).levels[level]
    if L.level != level:
        raise Exception( 'Level number wrong', L.level ,'!=',level )
    else:
        pass
        #print "found level",L.level
    return L
