from collections import namedtuple
import unittest, copy

# ------------------------------------------------
# Simple C4 containers
# ------------------------------------------------

C4Point = namedtuple( 'C4Point', 'projectile target targetMetastableState MF MT productMetastableState status cmFlag energy  dEnergy  data  dData   cosMuOrLegendreOrder   dCosMuOrLegendreOrder   eLevelOrHalflife  dELevelOrHalflife idOf78 reference exforEntry exforSubEntry multiDimFlag' )

C4DataSet = namedtuple( 'C4DataSet', 'dataSet date reaction projectile target MF MT c4Begin numData data' )

C4Entry = namedtuple( 'C4Entry', 'entry author1 year institute title authors refCode reference numDataSets dataSets' )

# ------------------------------------------------
# C4 Parsers
# ------------------------------------------------
def emptyStringToNone( s ): 
    if s.strip() == '': return None
    return s
    
def NoneToEmptyString( n ):
    if n == None: return " "
    return n
    
def readFunkyFloat( value ):
    value = emptyStringToNone( value )
    if value == None: return
    try :
        f = float( value )
    except :
        value = value.replace( ' ', '' )
        i = value.find( '.' )
        if( i != -1 ) :
            j = value[i:].find( '+' )
            if( j == -1 ) : j = value[i:].find( '-' )
            if( j == -1 ) : raise ValueError( "Float value %s is not funky enough <%s>" % ( value ) )
            value = value[:i+j] + 'e' + value[i+j:]
        try :
            f = float( value )
        except :
            if( value == len( value ) * ' ' ) : return( 0 )
            raise ValueError( 'Could not convert value "%s"' % ( value ) )
    return( f )

def writeFunkyFloat( f, width=11, significant_digits=None ): 
    '''
    Write one of Red's Funky Floats...
    
    Anatomy of a funky float:
    
    Consider the number 3124.5611 written to a field with width=10
        ' 3.12456+3'
        '_X.XXXXX+E'
    where '_' is either a space (' ') or a minus sign ('-'), 
    'X.XXXXX' is the number itself, in exponential notation, with exponent 'E'.  
    The floating point 'E' is dropped, saving a character.
    '''
    if f == None: return width*' '
    
    # compute the number of sig figs the user "really" wants
    if significant_digits == None or significant_digits > width-1: significant_digits=width-1
    
    # compute the possible output formats
    g_frmt='%'+str(width)+'.'+str(significant_digits)+'g'
    e_frmt='%'+str(width)+'.'+str(significant_digits-4)+'e'
    
 #   print g_frmt, e_frmt
    
    # try with the 'g' format, we may get lucky
    s = (g_frmt) % round(f,significant_digits)
    
    # OK, number too big, so must use explicit exponential notation, but with E-less representation
    if 'e' in s or len(s) > width:
        s = (e_frmt) % round(f,significant_digits)
        s = s.replace( 'e', '' )
    
    # Check to make sure the number isn't too long
    if len(s) > width: print 'too big, len(%s)=%i<%i'%(s,len(s),width)
    return s.strip().ljust(width)

def readC4File( fList, asPointList = True ):
    if asPointList: return map( readC4Point, filter( lambda x: not ( x.startswith( '#' ) or x.strip()=='' ), fList ) )
    newList = []
    for line in fList:
        if line.strip() in [ '', '#' ]: pass
        elif line.startswith( '#ENTRY' ): newList.append( [ line ] )
        else: newList[-1].append( line )
    return map( readC4Entry, newList )

def readC4DataSet( fList ): 
    '''
    An example dataset header taken from a c4 file:

        #DATASET    40617007
        #DATE       19850305
        #REACTION   40-ZR-92(N,TOT),,SIG
        #PROJ       1
        #TARG       40092
        #MF         3
        #MT         1
        #C4BEGIN    [    1 40092   3   1 A ]
        #DATA       4
        # Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
        #---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
            1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
        #/DATA      4
        #/DATASET
    '''    
    dataSet     =   ''
    date        =   ''
    reaction    =   ''
    projectile  =   0
    target      =   0
    MF          =   0
    MT          =   0
    c4Begin     =   '' 
    numData     =   0
    data        =   []
    flist = copy.copy( fList )
    line = flist.pop(0)
    while flist:
        if line.startswith("#DATASET"):     dataSet = line[12:].strip()
        elif line.startswith("#DATE"):      date = line[12:].strip()
        elif line.startswith("#REACTION"):  reaction = line[12:].strip()
        elif line.startswith("#PROJ"):      projectile = int( line[12:].strip() )
        elif line.startswith("#TARG"):      target = int( line[12:].strip() )
        elif line.startswith("#MF"):        MF = int( line[12:].strip() )
        elif line.startswith("#MT"):        MT = int( line[12:].strip() )
        elif line.startswith("#C4BEGIN"):   c4Begin = line[12:].strip()
        elif line.startswith("#DATA"):   
            numData = int( line[12:].strip() )
            line = flist.pop(0)
            line = flist.pop(0)
            line = flist.pop(0)
            while not line.startswith( "#/DATA" ): 
                data.append( readC4Point( line ) )
                line = flist.pop(0)
        elif line.strip() == "#": pass
        elif line.startswith( "#/DATASET" ): pass
        else: pass
        line = flist.pop(0)
    return C4DataSet( dataSet=dataSet, date=date, reaction=reaction, projectile=projectile, target=target, MF=MF, MT=MT, c4Begin=c4Begin, numData=numData, data=data )

def readC4Entry( fList ):
    '''
    An example header taken from a c4 file:

        #ENTRY      40617
        #AUTHOR1    M.V.Pasechnik+
        #YEAR       1980
        #INSTITUTE  (4CCPIJI)
        #TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
        #+          AND ZYRCONIUM AT LOW ENERGIES
        #AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
        #+          G.A.Smetanin, T.I.Jakovenko
        #REF-CODE   (C,80KIEV,1,304,8009)
        #REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
        #+          Vol.1, p.304, 1980
        #DATASETS   7
        #
        #DATASET    40617007
        #DATE       19850305
        #REACTION   40-ZR-92(N,TOT),,SIG
        #PROJ       1
        #TARG       40092
        #MF         3
        #MT         1
        #C4BEGIN    [    1 40092   3   1 A ]
        #DATA       4
        # Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
        #---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
            1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
            1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
        #/DATA      4
        #/DATASET
        #/ENTRY

    '''
    entry       =   ''
    author1     =   'AUTHOR1' 
    year        =   0
    institute   =   ''
    title       =   ''
    authors     =   ''
    refCode     =   ''
    reference   =   ''
    numDataSets =   0
    dataSets    =   []
    flist = copy.copy( fList )
    line = flist.pop(0)
    while flist:
        if line.startswith("#ENTRY"):       entry = line[12:].strip()
        elif line.startswith("#AUTHOR1"):   author1 = line[12:].strip()
        elif line.startswith("#YEAR"):      year = int( line[12:].strip() )
        elif line.startswith("#INSTITUTE"): institute = line[12:].strip()
        elif line.startswith("#TITLE"):     
            title = line[12:].strip() 
            while flist[0].startswith( '#+' ):
                line = flist.pop( 0 )
                title += ' ' + line[12:].strip()
        elif line.startswith("#AUTHOR(S)"): 
            authors = line[12:].strip()
            while flist[0].startswith( '#+' ):
                line = flist.pop( 0 )
                authors += ' ' + line[12:].strip()
        elif line.startswith("#REF-CODE"):  refCode = line[12:].strip()
        elif line.startswith("#REFERENCE"): 
            reference = line[12:].strip()
            while flist[0].startswith( '#+' ):
                line = flist.pop( 0 )
                reference += ' ' + line[12:].strip()
        elif line.startswith("#DATASETS"):  numDataSets = int( line[12:].strip() )
        elif line.startswith("#DATASET"):   
            sublist = [ line ]
            while not flist[0].startswith( '#/DATASET' ):
                sublist.append( flist.pop( 0 ) )
            sublist.append( flist.pop(0) )
            dataSets.append(  readC4DataSet( sublist ) )
        elif line.strip() == "#": pass  
        elif line.startswith( "/ENTRY" ): pass  
        else: pass
        line = flist.pop(0)
    return C4Entry( 
        entry = entry,
        author1 = author1,
        year = year,
        institute = institute,
        title = title,
        authors = authors,
        refCode = refCode,
        reference = reference,
        numDataSets = numDataSets,
        dataSets = dataSets)

def writeC4Point( pt ):
    return \
        str(pt.projectile).rjust(5) + \
        str(pt.target).rjust(6) +  \
        str(NoneToEmptyString(pt.targetMetastableState))[0] + \
        str(pt.MF).rjust(3) + \
        str(pt.MT).rjust(4) + \
        str(NoneToEmptyString(pt.productMetastableState))[0] + \
        str(NoneToEmptyString(pt.status))[0] + \
        str(NoneToEmptyString(pt.cmFlag))[0] + \
        writeFunkyFloat(pt.energy,9) + \
        writeFunkyFloat(pt.dEnergy,9) + \
        writeFunkyFloat(pt.data,9) + \
        writeFunkyFloat(pt.dData,9) + \
        writeFunkyFloat(pt.cosMuOrLegendreOrder,9) + \
        writeFunkyFloat(pt.dCosMuOrLegendreOrder,9) + \
        writeFunkyFloat(pt.eLevelOrHalflife,9) + \
        writeFunkyFloat(pt.dELevelOrHalflife,9) + \
        str(NoneToEmptyString(pt.idOf78)).rjust(3) + \
        str(pt.reference).ljust(25) + \
        str(pt.exforEntry).rjust(4) + \
        str(pt.exforSubEntry).rjust(3) + \
        str(NoneToEmptyString(pt.multiDimFlag))

def readC4Point( fline ):
    '''
    Here is an example taken from a c4 file of a list of points:

        # Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
        #---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
            1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
            1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
            1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20
            .
            .
            .
    
    Note, the PXC field is really 3 one character fields.  From the x4toc4 manual: 
    
          Columns   Description
          -------   -----------
            1-  5   Projectile ZA (e.g. neutron =1, proton =1001)
                    (defined by reaction dictionary).
            6- 11   Target ZA (e.g. 26-Fe-56 =  26056)
                    (defined by EXFOR reaction).
               12   Target metastable state (e.g. 26-FE-56m = M)
                    (defined by EXFOR reaction).
           13- 15   MF (ENDF conventions, plus additions)
                    (defined by reaction dictionary).
           16- 19   MT (ENDF conventions, plus additions)
                    (defined by reaction dictionary).
               20   Product metastable state (e.g. 26-FE-56M = M)
                    (defined by EXFOR reaction).
               21   EXFOR status
                    (defined by EXFOR keyword status).
               22   Center-of-mass flag (C=center-of-mass, blank=lab)
                    (defined by EXFOR title dictionary).
           23- 94   8 data fields (each in E9.3 format defined below)
                    (defined by MF and title dictionary).
           95- 97   Identification of data fields 7 and 8
                    (e.g., LVL=level, HL=half-life, etc.).
                    For a complete list of codes see title dictionary
                    (defined by MF and title dictionary).
           98-122   Reference (first author and year)
                    (defined by EXFOR keywords title and reference).
          123-127   EXFOR accession number
                    (defined by EXFOR format).
          128-130   EXFOR sub-accession number
                    (defined by EXFOR format).
              131   Multi-dimension table flag
                    (defined by EXFOR keyword reaction or common fields).        
    '''     
    try:               multiDimFlag = emptyStringToNone( fline[130] )
    except IndexError: multiDimFlag = None
    return C4Point( 
        projectile              =   int(fline[0:5]), 
        target                  =   int(fline[5:11]), 
        targetMetastableState   =   emptyStringToNone( fline[11] ), 
        MF                      =   int(fline[11:15]), 
        MT                      =   int(fline[15:19]), 
        productMetastableState  =   emptyStringToNone( fline[19] ), 
        status                  =   emptyStringToNone( fline[20] ), 
        cmFlag                  =   emptyStringToNone( fline[21] ), 
        energy                  =   readFunkyFloat( fline[22:31] ), 
        dEnergy                 =   readFunkyFloat( fline[31:40] ), 
        data                    =   readFunkyFloat( fline[40:49] ), 
        dData                   =   readFunkyFloat( fline[49:58] ), 
        cosMuOrLegendreOrder    =   readFunkyFloat( fline[58:67] ), 
        dCosMuOrLegendreOrder   =   readFunkyFloat( fline[67:76] ), 
        eLevelOrHalflife        =   readFunkyFloat( fline[76:85] ), 
        dELevelOrHalflife       =   readFunkyFloat( fline[85:94] ), 
        idOf78                  =   emptyStringToNone( fline[94:97] ), 
        reference               =   fline[97:122].strip(), 
        exforEntry              =   fline[122:127], 
        exforSubEntry           =   int(fline[127:130]), 
        multiDimFlag            =   multiDimFlag )
  
  
  
        
# ------------------------------------------------
# Unit tests
# ------------------------------------------------

class TestFunkyFloats( unittest.TestCase ):
    def setUp(self):
        self.nums={ 0.314159e14:"3.141590+13", -3.14159e13:"-3.141590+13", 51450.00:'51450      ', -51450.00:'-51450     ', 0.00:'0          ', 504400.:"504400     ", 4961.493:"4961.493   ", 8.146:"8.146      ", 0.4042:"0.4042     " }
    def test_write_read_consistency( self ):
        for x in self.nums:
            a = writeFunkyFloat( x )
            b = readFunkyFloat( a )
            self.assertEqual( x, b )
    def test_write( self ):
        for k,v in self.nums.items():
            self.assertEqual( v, writeFunkyFloat( k ) )

class TestReadC4Point( unittest.TestCase ):
    def setUp( self ): 
        self.a = C4Point( projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None, status='A', cmFlag=None, energy=504400.0, dEnergy=4961.493, data=8.146000, dData=0.404200, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='L.GREEN,ET.AL. (73)', exforEntry='10225', exforSubEntry=20, multiDimFlag=None )
    def test_a( self ): 
        txt = '    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20 '
        b = readC4Point( txt ) 
        self.assertEqual( self.a, b )
        c = writeC4Point( b )
        self.assertEqual( txt, c )
        self.assertEqual( self.a, readC4Point(c) )

class TestReadC4Entry( unittest.TestCase ):
    def setUp( self ): 
        self.a = C4Entry( 
            entry='40617', author1 = 'M.V.Pasechnik+', year = 1980, institute = '(4CCPIJI)', 
            title = "TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM AND ZYRCONIUM AT LOW ENERGIES", 
            authors = "M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko, G.A.Smetanin, T.I.Jakovenko",
            refCode = "(C,80KIEV,1,304,8009)",
            reference = "Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980 Vol.1, p.304, 1980",
            numDataSets = 7, 
            dataSets =[ 
               readC4DataSet( '''#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET'''.split('\n') ) ] )
    def test_a( self ): 
        b = readC4Entry( '''#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
'''.split('\n') )
        self.assertEqual( self.a, b )
     
class TestReadC4DataSet( unittest.TestCase ):
    def setUp( self ): 
        self.a = C4DataSet( 
            dataSet = '40617007', date = '19850305', 
            reaction = '40-ZR-92(N,TOT),,SIG', projectile = 1, target = 40092, MF = 3, MT = 1, 
            c4Begin = "[    1 40092   3   1 A ]", numData = 4,
            data = [
                C4Point( projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None, status='A', cmFlag=None, energy=442000.0, dEnergy=None, data=12.74000, dData=1.700000, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617', exforSubEntry=7, multiDimFlag=None ),
                C4Point( projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None, status='A', cmFlag=None, energy=507000.0, dEnergy=None, data=8.790000, dData=0.570000, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617', exforSubEntry=7, multiDimFlag=None ),
                C4Point( projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None, status='A', cmFlag=None, energy=572000.0, dEnergy=None, data=9.520000, dData=0.200000, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617', exforSubEntry=7, multiDimFlag=None ),
                C4Point( projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None, status='A', cmFlag=None, energy=637000.0, dEnergy=None, data=8.480000, dData=0.110000, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617', exforSubEntry=7, multiDimFlag=None ) ] )
    def test_a( self ): 
        b = readC4DataSet( '''#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
'''.split('\n') )
     
class TestReadC4File( unittest.TestCase ):
    def setUp( self ): 
        self.testData = '''
#
#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
#
#
#ENTRY      10225
#AUTHOR1    L.Green+
#YEAR       1973
#INSTITUTE  (1USABET)
#TITLE      Total cross section measurements with a 252Cf time-of-
#+          flight spectrometer.
#AUTHOR(S)  L.Green, J.A.Mitchell
#REF-CODE   (R,WAPD-TM-1073,197304)
#REFERENCE  Rept. Westinghouse Atomic Power Div.(Bettis) Reports
#+          No.1073, 1973
#DATASETS   27
#
#DATASET    10225020
#DATE       20010305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       3
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20
#/DATA      3
#/DATASET
#/ENTRY
#
'''
        self.a = [ 
            readC4Entry('''
#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
            '''.split('\n') ), 
            readC4Entry('''
#ENTRY      10225
#AUTHOR1    L.Green+
#YEAR       1973
#INSTITUTE  (1USABET)
#TITLE      Total cross section measurements with a 252Cf time-of-
#+          flight spectrometer.
#AUTHOR(S)  L.Green, J.A.Mitchell
#REF-CODE   (R,WAPD-TM-1073,197304)
#REFERENCE  Rept. Westinghouse Atomic Power Div.(Bettis) Reports
#+          No.1073, 1973
#DATASETS   27
#
#DATASET    10225020
#DATE       20010305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       3
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20
#/DATA      4
#/DATASET
#/ENTRY
            '''.split('\n') ) ]
        self.b = map( readC4Point, '''    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20'''.split('\n') )
    def test_a( self ): 
        self.assertEqual( readC4File( self.testData.split( '\n' ), asPointList = False ), self.a ) 
    def test_b( self ): 
        self.assertEqual( readC4File( self.testData.split( '\n' ), asPointList = True ), self.b ) 
     
# ------------------------------------------------
# Main !!
# ------------------------------------------------

if __name__ == "__main__": 

    try:
        import xmlrunner
        unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-results'))
    except ImportError:
        unittest.main()
        print 
    