import os
import talysMisc, talysStructure
from fudge.legacy.endl import endl_Z

class Input( list ) :
    """
Base class for TALYS input files.
All custom talys input files should inheret from this class,
and then add inputlines as InputLine classes.
Knows how to add new input lines or groups of input lines together.
"""
    def __init__( self ) :
        list.__init__( self )
    
    def __add__( self, other ) :
        if isinstance( other, InputLine ) :
            self.append( other )
            return self
        elif isinstance( other, Input ) :
            self.append( '' )
            for i in other :
                self.append( i )
            return self
        else : raise TypeError
   
    def __iadd__( self, other ) :
        return self.__add__( other )
    
    def __repr__( self ):
        return '\n'.join( map( str, self ) )+'\n'

    def toFile( self, dirname='.', filename='input' ) :
        """Writes the input file to dirname/filename. If it exists raises IOError."""
        name = os.path.join( dirname, filename )
        if os.path.isfile( name ) : raise IOError
        f = open( str( name ), 'w' )
        f.write( str( self ) )
        f.close()

    def addKeyword( self, *args, **kw ) :
        try : 
            for keyword, value in args :
                if not self.fissionBarrier( keyword, value ) :
                    self += InputLine( keyword, value )
        except :
            for keyword, value in args.iteritems() :
                if not self.fissionBarrier( keyword, value ) :
                    self += InputLine( keyword, value )
        for keyword in kw :
            if not self.fissionBarrier( keyword, kw[keyword] ) :
                self += InputLine( keyword, kw[keyword] )

    def fissionBarrier( self, keyword, value ) :
        barrierParameters = ['fisbar','fishw','lambda','pshift','deltaW','Nlow','Ntop','Exmatch',\
                             'T','E0','axtype','Rtransmom','Rclass2mom','class2width',]
        keywords = keyword.split()
        for p in barrierParameters:
            if p in keywords[0] and keywords[0] != p :
                barrier = keywords[0].lstrip( p )
                keywords[0] = keywords[0].rstrip( barrier )
                keyword = ' '.join( keywords )
                self += InputLine( keyword, value, barrier )
                return True

class InputLine:
    """
Class contains the information for one line of an TALYS input file.
name = talys keyword
If no additional arguments are given then 'y' is assumed,
otherwise the additional arguments are added to the line.
"""
    def __init__( self, name, *args ):
        self.inputline = [ str( name ) ]
        if not args :
            self.inputline.append( 'y' )
        for value in args :
            self.inputline.append( str( value ) )

    def __repr__( self ):
        return self.inputline[0].ljust(max(14,len(self.inputline[0])+1))+' '.join( self.inputline[1:] )

class InputHeader( Input ) :
    """
Contructs the TALYS default header that provides basic reaction
information that cannot be omitted.
talysZA is passed as an argument and the information extracted from it.
"""
    def __init__( self, talysZA ) :
        Input.__init__( self )
        self += InputLine( 'projectile', talysMisc.talysYTags( talysZA.yi ) )
        self += InputLine( 'element', endl_Z.endl_ZSymbol( talysZA.Z ) )
        self += InputLine( 'mass', talysZA.A )
        self += InputLine( 'energy', 'energies' )
        if talysZA.suffix :
            levelScheme = talysStructure.levels().levelScheme( talysZA.Z, talysZA.A )
            for level,l in enumerate( levelScheme.levels ) :
                if l.energy == talysZA.targetELevel : break
            else :
                raise RuntimeError, 'Isomer level energy %s not found in level scheme.' % talysZA.targetELevel
            self += InputLine( 'Ltarget', level )

class OutputFlags( Input ) :
    """Class gives the full list of output flags needed to provide a complete evaluation."""
    def __init__( self, **kw ) :
        Input.__init__( self )
        self += InputLine( '# Output Flags', '' )
        outputFlags = ()
        if kw.pop('endf',None) : outputFlags += ( 'endf', 'endfdetail' )
        if kw.pop('channels',None) : outputFlags += ( 'channels', 'filechannels' )
        if kw.pop('elastic',None) : outputFlags += ( 'outangle', 'fileelastic', 'outdiscrete' )
        if kw.pop('spectra',None) :
            if 'channels' not in outputFlags : outputFlags += ( 'channels', )
            outputFlags += ( 'outspectra', 'filespectrum' )
        if kw.pop('gammas',None) : outputFlags += ( 'outgamdis', 'filegamdis' )
        for i in outputFlags : self += InputLine( i )
        if 'endf' not in outputFlags and 'outspectra' in outputFlags : self.addKeyword( adddiscrete='n', addelastic='n' )
        if kw.pop('discrete',None) :
            if 'outangle' not in outputFlags :
                self.addKeyword( outangle='y', outdiscrete='y' )
            for i in range( 21 ) :
                self += InputLine( 'filediscrete', i )
                self += InputLine( 'fileangle', i )

class energyList( Input ) :
    """
Class to construct a energy list which can be read into TALYS.
Default is to give a list which starts with 0.001 MeV,
has rough log dependence upto 0.1 MeV, then steps of 0.2 MeV upto 20 MeV.
Can pass a custom energyList, or specifiy the start and end and step as keywords.
If endf=True then start energyList at 10e-11 and increase in factors of 10 upto 10e-3.
"""
    def __init__( self, energyList=None, start=0.001, step=0.2, end=20., endf=False ) :
        Input.__init__( self )
        lowEnergyList = [ 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1 ]
        if endf :
            start = 0
            endfLowEnergies = map(lambda x: pow(10,-x),range(11,3,-1))
            lowEnergyList = endfLowEnergies + lowEnergyList
        if energyList == None :
            number = int(end/step)
            energyList = [energy for energy in lowEnergyList if energy >= start]
            for energy in map( lambda x: float(x+1)*step, range( number ) ) :
                if energy >= start : energyList.append( energy )
        for energy in energyList : self += InputLine( energy, '' )

    def toFile( self, dirname='.', filename='energies' ) :
        """Writes the energylist to dirname/filename. If it exists raises IOError."""
        name = os.path.join( dirname, filename )
        if os.path.isfile( name ) : raise IOError
        f = open( str( name ), 'w' )
        f.write( str( self ) )
        f.close()

