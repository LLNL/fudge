# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Convert between Z, symbol and element name.
adapted from endl_Z
cmattoon, March 2011
"""

# how many Zs are tabulated?
nZs = 119

ZLabels = (
    (   0, "n",  "Neutron" ),
    (   1, "H",  "Hydrogen" ),
    (   2, "He", "Helium" ),
    (   3, "Li", "Lithium" ),
    (   4, "Be", "Beryllium" ),
    (   5, "B",  "Boron" ),
    (   6, "C",  "Carbon" ),
    (   7, "N",  "Nitrogen" ),
    (   8, "O",  "Oxygen" ),
    (   9, "F",  "Fluorine" ),
    (  10, "Ne", "Neon" ),
    (  11, "Na", "Sodium" ),
    (  12, "Mg", "Magnesium" ),
    (  13, "Al", "Aluminium" ),
    (  14, "Si", "Silicon" ),
    (  15, "P",  "Phosphorus" ),
    (  16, "S",  "Sulphur" ),
    (  17, "Cl", "Chlorine" ),
    (  18, "Ar", "Argon" ),
    (  19, "K",  "Potassium" ),
    (  20, "Ca", "Calcium" ),
    (  21, "Sc", "Scandium" ),
    (  22, "Ti", "Titanium" ),
    (  23, "V",  "Vanadium" ),
    (  24, "Cr", "Chromium" ),
    (  25, "Mn", "Manganese" ),
    (  26, "Fe", "Iron" ),
    (  27, "Co", "Cobalt" ),
    (  28, "Ni", "Nickel" ),
    (  29, "Cu", "Copper" ),
    (  30, "Zn", "Zinc" ),
    (  31, "Ga", "Gallium" ),
    (  32, "Ge", "Germanium" ),
    (  33, "As", "Arsenic" ),
    (  34, "Se", "Selenium" ),
    (  35, "Br", "Bromine" ),
    (  36, "Kr", "Krypton" ),
    (  37, "Rb", "Rubidium" ),
    (  38, "Sr", "Strontium" ),
    (  39, "Y",  "Yttrium" ),
    (  40, "Zr", "Zirconium" ),
    (  41, "Nb", "Niobium" ),
    (  42, "Mo", "Molybdenum" ),
    (  43, "Tc", "Technetium" ),
    (  44, "Ru", "Ruthenium" ),
    (  45, "Rh", "Rhodium" ),
    (  46, "Pd", "Palladium" ),
    (  47, "Ag", "Silver" ),
    (  48, "Cd", "Cadmium" ),
    (  49, "In", "Indium" ),
    (  50, "Sn", "Tin" ),
    (  51, "Sb", "Antimony" ),
    (  52, "Te", "Tellurium" ),
    (  53, "I",  "Iodine" ),
    (  54, "Xe", "Xenon" ),
    (  55, "Cs", "Cesium" ),
    (  56, "Ba", "Barium" ),
    (  57, "La", "Lanthanum" ),
    (  58, "Ce", "Cerium" ),
    (  59, "Pr", "Praseodymium" ),
    (  60, "Nd", "Neodymium" ),
    (  61, "Pm", "Promethium" ),
    (  62, "Sm", "Samarium" ),
    (  63, "Eu", "Europium" ),
    (  64, "Gd", "Gadolinium" ),
    (  65, "Tb", "Terbium" ),
    (  66, "Dy", "Dysprosium" ),
    (  67, "Ho", "Holmium" ),
    (  68, "Er", "Erbium" ),
    (  69, "Tm", "Thulium" ),
    (  70, "Yb", "Ytterbium" ),
    (  71, "Lu", "Lutetium" ),
    (  72, "Hf", "Hafnium" ),
    (  73, "Ta", "Tantalum" ),
    (  74, "W",  "Tungsten" ),
    (  75, "Re", "Rhenium" ),
    (  76, "Os", "Osmium" ),
    (  77, "Ir", "Iridium" ),
    (  78, "Pt", "Platinum" ),
    (  79, "Au", "Gold" ),
    (  80, "Hg", "Mercury" ),
    (  81, "Tl", "Thallium" ),
    (  82, "Pb", "Lead" ),
    (  83, "Bi", "Bismuth" ),
    (  84, "Po", "Polonium" ),
    (  85, "At", "Astatine" ),
    (  86, "Rn", "Radon" ),
    (  87, "Fr", "Francium" ),
    (  88, "Ra", "Radium" ),
    (  89, "Ac", "Actinium" ),
    (  90, "Th", "Thorium" ),
    (  91, "Pa", "Protactinium" ),
    (  92, "U",  "Uranium" ),
    (  93, "Np", "Neptunium" ),
    (  94, "Pu", "Plutonium" ),
    (  95, "Am", "Americium" ),
    (  96, "Cm", "Curium" ),
    (  97, "Bk", "Berkelium" ),
    (  98, "Cf", "Californium" ),
    (  99, "Es", "Einsteinium" ),
    ( 100, "Fm", "Fermium" ),
    ( 101, "Md", "Mendelevium" ),
    ( 102, "No", "Nobelium" ),
    ( 103, "Lr", "Lawrencium" ),
    ( 104, "Rf", "Rutherfordium" ),
    ( 105, "Db", "Dubnium" ),
    ( 106, "Sg", "Seaborgium" ),
    ( 107, "Bh", "Bohrium" ),
    ( 108, "Hs", "Hassium" ),
    ( 109, "Mt", "Meitnerium" ),
    ( 110, "Ds", "Darmstadtium" ),
    ( 111, "Rg", "Roentgenium" ),
    ( 112, "Cn", "Copernicium" ),
    ( 113, "Uut", "Ununtrium" ),
    ( 114, "Fl",  "Flerovium" ),
    ( 115, "Uup", "Ununpentium" ),
    ( 116, "Lv",  "Livermorium" ),
    ( 117, "Uus", "Ununseptium" ),
    ( 118, "Uuo", "Ununoctium" ) )

def ZToSymbol( Z ) :
    """Returns the symbol for the specified Z or 'None' if Z is out-of-bounds."""

    for i in ZLabels:
        if ( i[0] == Z ) : return i[1]
    return None

def ZToLabel( Z ) :
    """Returns the label (i.e., name) for the specified Z or 'None' if Z is out-of-bounds."""
    
    for i in ZLabels:
        if ( i[0] == Z ) : return i[2]
    return None

def SymbolToZ( symbol ) :
    """Returns the Z for the specified symbol or 'None' if no match for symbol."""

    for i in ZLabels :
        if( i[1] == symbol ) : return i[0]
    return None

def LabelToZ( label ) :
    """Returns the Z for the specified label or 'None' if no match for label."""

    for i in ZLabels :
        if( i[2] == label ) : return i[0]
    return None

def gndNameToZ_A_Suffix( name ):
    """Returns the tuple (Z, A, suffix, ZA) for an gnd isotope name (e.g., gnd name = 'Am242_m1' 
    returns ( 95, 242, 'm1', 95242 ). Replaces the endl2.py function gndNameToEndlZ_A_Suffix."""

    if( name == 'n' ) : return( 0, 1, '', 1 )
    if( name == 'gamma' ) : return( 0, 0, '', 0 )
    if( name[:19] == 'FissionProduct_ENDL' ) :
        ZA = int( name[19:] )
        Z = ZA / 1000
        A = 1000 * Z - ZA
        return( Z, A, '', ZA )
    if( '__' in name ) : raise Exception ( "Name = %s" % name )
    realname = name.split( '__' )[0]
    naturalSuffix = ''
    if( '_' in realname ) :         # Isotope names can have level designator (e.g., 'O16_e3') and naturals are of the form 'S_natural' or 'S_natural_l'
        s = realname.split( '_' )   # where S is element's symbol and l is level designator (e.g., 'Xe_natural' or 'Xe_natural_c').
        sZA, suffix = s[:2]
        if( len( s ) > 2 ) :
            if( ( len( s ) > 3 ) or ( suffix != 'natural' ) ) : raise Exception( 'Invalid name for endl ZA particle = %s' % name )
            naturalSuffix = s[2]
    else :
        sZA = realname
        suffix = ''
    for i, c in enumerate( sZA ) :
        if( c.isdigit( ) ) : break
    if( not c.isdigit( ) ) : i += 1
    sZ, sA = sZA[:i], sZA[i:]
    Z = SymbolToZ( sZ )
    if( Z == None ) : raise Exception( 'No element symbol for particle named %s' % name )
    if( sA == '' ) :
        if( suffix == 'natural' ) : return( Z, 0, naturalSuffix, 1000 * Z )
        raise Exception( 'No A for particle named %s' % name )
    elif( suffix == 'natural' ) :
        raise Exception( 'Natural element also has A defined for particle named %s' % name )
    else :
        try :
            A = int( sA )
        except :
            raise Exception( 'Could not convert A to an integer for particle named %s' % name )
    ZA = 1000 * Z + A
    return( Z, A, suffix, ZA )

def ZAToGNDName( ZA ):
    """ use this instead of endl2.endlToGNDName"""

    if( ZA == 1 ) : return( 'n' )
    if type(ZA) in (int,float):
        Z,A = divmod(ZA,1000)
        sym = ZToSymbol(Z)
        if A==0: A = '_natural'
        return '%s%s' % (sym,A)
    elif type(ZA) is str:
        # should be of form 'zaZZZAAA_suffix'
        Z,A = divmod( int(ZA[2:8]),1000 )
        sym = ZToSymbol(Z)
        suffix = ZA[8:].strip()
        if A==0: A = '_natural'
        if suffix:
            if suffix=='m': suffix='m1'
            return '%s%s_%s' % (sym,A,suffix)
        return '%s%s' % (sym,A)
