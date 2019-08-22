# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

MaxA = 400
continuumID = 1000000
sumID = continuumID + 1

from .. import misc as miscModule

chemicalElementZSymbolNames = (
    (   1, "H",  "Hydrogen" ),      (   2, "He", "Helium" ),        (   3, "Li", "Lithium" ),
    (   4, "Be", "Beryllium" ),     (   5, "B",  "Boron" ),         (   6, "C",  "Carbon" ),
    (   7, "N",  "Nitrogen" ),      (   8, "O",  "Oxygen" ),        (   9, "F",  "Fluorine" ),
    (  10, "Ne", "Neon" ),          (  11, "Na", "Sodium" ),        (  12, "Mg", "Magnesium" ),
    (  13, "Al", "Aluminium" ),     (  14, "Si", "Silicon" ),       (  15, "P",  "Phosphorus" ),
    (  16, "S",  "Sulphur" ),       (  17, "Cl", "Chlorine" ),      (  18, "Ar", "Argon" ),
    (  19, "K",  "Potassium" ),     (  20, "Ca", "Calcium" ),       (  21, "Sc", "Scandium" ),
    (  22, "Ti", "Titanium" ),      (  23, "V",  "Vanadium" ),      (  24, "Cr", "Chromium" ),
    (  25, "Mn", "Manganese" ),     (  26, "Fe", "Iron" ),          (  27, "Co", "Cobalt" ),
    (  28, "Ni", "Nickel" ),        (  29, "Cu", "Copper" ),        (  30, "Zn", "Zinc" ),
    (  31, "Ga", "Gallium" ),       (  32, "Ge", "Germanium" ),     (  33, "As", "Arsenic" ),
    (  34, "Se", "Selenium" ),      (  35, "Br", "Bromine" ),       (  36, "Kr", "Krypton" ),
    (  37, "Rb", "Rubidium" ),      (  38, "Sr", "Strontium" ),     (  39, "Y",  "Yttrium" ),
    (  40, "Zr", "Zirconium" ),     (  41, "Nb", "Niobium" ),       (  42, "Mo", "Molybdenum" ),
    (  43, "Tc", "Technetium" ),    (  44, "Ru", "Ruthenium" ),     (  45, "Rh", "Rhodium" ),
    (  46, "Pd", "Palladium" ),     (  47, "Ag", "Silver" ),        (  48, "Cd", "Cadmium" ),
    (  49, "In", "Indium" ),        (  50, "Sn", "Tin" ),           (  51, "Sb", "Antimony" ),
    (  52, "Te", "Tellurium" ),     (  53, "I",  "Iodine" ),        (  54, "Xe", "Xenon" ),
    (  55, "Cs", "Cesium" ),        (  56, "Ba", "Barium" ),        (  57, "La", "Lanthanum" ),
    (  58, "Ce", "Cerium" ),        (  59, "Pr", "Praseodymium" ),  (  60, "Nd", "Neodymium" ),
    (  61, "Pm", "Promethium" ),    (  62, "Sm", "Samarium" ),      (  63, "Eu", "Europium" ),
    (  64, "Gd", "Gadolinium" ),    (  65, "Tb", "Terbium" ),       (  66, "Dy", "Dysprosium" ),
    (  67, "Ho", "Holmium" ),       (  68, "Er", "Erbium" ),        (  69, "Tm", "Thulium" ),
    (  70, "Yb", "Ytterbium" ),     (  71, "Lu", "Lutetium" ),      (  72, "Hf", "Hafnium" ),
    (  73, "Ta", "Tantalum" ),      (  74, "W",  "Tungsten" ),      (  75, "Re", "Rhenium" ),
    (  76, "Os", "Osmium" ),        (  77, "Ir", "Iridium" ),       (  78, "Pt", "Platinum" ),
    (  79, "Au", "Gold" ),          (  80, "Hg", "Mercury" ),       (  81, "Tl", "Thallium" ),
    (  82, "Pb", "Lead" ),          (  83, "Bi", "Bismuth" ),       (  84, "Po", "Polonium" ),
    (  85, "At", "Astatine" ),      (  86, "Rn", "Radon" ),         (  87, "Fr", "Francium" ),
    (  88, "Ra", "Radium" ),        (  89, "Ac", "Actinium" ),      (  90, "Th", "Thorium" ),
    (  91, "Pa", "Protactinium" ),  (  92, "U",  "Uranium" ),       (  93, "Np", "Neptunium" ),
    (  94, "Pu", "Plutonium" ),     (  95, "Am", "Americium" ),     (  96, "Cm", "Curium" ),
    (  97, "Bk", "Berkelium" ),     (  98, "Cf", "Californium" ),   (  99, "Es", "Einsteinium" ),
    ( 100, "Fm", "Fermium" ),       ( 101, "Md", "Mendelevium" ),   ( 102, "No", "Nobelium" ),
    ( 103, "Lr", "Lawrencium" ),    ( 104, "Rf", "Rutherfordium" ), ( 105, "Db", "Dubnium" ),
    ( 106, "Sg", "Seaborgium" ),    ( 107, "Bh", "Bohrium" ),       ( 108, "Hs", "Hassium" ),
    ( 109, "Mt", "Meitnerium" ),    ( 110, "Ds", "Darmstadtium" ),  ( 111, "Rg", "Roentgenium" ),
    ( 112, "Cn", "Copernicium" ),   ( 113, "Nh", "Nihonium" ),      ( 114, "Fl", "Flerovium" ),
    ( 115, "Mc", "Moscovium" ),     ( 116, "Lv", "Livermorium" ),   ( 117, "Ts", "Tennessine" ),
    ( 118, "Og", "Oganesson" ) )

symbolFromZ = {}
for Z, symbol, name in chemicalElementZSymbolNames : symbolFromZ[Z] = symbol

nameFromZ = {}
for Z, symbol, name in chemicalElementZSymbolNames : nameFromZ[Z] = name

ZFromSymbol = {}
for Z, symbol, name in chemicalElementZSymbolNames : ZFromSymbol[symbol] = Z

nameFromSymbol = {}
for Z, symbol, name in chemicalElementZSymbolNames : nameFromSymbol[symbol] = name

ZFromName = {}
for Z, symbol, name in chemicalElementZSymbolNames : ZFromName[name] = Z

symbolFromName = {}
for Z, symbol, name in chemicalElementZSymbolNames : symbolFromName[name] = symbol

def checkZ( Z ) :

    if( not( isinstance( Z, int ) ) ) : raise TypeError( 'Z not an int object.' )
    if( 0 < Z <= chemicalElementZSymbolNames[-1][0] ) : return( Z )
    raise ValueError( 'Z = "%s" out of range.' % Z )

def checkSymbol( symbol ) :

    if( not( isinstance( symbol, str ) ) ) : raise TypeError( 'symbol not a str object.' )
    if( symbol not in ZFromSymbol ) : raise ValueError( 'Invalid symbol: %s.' % miscModule.toLimitedString( symbol ) )

def checkA( A ) :

    if( not( isinstance( A, int ) ) ) : raise TypeError( 'A not an int' )
    if( not( 0 <= A < MaxA ) ) : raise ValueError( 'integer A out of range: A = %s' % A )
    return( A )

def checkIndex( index ) :

    if( not( isinstance( index, int ) ) ) : raise TypeError( 'index attribute not int object.' )
    if( index < 0 ) : ValueError( 'Negative level index = %s' % index )
    return( index )

def chemicalElementALevelIDsAndAnti( id, qualifierAllowed = False ) :

    baseID, anti, qualifier = miscModule.baseAntiQualifierFromID( id, qualifierAllowed = qualifierAllowed )

    chemicalElementSymbol = None
    A = None
    levelID = None
    isNucleus = False
    if( baseID.count( '_e' ) < 2 ) :        # Can have 0 or 1 '_e' (e.g., 'O16', 'O16_e3') if chemicalElement derived.
        isotopeID, sep, levelID = baseID.partition( '_e' )
        if( sep == '_e' ) :
            try :
                levelID = int( levelID )
                checkIndex( levelID )
            except :
                levelID = None
                isotopeID = ''
        else :
            levelID = None

        if( len( isotopeID ) > 0 ) :
            for i1, character in enumerate( isotopeID ) :
                if( character.isdigit( ) ) : break
            isotopeIDUpper = isotopeID[0].upper( ) + isotopeID[1:]
            if( character.isdigit( ) ) :
                chemicalElementSymbol = isotopeIDUpper[:i1]
                try :
                    A = int( isotopeIDUpper[i1:] )
                    checkSymbol( chemicalElementSymbol )
                    checkA( A )
                    isNucleus = isotopeIDUpper != isotopeID
                    if( levelID is None ) : levelID = 0
                except :
                    chemicalElementSymbol = None
            else :
                try :
                    checkSymbol( isotopeID )
                    chemicalElementSymbol = isotopeID
                except :
                    chemicalElementSymbol = None

    if( chemicalElementSymbol is None ) : A, levelID, isNucleus = None, None, None
    return( baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier )

def ZAInfo( particle ) :

    from .. import IDs as IDsModule
    from ..families import nuclide as nuclideModule
    from ..families import nucleus as nucleusModule
    from ..families import baryon as baryonModule
    from . import isotope as isotopeModule

    level = 0
    Z = 0
    A = 0
    if( isinstance( particle, ( isotopeModule.isotope, ) ) ) :
        Z = particle.Z
        A = particle.A
    elif( isinstance( particle, ( nuclideModule.particle, ) ) ) :
        Z = particle.Z
        A = particle.A
        level = particle.index
    elif( isinstance( particle, ( nucleusModule.particle, ) ) ) :
        Z = particle.Z
        A = particle.A
        level = particle.index
    elif( isinstance( particle, ( baryonModule.particle, ) ) ) :
        if( particle.id == IDsModule.neutron ) :
            A = 1
        if( particle.id == IDsModule.proton ) :
            Z = 1
            A = 1

    try :
        return( Z, A, 1000 * Z + A, level )
    except :
        raise

def ZA( particle ) :

    return( ZAInfo( particle )[2] )

def idFromZAndA( Z, A ) :

    from .. import IDs as IDsModule

    if( ( Z == 0 ) and ( A == 1 ) ) : return( IDsModule.neutron )
    return( isotopeSymbolFromChemicalElementIDAndA( symbolFromZ[Z], A ) )

def idFromZA( ZA ) :

    return( idFromZAndA( ZA / 1000, ZA % 1000 ) )

def nucleusIDFromZAndA( Z, A ) :

    nucleusID = idFromZAndA( Z, A )
    return( nucleusID[0].lower( ) + nucleusID[1:] )

def hasNucleus( particle, nucleusReturnsTrue = False ) :

    from . import isotope as isotopeModule
    from ..families import nuclide as nuclideModule
    from ..families import nucleus as nucleusModule

    if( isinstance( particle, ( isotopeModule.isotope, nuclideModule.particle ) ) ) : return( True )
    if( nucleusReturnsTrue and isinstance( particle, nucleusModule.particle ) ) : return( True )
    return( False )

def isotopeSymbolFromChemicalElementIDAndA( elementID, A ) :

    checkSymbol( elementID )
    checkA( A )
    return( "%s%s" % ( elementID, A ) )

def nuclideIDFromIsotopeSymbolAndIndex( isotopeSymbol, index ) :

    checkIndex( index )
    if( index == 0 ) : return( isotopeSymbol )
    return( "%s_e%s" % ( isotopeSymbol, index ) )

def nucleusIDFromNuclideID( nuclideID ) :

    return( nuclideID[0].lower( ) + nuclideID[1:] )

def nuclideIDFromNucleusID( nucleusID ) :

    return( nucleusID[0].upper( ) + nucleusID[1:] )
