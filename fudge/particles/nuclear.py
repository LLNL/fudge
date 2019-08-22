# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

elementsZSymbolName = {
      0 : ( "n",  "Neutron" ),        1 : ( "H",  "Hydrogen" ),       2 : ( "He", "Helium" ),         3 : ( "Li", "Lithium" ),       4 : ( "Be", "Beryllium" ), 
      5 : ( "B",  "Boron" ),          6 : ( "C",  "Carbon" ),         7 : ( "N",  "Nitrogen" ),       8 : ( "O",  "Oxygen" ),        9 : ( "F",  "Fluorine" ), 
     10 : ( "Ne", "Neon" ),          11 : ( "Na", "Sodium" ),        12 : ( "Mg", "Magnesium" ),     13 : ( "Al", "Aluminium" ),    14 : ( "Si", "Silicon" ), 
     15 : ( "P",  "Phosphorus" ),    16 : ( "S",  "Sulphur" ),       17 : ( "Cl", "Chlorine" ),      18 : ( "Ar", "Argon" ),        19 : ( "K",  "Potassium" ),
     20 : ( "Ca", "Calcium" ),       21 : ( "Sc", "Scandium" ),      22 : ( "Ti", "Titanium" ),      23 : ( "V",  "Vanadium" ),     24 : ( "Cr", "Chromium" ), 
     25 : ( "Mn", "Manganese" ),     26 : ( "Fe", "Iron" ),          27 : ( "Co", "Cobalt" ),        28 : ( "Ni", "Nickel" ),       29 : ( "Cu", "Copper" ), 
     30 : ( "Zn", "Zinc" ),          31 : ( "Ga", "Gallium" ),       32 : ( "Ge", "Germanium" ),     33 : ( "As", "Arsenic" ),      34 : ( "Se", "Selenium" ), 
     35 : ( "Br", "Bromine" ),       36 : ( "Kr", "Krypton" ),       37 : ( "Rb", "Rubidium" ),      38 : ( "Sr", "Strontium" ),    39 : ( "Y",  "Yttrium" ),
     40 : ( "Zr", "Zirconium" ),     41 : ( "Nb", "Niobium" ),       42 : ( "Mo", "Molybdenum" ),    43 : ( "Tc", "Technetium" ),   44 : ( "Ru", "Ruthenium" ), 
     45 : ( "Rh", "Rhodium" ),       46 : ( "Pd", "Palladium" ),     47 : ( "Ag", "Silver" ),        48 : ( "Cd", "Cadmium" ),      49 : ( "In", "Indium" ), 
     50 : ( "Sn", "Tin" ),           51 : ( "Sb", "Antimony" ),      52 : ( "Te", "Tellurium" ),     53 : ( "I",  "Iodine" ),       54 : ( "Xe", "Xenon" ), 
     55 : ( "Cs", "Cesium" ),        56 : ( "Ba", "Barium" ),        57 : ( "La", "Lanthanum" ),     58 : ( "Ce", "Cerium" ),       59 : ( "Pr", "Praseodymium" ),
     60 : ( "Nd", "Neodymium" ),     61 : ( "Pm", "Promethium" ),    62 : ( "Sm", "Samarium" ),      63 : ( "Eu", "Europium" ),     64 : ( "Gd", "Gadolinium" ), 
     65 : ( "Tb", "Terbium" ),       66 : ( "Dy", "Dysprosium" ),    67 : ( "Ho", "Holmium" ),       68 : ( "Er", "Erbium" ),       69 : ( "Tm", "Thulium" ), 
     70 : ( "Yb", "Ytterbium" ),     71 : ( "Lu", "Lutetium" ),      72 : ( "Hf", "Hafnium" ),       73 : ( "Ta", "Tantalum" ),     74 : ( "W",  "Tungsten" ), 
     75 : ( "Re", "Rhenium" ),       76 : ( "Os", "Osmium" ),        77 : ( "Ir", "Iridium" ),       78 : ( "Pt", "Platinum" ),     79 : ( "Au", "Gold" ),
     80 : ( "Hg", "Mercury" ),       81 : ( "Tl", "Thallium" ),      82 : ( "Pb", "Lead" ),          83 : ( "Bi", "Bismuth" ),      84 : ( "Po", "Polonium" ), 
     85 : ( "At", "Astatine" ),      86 : ( "Rn", "Radon" ),         87 : ( "Fr", "Francium" ),      88 : ( "Ra", "Radium" ),       89 : ( "Ac", "Actinium" ), 
     90 : ( "Th", "Thorium" ),       91 : ( "Pa", "Protactinium" ),  92 : ( "U",  "Uranium" ),       93 : ( "Np", "Neptunium" ),    94 : ( "Pu", "Plutonium" ), 
     95 : ( "Am", "Americium" ),     96 : ( "Cm", "Curium" ),        97 : ( "Bk", "Berkelium" ),     98 : ( "Cf", "Californium" ),  99 : ( "Es", "Einsteinium" ),
    100 : ( "Fm", "Fermium" ),      101 : ( "Md", "Mendelevium" ),  102 : ( "No", "Nobelium" ),     103 : ( "Lr", "Lawrencium" ),  104 : ( "Rf", "Rutherfordium" ), 
    105 : ( "Db", "Dubnium" ),      106 : ( "Sg", "Seaborgium" ),   107 : ( "Bh", "Bohrium" ),      108 : ( "Hs", "Hassium" ),     109 : ( "Mt", "Meitnerium" ), 
    110 : ( "Ds", "Darmstadtium" ), 111 : ( "Rg", "Roentgenium" ),  112 : ( "Cn", "Copernicium" ),  113 : ( "Uut", "Ununtrium" ),  114 : ( "Uuq", "Ununquadium" ), 
    115 : ( "Uup", "Ununpentium" ), 116 : ( "Uuh", "Ununhexium" ),  117 : ( "Uus", "Ununseptium" ), 118 : ( "Uuo", "Ununoctium" ) }

elementsSymbolZ = {}
for Z in elementsZSymbolName : elementsSymbolZ[elementsZSymbolName[Z][0]] = Z

def elementSymbolFromZ( Z ) :

    Z = int( Z )
    return( elementsZSymbolName[Z][0] )

def elementNameFromZ( Z ) :

    Z = int( Z )
    return( elementsZSymbolName[Z][1] )

def nucleusNameFromZAndA( Z, A ) :

    A = int( A )
    if( A < 1 ) : raise Exception( 'Invalid A = %s for Z = %s' % ( A, Z ) )
    AStr = str( A )

    Z = int( Z )
    if( Z == 0 ) :
        if( A != 1 ) : raise Exception( 'Invalid A = %s for neutron (Z = 0)' % A )
        AStr = ''

    return( elementSymbolFromZ( Z ) + AStr )

def nucleusNameFromZA( ZA ) :

    return( nucleusNameFromZAndA( ZA / 1000, ZA % 1000 ) )

def getZAOrNameAs_xParticle( particle ) :

    from fudge.gnd import xParticle
    particle_ = particle
    if( type( particle_ ) == type( 1 ) ) : particle_ = nucleusNameFromZA( particle_ )
    if( type( particle_ ) != type( "" ) ) : raise Exception( "Invalid particle = %s type" % particle )
    return( xParticle.xParticle( particle_, xParticle.particleType_Nuclear, getMassFromName( particle_ ) ) )

def getMassFromName( name ) :

    from fudge.structure import masses
    if( name == 'n' ) : return( masses.getMassWithUnitFromZA( 1 ) )
    if( '_natural' in name ) : 
        symbol, A = name.split( '_' )[0], 0
    else :
        for i, l in enumerate( name ) :
            if( l.isdigit( ) ) : break
        symbol, AStr = name[:i], name[i:]
        for i, l in enumerate( AStr ) :
            if( not( l.isdigit( ) ) ) : break
        A = int( AStr[:i+1] )
    Z = elementsSymbolZ[symbol]
    return( masses.getMassWithUnitFromZA( 1000 * Z + A ) )

def getZandAFromName( name ) :

    if( name == 'n' ) : return( ( 0, 1 ) )
    for i, l in enumerate( name ) :
        if( l.isdigit( ) ) : break
    symbol, AStr = name[:i], name[i:]
    Z = elementsSymbolZ[symbol]
    for i, l in enumerate( AStr ) :
        if( not( l.isdigit( ) ) ) : break
    A = int( AStr[:i+1] )
    return (Z,A)
