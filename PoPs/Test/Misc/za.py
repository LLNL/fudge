# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xml.etree import cElementTree

from PoPs import database as databaseModule
from PoPs import alias as aliasModule
from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nucleus as nucleusModule
from PoPs.groups import chemicalElement as chemicalElementModule

from PoPs import misc as miscModule

pops = databaseModule.database.readFile( 'p.xml' )

ids_ZAs = []
for particle in pops :
    ids_ZAs.append( [ miscModule.ZA( particle ), particle.id ] )

ids_ZAs.sort( )
for ZA, _id in ids_ZAs : print("%-8s %s" % (ZA, _id))

ids_ZAs = []
for particle in pops :
    ids_ZAs.append( [ miscModule.ZAInfo( particle ), particle.id ] )

ids_ZAs.sort( )
for ZA, _id in ids_ZAs : print("%-8s %s" % (ZA, _id))
