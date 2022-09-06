#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import database as databaseModule

pops = databaseModule.read('Answers/database4.py.out')

for particle in pops : print( "%-12s %s" % ( particle.moniker, particle.key ) )
