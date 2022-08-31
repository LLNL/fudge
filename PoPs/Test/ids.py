#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import misc as miscModule

parts = miscModule.baseAntiQualifierFromID( 'e-' )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'e-_anti' )
print( parts )

parts = miscModule.baseAntiQualifierFromID( 'O' )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O_anti' )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O{abc}', qualifierAllowed = True )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O_anti{abc}', qualifierAllowed = True )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O16' )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O16_anti' )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O16{abc}', qualifierAllowed = True )
print( parts )
parts = miscModule.baseAntiQualifierFromID( 'O16_anti{abc}', qualifierAllowed = True )
print( parts )
