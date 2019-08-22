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

"""
This module contains the nuclear level classes.
"""

continuumID = 1000000
sumID = continuumID + 1

from .. import misc as miscModule

from ..quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from . import particle as particleModule
from . import nucleus as nucleusModule

class alias( particleModule.alias ) :

    moniker = 'nuclearLevelAlias'

    @property
    def chemicalElement( self ) :

        return( self.__particle.chemicalElement )

    @property
    def Z( self ) :

        return( self.__particle.Z )

    @property
    def A( self ) :

        return( self.__particle.A )

    @property
    def index( self ) :

        return( self.__particle.index )

    @property
    def intIndex( self ) :

        return( self.__particle.intIndex )

    @property
    def energy( self ) :

        return( self.__particle.energy )

class particle( particleModule.particle ) :

    moniker = 'nuclearLevel'
    alias = alias

    def __init__( self, id, nucleus = None ) :

        particleModule.particle.__init__( self, id )

        self.nucleus = nucleus

    def __eq__( self, other ) :

        from ..groups import isotope as isotopeModule

        if(   isinstance( other, particle ) ) :
            return( self.id == other.id )
        elif( isinstance( other, isotopeModule.suite ) ) :
            _particle = other.particle
            return( self.id == _particle.id )
        else :
            return( False )

    @property
    def nucleus( self ) :

        return( self.__nucleus )

    @nucleus.setter
    def nucleus( self, _nucleus ) :

        if( not( isinstance( _nucleus, nucleusModule.particle ) ) ) : TypeError( 'Invalid nucleus' )
        self.__nucleus = _nucleus
        if( _nucleus is not None ) :
            nucleusID = self.id[0].lower( ) + self.id[1:]
            if( nucleusID != _nucleus.id ) :
                ValueError( '''nuclearLevel id = "%s" does not agree with the nucleus id = "%s"''' % ( id, _nucleus.id ) )
            _nucleus.setAncestor( self )

    @property
    def A( self ) :

        return( self.nucleus.__A )

    @property
    def intA( self ) :

        return( self.nucleus.__intA )

    @property
    def chemicalElement( self ) :

        return( self.nucleus.__chemicalElement )

    @property
    def index( self ) :

        return( self.nucleus.__index )

    @property
    def intIndex( self ) :

        return( self.nucleus.__intIndex )

    @property
    def energy( self ) :

        return( self.nucleus.energy )

    @property
    def Z( self ) :

        return( self.nucleus.__Z )

    def convertUnits( self, unitMap ) :

        particleModule.particle.convertUnits( self, unitMap )
        self.__nucleus.convertUnits( unitMap )

    def copy( self ) :

        _particle = particle( self.id )
        self.__copyStandardQuantities( _particle )
        _particle.nucleus = self.nucleus.copy( )
        return( _particle )

    def extraXMLElements( self, indent, **kwargs ) :

        if( self.nucleus is None ) : return( [] )
        return( self.nucleus.toXMLList( indent, **kwargs ) )

    def getMass( self, unit ) :

        if( len( self.mass ) > 0 ) : return( self.mass[0].float( unit ) )
        if( self.intIndex == 0 ) : raise Exception( 'Recursion detected as group-state does not have a mass: ID = %s.' % self.id )
        return( self.getAncestor( )[0].mass[0].float( unit ) + self.energy[0].float( unit + ' * c**2' ) )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        if( element.tag == nucleusModule.particle.moniker ) :
            self.nucleus = nucleusModule.particle.parseXMLNodeAsClass( element, xPath, linkData )
            return( True )

        return( False )

    def sortCompare( self, other ) :

        if( not( isinstance( other, particle ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.intIndex - other.intIndex )
