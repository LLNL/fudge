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
This module contains the class for the suite chemicalElements.
"""

from .. import suite as suiteModule

from ..groups import misc as chemicalElementMiscModule

from . import isotope as isotopeModule
from . import chemicalElement as chemicalElementModule

from ..families import nuclide as nuclideModule
from ..families import nucleus as nucleusModule

class chemicalElements( suiteModule.sortedSuite ) :
    """
    This class stores a list of chemicalElement instances.
    """

    moniker = 'chemicalElements'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, [ chemicalElementModule.chemicalElement ], replace = True )

    def __contains__( self, key ) :

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( chemicalElementSymbol is None ) : return( False )
        if( suiteModule.sortedSuite.__contains__( self, chemicalElementSymbol ) ) :
            return( key in suiteModule.sortedSuite.__getitem__( self, chemicalElementSymbol ) )
        return( False )

    def __getitem__( self, key ) :
        """
        :param key: string or integer
        :return: chemicalElement, nuclide or nucleus
        """
        try:
            baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        except:
            chemicalElementSymbol = None

        if( chemicalElementSymbol is not None ) :
            chemicalElement = suiteModule.sortedSuite.__getitem__( self, chemicalElementSymbol )
            if( baseID == chemicalElementSymbol ) : return( chemicalElement )
            return( chemicalElement[key] )

        return suiteModule.sortedSuite.__getitem__( self, key )

    def add( self, item ) :

        if( isinstance( item, chemicalElementModule.chemicalElement ) ) :
            suiteModule.sortedSuite.add( self, item )
        elif( isinstance( item, ( isotopeModule.isotope, nuclideModule.particle, nucleusModule.particle ) ) ) :
            chemicalElementSymbol = item.chemicalElementSymbol
            if( self.hasSymbol( chemicalElementSymbol ) ) :
                suiteModule.sortedSuite.__getitem__( self, chemicalElementSymbol ).add( item )
            else :
                chemicalElement = chemicalElementModule.chemicalElement( chemicalElementSymbol, item.Z, chemicalElementMiscModule.nameFromZ[item.Z] )
                suiteModule.sortedSuite.add( self, chemicalElement )
                chemicalElement.add( item )
        else :
            raise TypeError( "Particle not a nuclide or nucleus object." )

    def convertUnits( self, unitMap ) :

        for chemicalElement in self : chemicalElement.convertUnits( unitMap )

    def getSymbol( self, key ) :

        return( suiteModule.sortedSuite.__getitem__( self, key ) )

    def hasSymbol( self, key ) :

        return( suiteModule.sortedSuite.__contains__( self, key ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )
        for child in element :
            self.add( chemicalElementModule.chemicalElement.parseXMLNodeAsClass( child, xPath, linkData ) )
        xPath.pop()

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )
        self = cls( )
        xPath.pop()
        self.parseXMLNode( element, xPath, linkData )

        return( self )
