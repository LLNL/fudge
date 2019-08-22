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

"""Angular/Energy double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

import xData.multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import xs_pdf_cdf as xs_pdf_cdfModule
from . import angular as angularModule

__metaclass__ = type

class xs_pdf_cdf1d( xs_pdf_cdfModule.xs_pdf_cdf1d ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( xs_pdf_cdf1d, ) )

class XYs3d( multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class subform( baseModule.subform ) :
    """Abstract base class for angularEnergyMC subforms."""

    def __init__( self, data ) :

        baseModule.subform.__init__( self )
        if( not( isinstance( data, self.allowedSubElements ) ) ) : raise TypeError( 'Invalid instance: %s' % type( data ) )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.data.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class angular( subform ) :

    moniker = 'angular'
    allowedSubElements = ( angularModule.XYs2d, )
    ancestryMembers = ( 'angular', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    angularModule.XYs2d.moniker       : angularModule.XYs2d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _angular= angular( angularSubform )

        xPath.pop( )
        return( _angular )

class angularEnergy( subform ) :

    moniker = 'angularEnergy'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'angularEnergy', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    XYs3d.moniker       : XYs3d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angularEnergy subform "%s"' % subformElement.tag )
        angularEnergySubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _angularEnergy = angularEnergy( angularEnergySubform )

        xPath.pop( )
        return( _angularEnergy )

class form( baseModule.form ) :

    moniker = 'angularEnergyMC'
    subformAttributes = ( 'angular', 'angularEnergy' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angular, _angularEnergy ) :

        if( not( isinstance( _angular, angular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angular ) )
        if( not( isinstance( _angularEnergy, angularEnergy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angularEnergy ) )
        baseModule.form.__init__( self, label, productFrame, ( _angular, _angularEnergy ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angular.convertUnits( unitMap )
        self.angularEnergy.convertUnits( unitMap )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <angularEnergyMC> element from xml."""

        xPath.append( element.tag )

        _angular = None
        _angularEnergy = None
        for child in element :
            if( child.tag == angular.moniker ) :
                _angular = angular.parseXMLNode( child, xPath, linkData )
            elif( child.tag == angularEnergy.moniker ) :
                _angularEnergy = angularEnergy.parseXMLNode( child, xPath, linkData )
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % subformElement.tag )

        angularEnergyMC = form( element.get( "label" ), element.get('productFrame'), 
                _angular, _angularEnergy )

        xPath.pop( )
        return( angularEnergyMC ) 
