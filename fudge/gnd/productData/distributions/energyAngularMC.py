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

"""Energy/angular double differential distribution classes for Monte Carlo with xs, pdf, and cdf 1d distribtions."""

import xData.multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import xs_pdf_cdf as xs_pdf_cdfModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule

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
    """Abstract base class for energyAngular subforms."""

    def __init__( self, data ) :

        baseModule.subform.__init__( self )
        if( not( isinstance( data, self.allowedSubElements ) ) ) : raise TypeError( 'Invalid instance: %s' % type( data ) )
        self.data = data
        self.data.setAncestor( self )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.data.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class energy( subform ) :

    moniker = 'energy'
    allowedSubElements = ( energyModule.XYs2d, )
    ancestryMembers = ( 'energy', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    energyModule.XYs2d.moniker       : energyModule.XYs2d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown energy subform "%s"' % subformElement.tag )
        energySubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _energy = energy ( energySubform )

        xPath.pop( )
        return( _energy )

class energyAngular( subform ) :

    moniker = 'energyAngular'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'energyAngular', )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {    XYs3d.moniker       : XYs3d,
                        }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown energyAngular subform "%s"' % subformElement.tag )
        energyAngularSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )

        _energyAngular = energyAngular( energyAngularSubform )

        xPath.pop( )
        return( _energyAngular )

class form( baseModule.form ) :

    moniker = 'energyAngularMC'
    subformAttributes = ( 'energy', 'energyAngular' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _energy, _energyAngular ) :

        if( not( isinstance( _energy, energy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energy ) )
        if( not( isinstance( _energyAngular, energyAngular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energyAngular ) )
        baseModule.form.__init__( self, label, productFrame, ( _energy, _energyAngular ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <energyAngularMC> element from xml."""

        xPath.append( element.tag )

        _energy = None
        _energyAngular = None
        for child in element :
            if( child.tag == energy.moniker ) :
                _energy = energy.parseXMLNode( child, xPath, linkData )
            elif( child.tag == energyAngular.moniker ) :
                _energyAngular = energyAngular.parseXMLNode( child, xPath, linkData )
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % subformElement.tag )

        energyAngularMC = form( element.get( "label" ), element.get( "productFrame" ), _energy, _energyAngular )

        xPath.pop( )
        return( energyAngularMC )
