# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon coherent and incoherent base function sub-nodes.
"""

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule

from fudge import abstractClasses as abstractClassesModule

class Form( abstractClassesModule.Form ) :

    def __init__( self, pid, label, productFrame, subForms, identicalParticles = False ) :

        abstractClassesModule.Form.__init__( self )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid must be a str instance' )
        self.__pid = pid

        self.label = label

        self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(productFrame)

        self.__identicalParticles = identicalParticles

        for i1, subform in enumerate( subForms ) :
            setattr( self, self.subformAttributes[i1], subform )
            if( subform is not None ) : subform.setAncestor( self )
        self.subforms = [ subform for subform in subForms ]

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'value must be a string' )
        self.__label = value

    @property
    def pid( self ) :

        return( self.__pid )

    @property
    def productFrame( self ) :

        return( self.__productFrame )

    @property
    def identicalParticles( self ) :

        return( self.__identicalParticles )

    def isThermalNeutronScatteringLaw( self ) :

        return( False )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        if( hasattr( self, 'subforms' ) ) :
            for subform in self.subforms :
                if( subform is None ) : continue           # FIXME, needed for photo-atomic coherentScattering data with no anomalousScattering subforms.
                subform.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr = ' label="%s"' % self.label
        attributeStr += ' pid="%s" productFrame="%s"' % ( self.__pid, self.__productFrame )
        if( self.identicalParticles ) : attributeStr += ' identicalParticles="true"'
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        for subform in self.subforms : 
            if( subform is not None ) : xmlString += subform.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class XYs1d( XYs1dModule.XYs1d ) :

    pass

class Regions1d( regionsModule.Regions1d ) :

    pass

def defaultAxes( factorLabel, energyUnit ) :

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( factorLabel, 0, "" )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

def check( self, info ) :

    return( [] )
