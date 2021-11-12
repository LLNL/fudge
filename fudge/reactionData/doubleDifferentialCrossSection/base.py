# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon coherent and incoherent base function sub-nodes.
"""

import xData.standards as standardsModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule

from fudge import abstractClasses as abstractClassesModule

class form( abstractClassesModule.form ) :

    def __init__( self, pid, label, productFrame, subForms, identicalParticles = False ) :

        abstractClassesModule.form.__init__( self )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid must be a str instance' )
        self.__pid = pid

        self.label = label

        if( productFrame not in standardsModule.frames.allowedFrames ) :
            raise TypeError( 'Invalid productFrame = "%s"' % productFrame )
        self.__productFrame = productFrame

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

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr = ' label="%s"' % self.label
        attributeStr += ' pid="%s" productFrame="%s"' % ( self.__pid, self.__productFrame )
        if( self.identicalParticles ) : attributeStr += ' identicalParticles="true"'
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        for subform in self.subforms : 
            if( subform is not None ) : xmlString += subform.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class XYs1d( XYsModule.XYs1d ) :

    pass

class regions1d( regionsModule.regions1d ) :

    pass

def defaultAxes( factorLabel, energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( factorLabel, 0, "" )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

def check( self, info ) :

    return( [] )
