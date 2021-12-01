# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU
from xData import ancestry as ancestryModule, XYs as XYsModule, constant as constantModule

__metaclass__ = type

class scatteringRadius( ancestryModule.ancestry ):
    """
    The scatteringRadius determines scattering length. May be constant or energy-dependent.
    Each resonances section contains a scatteringRadius. The resolved / unresolved regions may also
    provide their own copy that overrides the top-level definition.
    The RMatrix class also provides a way of defining channel-specific radii (deprecated but necessary to handle
    ENDF evaluations)
    """

    moniker = 'scatteringRadius'

    def __init__( self, form=None ) :

        ancestryModule.ancestry.__init__( self )
        form.setAncestor( self )
        self.form = form

    def __eq__(self, other):
        if not isinstance(other, scatteringRadius): return False
        return( self.form==other.form )

    def __str__(self): return self.form.moniker

    def __bool__( self ) :

        return bool( self.form )

    __nonzero__ = __bool__

    def toString( self, simpleString = False ) :
        """Returns a string representation of self. If simpleString is True,
        the string contains only an overview without listing resonances"""
        if simpleString: return str(self)
        return str(self)

    def copy( self ):

        return self.__class__( self.form.copy() )

    def check( self, info ):
        """
        Checks that the scattering radius is within a factor of 3 of the expected scattering radius
        of AP = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08

        This default AP is roughly the nuclear radius, so getting within a factor of 3 of the default
        shouldn't be a problem.

        A similar test is included in PSYCHE, but PSYCHE's cannot handle LRF=7
        :param info:
        :return:
        """
        from fudge import warning
        warnings = []
        target = info['reactionSuite'].PoPs[ info['reactionSuite'].target ]
        if target.id in info['reactionSuite'].PoPs.aliases:
            target = info['reactionSuite'].PoPs[ target.pid ]
        expectedAP = 10.0*( 0.123 * target.getMass('amu')**(1./3.) + 0.08 ) # expected radius in fm
        factor=3.0
        if self.isEnergyDependent():
            egrid=self.form.domainGrid
            APs = self.getValueAs( 'fm', energy_grid=egrid )
            for iE,AP in enumerate(APs):
                if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                    warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP, E=egrid[iE])
        else:
            AP = self.form.value
            if self.form.rangeUnit != 'fm':
                AP *= PQU.PQU(1, self.form.rangeUnit).getValueAs('fm')
            if AP/expectedAP > factor or AP/expectedAP < 1./factor:
                warning.badScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP)
        return warnings

    def convertUnits( self, unitMap ):

        self.form.convertUnits( unitMap )

    def isEnergyDependent(self):
        return isinstance( self.form, XYsModule.XYs1d )

    def getValueAs( self, unit, energy_grid=None ):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise NameError("Missing: energy_grid to evaluate E-dependent scattering radius")
            energy_unit = self.form.axes[-1].unit
            xScale = self.form.domainUnitConversionFactor( energy_unit )
            yScale = self.form.rangeUnitConversionFactor( unit )
            return [ yScale * self.form.evaluate( xScale * e ) for e in energy_grid ]
        else:
            oldUnit = self.form.rangeUnit
            factor = PQU.PQU(1, oldUnit).getValueAs( unit )
            return self.form.value * factor

    def toXMLList( self, indent = '', **kwargs ):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.form.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append( element.tag )
        form = {
            constantModule.constant1d.moniker: constantModule.constant1d,
            XYsModule.XYs1d.moniker: XYsModule.XYs1d,
        }[ element[0].tag ].parseXMLNode( element[0], xPath, linkData )
        SR = cls( form )
        xPath.pop()
        return SR

class hardSphereRadius(scatteringRadius):

    moniker = 'hardSphereRadius'
