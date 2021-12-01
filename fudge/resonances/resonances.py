# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Container for all resonance parameters (resolved, unresolved and/or scattering radius)
"""

from pqu import PQU as PQUModule
import xData.ancestry as ancestryModule

__metaclass__ = type

class resonances( ancestryModule.ancestry ) :
    """ 
    This is the top-level class for storing resonance parameters.
    For light targets it may contain only a scattering radius. For heavier targets it typically
    contains a resolved and/or unresolved section.
    
    resonances also has a boolean flag 'reconstructCrossSection'. If False, either the cross section
    has already been reconstructed, or the parameters are given for information only and no reconstruction
    should be performed.
    """

    moniker = 'resonances'
    ancestryMembers = ('scatteringRadius', 'resolved', 'unresolved')
    children = ancestryMembers                      # This should be deprecated as ancestryMembers suffices.

    def __init__(self, scatteringRadius_, resolved_=None, unresolved_=None):

        ancestryModule.ancestry.__init__( self )

        self.__scatteringRadius = scatteringRadius_
        if( self.__scatteringRadius is not None ) : self.__scatteringRadius.setAncestor( self )
        self.__resolved = resolved_
        self.__unresolved = unresolved_

        for child in (self.scatteringRadius, self.resolved, self.unresolved):
            if child is not None: child.setAncestor( self )
        
    def  __str__( self ) :
        """ string representation """
        return( self.toString( simpleString = False ) )

    @property
    def scatteringRadius( self ) :
        """Returns a reference to self's scatteringRadius."""

        return( self.__scatteringRadius )

    @property
    def resolved( self ) :
        """Returns a reference to self's resolved."""

        return( self.__resolved )

    @property
    def unresolved( self ) :
        """Returns a reference to self's unresolved."""

        return( self.__unresolved )

    def convertUnits( self, unitMap ):
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        for child in self.children:
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits( unitMap )

    def check( self, info ):
        from fudge import warning
        warnings = []
        for child in self.children:
            section = getattr(self,child)
            if section is not None:
                warningList = section.check(info)
                if warningList:
                    warnings.append( warning.context( section.moniker, warningList ) )
        return warnings
    
    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
    
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for child in self.children:
            section = getattr(self, child)
            if section is not None:
                xmlString += section.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def domain( self, unitTo = None, asPQU = False ):
        """ Return resonance region domain as a tuple of floats: (lowest edge, highest edge).

        options:
          unitTo: convert output to specified unit (given as a string).
          asPQU = True: return a tuple of PhysicalQuantityWithUncertainty instances instead of floats.
        """

        bounds = [ (self.scatteringRadius.domainMin, self.scatteringRadius.domainMax) ]
        if self.resolved:
            if self.resolved.multipleRegions:
                bounds += [(reg.domainMin, reg.domainMax) for reg in self.resolved.regions]
            else: bounds.append( (self.resolved.domainMin, self.resolved.domainMax) )
        if self.unresolved:
            bounds.append( (self.unresolved.domainMin, self.unresolved.domainMax) )

        for idx in range(len(bounds)-1):
            assert bounds[idx][1] == bounds[idx+1][0], "Resonance region boundaries don't match!"

        if( asPQU ):
            return (bounds[0][0], bounds[-1][1])
        elif unitTo:
            return (bounds[0][0].getValue(unitTo), bounds[-1][1].getValue(unitTo))
        else:
            return (bounds[0][0].value, bounds[-1][1].value)

    @property
    def reconstructCrossSection( self ):
        if self.resolved and not self.resolved.evaluated.useForSelfShieldingOnly: return True
        if self.unresolved and not self.unresolved.evaluated.useForSelfShieldingOnly: return True
        return False

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        from .resolved import resolved
        from .scatteringRadius import scatteringRadius
        from .unresolved import unresolved

        xPath.append( element.tag )

        scatRadius, RRR, URR = None,None,None
        for child in element:
            if child.tag== scatteringRadius.moniker:
                scatRadius = scatteringRadius.parseXMLNode(child, xPath, linkData)
            elif child.tag== resolved.moniker:
                RRR = resolved.parseXMLNode(child, xPath, linkData)
            elif child.tag== unresolved.moniker:
                URR = unresolved.parseXMLNode(child, xPath, linkData)
            else:
                raise Exception("unknown element '%s' encountered in resonances!" % child.tag)

        res = cls( scatteringRadius_ = scatRadius, resolved_=RRR, unresolved_=URR )
        xPath.pop()
        return res
    
    def toString( self, simpleString = False ) :
        """Returns a string representation of self. If simpleString is True, 
        the string contains only an overview without listing resonances"""
        s = 'Resonances:\n'
        if self.scatteringRadius:
            s += self.scatteringRadius.toString( simpleString = simpleString )
        if self.resolved:
            s += self.resolved.toString( simpleString = simpleString )
        if self.unresolved:
            s += self.unresolved.toString( simpleString = simpleString )
        return( s )
