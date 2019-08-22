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
This module contains a special type of Q-value, unique to fission reactions.
Fission releases energy several different ways (neutrons, gammas, etc.), and it's useful to subdivide the Q-value
into these different terms.
"""

# FIXME this is really a form of Q-value. Should it be defined in Q.py instead?

import xData.ancestry as ancestryModule
import xData.standards as standardsModule
import xData.axes as axesModule
import xData.series1d as series1dModule
import xData.XYs as XYsModule
import xData.gridded as griddedModule
import xData.uncertainties as uncertaintiesModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( 'energy_out', 0, energyUnit )
    axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
    return( axes )

class polynomial1d( series1dModule.polynomial1d ):

    def __init__( self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, sep = ' ',
            coefficientUncertainties = None ):

        series1dModule.polynomial1d.__init__( self, coefficients=coefficients, domainMin=domainMin, domainMax=domainMax,
            lowerIndex=lowerIndex, axes=axes, index=index, valueType=valueType, value=value, label=label, sep=sep )
        if coefficientUncertainties is not None:
            self.uncertainty = uncertaintiesModule.uncertainty( functional =
                    series1dModule.polynomial1d( coefficientUncertainties, domainMin, domainMax, axes=axes ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class XYs1d( XYsModule.XYs1d ) :

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        multiGroup = miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self )
        multiGroup.label = None
        return( multiGroup )

class gridded1d( griddedModule.gridded1d ) :

    pass

class fissionEnergyReleaseTerm( ancestryModule.ancestry ):
    """ Base class for all types of fission energy release. """

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        ancestryModule.ancestry.__init__( self )
        self.data = data

    @property
    def data(self):
        return self.__data

    @data.setter
    def data( self, value ) :

        if not isinstance( value, ( XYs1d, polynomial1d, gridded1d ) ) :
            raise TypeError( 'Invalid class "%s" for fissionEnergyReleaseTerm' % type(value) )
        value.setAncestor(self)
        self.__data = value

    def check(self, info):
        from fudge.gnds import warning

        warnings = []
        linearized = self.data.toPointwise_withLinearXYs(accuracy=1e-6)
        if linearized.rangeMin < 0:
            warnings.append( warning.badFissionEnergyRelease( worstCase=linearized.rangeMin, obj=self) )
        if not 0.1 < linearized.rangeMax / linearized.rangeMin < 10:
            warnings.append( warning.badFissionEnergyRelease( worstCase=linearized.rangeMax, obj=self) )
        return warnings

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.__class__( self.data.processMultiGroup( style, tempInfo, indent ) ) )

    def toXMLList( self, indent="", **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xmllist = ['%s<%s>' % (indent, self.moniker)]
        xmllist += self.data.toXMLList( indent2, **kwargs )
        xmllist[-1] += '</%s>' % self.moniker
        return xmllist

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        dataClass = {
            XYs1d.moniker: XYs1d,
            polynomial1d.moniker: polynomial1d,
            gridded1d.moniker: gridded1d
        }.get( element[0].tag, None )
        data = dataClass.parseXMLNode( element[0], xPath, linkData )
        FERT = cls( data )
        xPath.pop()
        return FERT

class promptProductKE( fissionEnergyReleaseTerm ):

    moniker = 'promptProductKE'

class promptNeutronKE( fissionEnergyReleaseTerm ):

    moniker = 'promptNeutronKE'

class delayedNeutronKE( fissionEnergyReleaseTerm ):

    moniker = 'delayedNeutronKE'

class promptGammaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'promptGammaEnergy'

class delayedGammaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'delayedGammaEnergy'

class delayedBetaEnergy( fissionEnergyReleaseTerm ):

    moniker = 'delayedBetaEnergy'

class neutrinoEnergy( fissionEnergyReleaseTerm ):

    moniker = 'neutrinoEnergy'

class nonNeutrinoEnergy( fissionEnergyReleaseTerm ):

    moniker = 'nonNeutrinoEnergy'

class totalEnergy( fissionEnergyReleaseTerm ):

    moniker = 'totalEnergy'

class field:
    """
    Descriptor to ensure ancestry is set when adding energy release terms
    to fissionEnergyReleased class.
    """

    def __init__(self, Class):
        self.Class = Class
        self.fname = '__' + Class.moniker

    def __get__(self, instance, owner):
        return getattr( instance, self.fname, None )

    def __set__(self, instance, value):
        if not isinstance( value, self.Class ):
            raise TypeError( "Incorrect type: expected %s, got %s" % (self.Class.moniker, type(value)) )
        value.setAncestor(instance)
        setattr( instance, self.fname, value )

class fissionEnergyReleased( ancestryModule.ancestry ) :
    """
    Store average energy released to different types of fission products.
    (prompt and delayed neutrons, prompt / delayed gammas, betas, neutrinos, etc.)
    Each term is currently (when translating from ENDF) stored as a polynomial expansion,
    although we expect to also see XYs1d representations in future evaluations
    """

    moniker = 'fissionEnergyReleased'

    ancestryMembers = ( 'promptProductKE',      'promptNeutronKE',      'delayedNeutronKE', 'promptGammaEnergy',
                        'delayedGammaEnergy',   'delayedBetaEnergy',    'neutrinoEnergy',   'nonNeutrinoEnergy',
                        'totalEnergy' )

    promptProductKE = field( promptProductKE )
    promptNeutronKE = field( promptNeutronKE )
    delayedNeutronKE = field( delayedNeutronKE )
    promptGammaEnergy = field( promptGammaEnergy )
    delayedGammaEnergy = field( delayedGammaEnergy )
    delayedBetaEnergy = field( delayedBetaEnergy )
    neutrinoEnergy = field( neutrinoEnergy )
    nonNeutrinoEnergy = field( nonNeutrinoEnergy )
    totalEnergy = field( totalEnergy )

    def __init__( self, label, **kwargs ) :

        ancestryModule.ancestry.__init__( self )
        self.label = label

        for key in kwargs.keys() :
            if key in self.ancestryMembers:
                setattr( self, key, kwargs.pop( key, None ) )
        if kwargs:
            raise TypeError("fissionEnergyReleased received unexpected argument(s): '%s'" % kwargs.keys())

    def __iter__( self ) :

        for term in self.ancestryMembers : yield getattr( self, term )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        for term in self.ancestryMembers :
            fwarnings = getattr( self, term ).convertUnits( unitMap )

    def check( self, info ) :

        from fudge.gnds import warning
        warnings = []
        for term in self.ancestryMembers:
            fwarnings = getattr(self, term).check( info )
            if fwarnings:
                warnings.append( warning.context('%s:' % term, fwarnings) )
        return warnings

    def processMultiGroup( self, style, tempInfo, indent ) :

        kwargs = {}
        for term in self.ancestryMembers :
            kwargs[term] = getattr( self, term ).processMultiGroup( style, tempInfo, indent )
        return( fissionEnergyReleased( style.label, **kwargs ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.promptProductKE.data.toPointwise_withLinearXYs( **kwargs ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent="", **kwargs ):

        indent2 = indent+"  "
        xmlList = ['%s<%s label="%s">' % ( indent, self.moniker, self.label )]
        for term in self.ancestryMembers:
            xmlList += getattr(self, term).toXMLList( indent2, **kwargs )
        xmlList[-1] += "</%s>" % self.moniker

        return xmlList

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """Parse <fissionEnergyReleased> from xml."""

        xPath.append( element.tag )
        children = {}
        childClasses = dict( [(Class.moniker, Class) for Class in
            (   promptProductKE, promptNeutronKE, delayedNeutronKE, promptGammaEnergy,
                delayedGammaEnergy, delayedBetaEnergy, neutrinoEnergy, nonNeutrinoEnergy,
                totalEnergy ) ] )
        for child in element:
            children[ child.tag ] = childClasses[ child.tag ].parseXMLNode( child, xPath, linkData )
        fer = fissionEnergyReleased( label = element.get( "label" ), **children )
        xPath.pop()
        return fer
