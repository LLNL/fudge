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

__metaclass__ = type

import math

from . import documentation as documentationModule, suites as suitesModule
from xData import XYs as XYsModule, physicalQuantity as physicalQuantityModule, gridded as griddedModule
from xData import ancestry as ancestryModule


class thermalScattering( ancestryModule.ancestry ):

    moniker = 'thermalScattering'
    ancestryMembers = ('mass', 'cutoffEnergy', 'coherentElastic', 'incoherentElastic', 'incoherentInelastic')

    def __init__(self, cutoffEnergy, material=None, MAT=None, mass=None, documentation=None, coherentElastic=None,
            incoherentElastic=None, incoherentInelastic=None):
        ancestryModule.ancestry.__init__( self )
        self.material = material
# BRB6 hardwired
        self.projectile='n'
        self.MAT = MAT
        self.mass = mass
        self.cutoffEnergy = cutoffEnergy
        self.documentation = documentation
# BRB6 need to check for correct class.
        self.coherentElastic = coherentElastic
        self.incoherentElastic = incoherentElastic
        self.incoherentInelastic = incoherentInelastic

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [ indent+'<%s material="%s" MAT="%s">' % (self.moniker,self.material,self.MAT) ]
        if self.documentation is not None: xml += self.documentation.toXMLList( indent2, **kwargs )
        xml += self.cutoffEnergy.toXMLList( indent2, **kwargs )
        xml += self.mass.toXMLList( indent2, **kwargs )
        for ts in (coherentElastic.moniker,incoherentElastic.moniker,incoherentInelastic.moniker):
            if getattr(self,ts) is not None:
                xml += getattr(self,ts).toXMLList( indent2, **kwargs )
        xml.append(indent+'</%s>' % self.moniker)
        return xml

    def check( self, **kwargs ) :
        """
        Check all data in the reactionSuite, returning a list of warnings.
        """

        from fudge.gnds import warning

        options = { }
        for key in kwargs:
            if key in options: options[key] = kwargs[key]
            else: raise KeyError, "check() received unknown keyword argument '%s'" % key

        # assemble some useful info, to be handed down to children's check() functions:
        info = { 'reactionSuite': self, }
        info.update( options )

        warnings = []


        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.material), warnings)
        result.info = info
        return result

    def saveToFile( self, filename, **kwargs ):

        with open(filename,'w') as fout:
# BRB6 hardwired
            fout.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            fout.write( '\n'.join( self.toXMLList(**kwargs) ) )

    @staticmethod
    def parseXMLNode( element ):

        xPath = [ thermalScattering.moniker ]
        linkData = {}
        try:
            kwargs = {
                'material': element.get('material'),
                'MAT': int( element.get('MAT') )
            }
            for child in element:
                _class = {
                    documentationModule.documentation.moniker: documentationModule.documentation,
                    mass.moniker: mass,
                    cutoffEnergy.moniker: cutoffEnergy,
                    coherentElastic.moniker: coherentElastic,
                    incoherentElastic.moniker: incoherentElastic,
                    incoherentInelastic.moniker: incoherentInelastic
                }.get( child.tag, None )
                if _class is None:
                    print("Warning: encountered unexpected element '%s' in thermalScattering!" % child.tag)
                kwargs[ child.tag ] = _class.parseXMLNode( child, xPath, linkData )
            TS = thermalScattering( **kwargs )
        except:
            print("Error encountered at xpath = /%s" % '/'.join(xPath))
            raise

        return TS


class mass( physicalQuantityModule.physicalQuantity ):
    moniker = 'mass'


class cutoffEnergy( physicalQuantityModule.physicalQuantity ):
    moniker = 'cutoffEnergy'


class coherentElastic( ancestryModule.ancestry ):

    moniker = 'coherentElastic'

    def __init__(self, S_table):
        ancestryModule.ancestry.__init__(self)
        self.S_table = S_table

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s>' % self.moniker]
        xml += self.S_table.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        Stab = S_table.parseXMLNode( element[0], xPath, linkData )
        cohEl = coherentElastic( Stab )
        xPath.pop()
        return cohEl


class S_table( ancestryModule.ancestry ):
    """ For elastic coherent, cumulative structure factor 'S' is given as function of incident energy and temp. """ 
    moniker = 'S_table'

    def __init__(self, gridded2d):
        self.gridded2d = gridded2d

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s>' % self.moniker]
        xml += self.gridded2d.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        g2d = griddedModule.gridded2d.parseXMLNode( element[0], xPath, linkData )
        Stab = S_table( g2d )
        xPath.pop()
        return Stab


class incoherentElastic( ancestryModule.ancestry ):
    moniker = 'incoherentElastic'

    def __init__(self, characteristicCrossSection, DebyeWaller):
        ancestryModule.ancestry.__init__(self)
        self.characteristicCrossSection = characteristicCrossSection
        self.DebyeWaller = DebyeWaller

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s>' % self.moniker]
        xml += self.characteristicCrossSection.toXMLList( indent2, **kwargs )
        xml += self.DebyeWaller.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        xsc = characteristicCrossSection.parseXMLNode(element.find(characteristicCrossSection.moniker), xPath, linkData)
        DW = DebyeWaller.parseXMLNode(element.find(DebyeWaller.moniker), xPath, linkData)
        incElastic = incoherentElastic( xsc, DW )

        xPath.pop()
        return incElastic

class characteristicCrossSection( physicalQuantityModule.physicalQuantity ):

    moniker = 'characteristicCrossSection'


class DebyeWaller( XYsModule.XYs1d ):
    """For incoherent elastic sections, all we need is a characteristic cross section and a
    temperature-dependent list of the Debye-Waller integral """

    mutableYUnit = False
    moniker = "DebyeWaller"

    def __init__(self, *args, **kwargs):
        XYsModule.XYs1d.__init__(self, *args, **kwargs)


class incoherentInelastic( ancestryModule.ancestry ):
    moniker = 'incoherentInelastic'

    def __init__(self, S_alpha_beta, calculatedAtThermal = False, asymmetric = False, atoms = None):
        ancestryModule.ancestry.__init__(self)
        self.scatteringAtoms = scatteringAtoms()
        if atoms is not None:
            for atom in atoms:
                self.scatteringAtoms.add( atom )
        self.S_alpha_beta = S_alpha_beta
        self.calculatedAtThermal = calculatedAtThermal
        self.asymmetric = asymmetric

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        attrs = ''
        for attr in ('calculatedAtThermal', 'asymmetric'):
            if getattr(self,attr): attrs += ' %s="true"' % attr
        xml = [indent+'<%s%s>' % (self.moniker, attrs)]
        xml += self.scatteringAtoms.toXMLList( indent2 )
        xml += self.S_alpha_beta.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( incoherentElastic.moniker )
        kwargs = {
            'calculatedAtThermal': element.get('calculatedAtThermal') == 'true',
            'asymmetric': element.get('asymmetric') == 'true'
        }
        atoms = [scatteringAtom.parseXMLNode( child, xPath, linkData)
                           for child in element.find('scatteringAtoms')]
        Sab = S_alpha_beta.parseXMLNode( element.find( S_alpha_beta.moniker ), xPath, linkData )
        incInelastic = incoherentInelastic( Sab, atoms = atoms, **kwargs )
        xPath.pop()
        return incInelastic


class scatteringAtoms( suitesModule.suite ):
    moniker = 'scatteringAtoms'

    def __init__( self ):
        suitesModule.suite.__init__( self, [scatteringAtom] )

class scatteringAtom( ancestryModule.ancestry ):
    """ Inelastic incoherent scattering requires a description of each atom that is scattered off """
    moniker = 'scatteringAtom'

    def __init__(self, label, numberPerMolecule, mass, freeAtomCrossSection, e_critical = None, e_max = None,
            functionalForm = None, T_effective = None ):
        ancestryModule.ancestry.__init__(self)
        self.label = label
        self.numberPerMolecule = numberPerMolecule
        self.mass = mass
        self.freeAtomCrossSection = freeAtomCrossSection
        self.e_critical = e_critical
        self.e_max = e_max
        self.functionalForm = functionalForm
        self.T_effective = T_effective

    def boundCrossSection(self):
        return self.freeAtomCrossSection.value * ((self.mass.value + 1) / self.mass.value)**2

    def shortCollisionTime( self, alpha, beta, T ):
        """
        Equation 7.8 in ENDF manual
        :param alpha:
        :param beta:
        :param T:
        :return:
        """

        Teff = self.T_effective.evaluate(T)
        #num = math.exp( -(alpha-abs(beta))**2 * T / (4*alpha*Teff) - (beta + abs(beta))/2)
        # reformulate using T/Teff = 1-f, should be more numerically stable
        f = 1 - T/Teff
        num = math.exp( (f * (alpha - abs(beta))**2 - (alpha + beta)**2) / (4*alpha) )
        den = math.sqrt( 4 * math.pi * alpha * Teff / T )
        return num/den

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ''
        if self.functionalForm is not None: attrs = ' functionalForm="%s"' % self.functionalForm
        xml = [indent+'<%s label="%s" numberPerMolecule="%i"%s>' % (self.moniker, self.label, self.numberPerMolecule, attrs)]
        for attribute in ('mass','freeAtomCrossSection','e_critical','e_max','T_effective'):
            if getattr(self, attribute) is not None:
                xml += getattr(self,attribute).toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        kwargs = {
            'label': element.get('label'),
            'numberPerMolecule': int(element.get('numberPerMolecule')),
            'functionalForm': element.get('functionalForm'),
        }
        for child in element:
            _class = {
                mass.moniker: mass,
                freeAtomCrossSection.moniker: freeAtomCrossSection,
                e_critical.moniker: e_critical,
                e_max.moniker: e_max,
                T_effective.moniker: T_effective
            }.get( child.tag, None )
            if _class is None:
                print("Warning: encountered unexpected element '%s' in scatteringAtom!" % child.tag)
            kwargs[ child.tag ] = _class.parseXMLNode( child, xPath, linkData )
        SA = scatteringAtom( **kwargs )
        xPath.pop()
        return SA


class freeAtomCrossSection( physicalQuantityModule.physicalQuantity ):
    moniker = 'freeAtomCrossSection'

class e_critical( physicalQuantityModule.physicalQuantity ):
    moniker = 'e_critical'

class e_max( physicalQuantityModule.physicalQuantity ):
    moniker = 'e_max'

class T_effective( XYsModule.XYs1d ):
    """ In incoherent inelastic sections, each scattering atom using the short collision time (SCT)
    approximation needs an effective temperature. """

    moniker = 'T_effective'
    mutableYUnit = False

    def __init__(self, *args, **kwargs):
        XYsModule.XYs1d.__init__(self, *args, **kwargs)


class S_alpha_beta( ancestryModule.ancestry ):
    """
    Inelastic incoherent section contains S as a function of (alpha,beta,T).
    Currently supports gridded3d, other options may be supported later.
    """
    moniker = 'S_alpha_beta'

    def __init__(self, gridded3d):
        ancestryModule.ancestry.__init__( self )
        self.gridded3d = gridded3d

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s>' % self.moniker]
        xml += self.gridded3d.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        g3d = griddedModule.gridded3d.parseXMLNode( element[0], xPath, linkData )
        Sab = S_alpha_beta( g3d )
        xPath.pop()
        return Sab


def readXML( gndsFile ):
    """
    Read a GNDS/xml file and create a new reactionSuite instance from the result.

    :param gndsFile: path to a GNDS file, as a string.
    :return: reactionSuite instance containing all data from the file.
    """

    from xml.etree import cElementTree
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode

    tsElement = cElementTree.parse( gndsFile ).getroot()
    tsElement = xmlNode( tsElement, xmlNode.etree )
    return thermalScattering.parseXMLNode( tsElement )