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

from xData import XYs as XYsModule

coherentElasticToken = 'coherentElastic'
incoherentElasticToken = 'incoherentElastic'
incoherentInelasticToken = 'incoherentInelastic'

class thermalScattering:
    def __init__(self, material=None, MAT=None, mass=None, emax=None, documentation=None, coherentElastic=None,
            incoherentElastic=None, incoherentInelastic=None):
        self.material = material
        self.MAT = MAT
        self.mass = mass
        self.emax = emax
        self.documentation = documentation
        self.coherentElastic = coherentElastic
        self.incoherentElastic = incoherentElastic
        self.incoherentInelastic = incoherentInelastic

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = ['<?xml version="1.0" encoding="UTF-8"?>']
        xml.append( indent+'<thermalScattering material="%s" MAT="%s" mass="%s" emax="%s">'
                    % (self.material,self.MAT,self.mass,self.emax) )
        if self.documentation is not None: xml += self.documentation.toXMLList( indent2, **kwargs )
        for ts in (coherentElasticToken,incoherentElasticToken,incoherentInelasticToken):
            if getattr(self,ts) is not None:
                xml += getattr(self,ts).toXMLList( indent2, **kwargs )
        xml.append(indent+'</thermalScattering>')
        return xml

    def check( self, **kwargs ) :
        """
        Check all data in the reactionSuite, returning a list of warnings. 
        """

        from fudge.gnd import warning

        options = { }
        for key in kwargs:
            if key in options: options[key] = kwargs[key]
            else: raise KeyError, "check() received unknown keyword argument '%s'" % key

        # assemble some useful info, to be handed down to children's check() functions:
        info = { 'reactionSuite': self, }
        info.update( options )

        warnings = []


        result = warning.context('ReactionSuite: %s + %s' % (self.projectile, self.target), warnings)
        result.info = info
        return result


    def saveToFile( self, filename ):
        with open(filename,'w') as fout:
            fout.write( '\n'.join( self.toXMLList() ) )


class coherentElastic:
    def __init__(self, S_table):
        self.S_table = S_table

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s>' % coherentElasticToken]
        xml += self.S_table.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % coherentElasticToken
        return xml


class S_table:
    """ For elastic coherent, cumulative structure factor 'S' is given as function of incident energy and temp. """ 
    def __init__(self, form=None):
        self.forms = []
        if form is not None: self.forms.append( form )

    def __getitem__(self, idx): return self.forms[idx]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<S_table>']
        for form in self:
            xml += form.toXMLList( indent2, **kwargs )
        xml[-1] += '</S_table>'
        return xml


class incoherentElastic:
    def __init__(self, characteristicCrossSection, DebyeWaller):
        self.characteristicCrossSection = characteristicCrossSection
        self.DebyeWaller = DebyeWaller

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<%s characteristicCrossSection="%s">' % (incoherentElasticToken, 
            self.characteristicCrossSection) ]
        xml += self.DebyeWaller.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % incoherentElasticToken
        return xml


class DebyeWaller( XYsModule.XYs1d ):
    """For incoherent elastic sections, all we need is a characteristic cross section and a
    temperature-dependent list of the Debye-Waller integral """

    mutableYUnit = False

    def __init__(self, *args, **kwargs):
        XYsModule.XYs1d.__init__(self, *args, **kwargs)
        self.moniker = "DebyeWaller"


class incoherentInelastic:
    def __init__(self, S_alpha_beta, calculatedAtThermal = False, asymmetric = False, atoms = []):
        self.scatteringAtoms = atoms
        self.S_alpha_beta = S_alpha_beta
        self.calculatedAtThermal = calculatedAtThermal
        self.asymmetric = asymmetric

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xml = [indent+'<%s' % incoherentInelasticToken]
        if self.calculatedAtThermal: xml[0] += ' calculatedAtThermal="true"'
        if self.asymmetric: xml[0] += ' asymmetric="true"'
        xml[0] += '>'
        xml += [indent+'  <scatteringAtoms>']
        for atom in self.scatteringAtoms:
            xml += atom.toXMLList( indent3, **kwargs )
        xml[-1] += '</scatteringAtoms>'
        xml += self.S_alpha_beta.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % incoherentInelasticToken
        return xml


class scatteringAtom:
    """ Inelastic incoherent scattering requires a description of each atom that is scattered off """
    def __init__(self, mass, numberPerMolecule, freeAtomCrossSection, e_critical = None, e_max = None,
            functionalForm = None, effectiveTemperature = None ):
        self.mass = mass
        self.numberPerMolecule = numberPerMolecule
        self.freeAtomCrossSection = freeAtomCrossSection
        self.e_critical = e_critical
        self.e_max = e_max
        self.functionalForm = functionalForm
        self.effectiveTemperature = effectiveTemperature

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<atom mass="%s" numberPerMolecule="%i" freeAtomCrossSection="%s"' % (self.mass,
            self.numberPerMolecule, self.freeAtomCrossSection)]
        for attribute in ('e_critical','e_max','functionalForm'):
            if getattr(self, attribute) is not None:
                xml[-1] += ' %s="%s"' % (attribute, getattr(self,attribute))
        if self.effectiveTemperature is not None:
            xml[-1] += '>'
            xml += self.effectiveTemperature.toXMLList( indent2, **kwargs )
            xml[-1] += '</atom>'
        else:
            xml[-1] += '/>'
        return xml

class T_effective( XYsModule.XYs1d ):
    """ In incoherent inelastic sections, each scattering atom using the short collision time (SCT)
    approximation needs an effective temperature. """

    mutableYUnit = False

    def __init__(self, *args, **kwargs):
        XYsModule.XYs1d.__init__(self, *args, **kwargs)
        self.moniker = "T_effective"


class S_alpha_beta:
    """ Inelastic incoherent section contains one S_ab table per temperature """
    def __init__(self, form=None):
        self.forms = []
        if form is not None: self.forms.append( form )

    def __getitem__(self, index): return self.forms[index]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [indent+'<S_alpha_beta>']
        for form in self:
            xml += form.toXMLList( indent2, **kwargs )
        xml[-1] += '</S_alpha_beta>'
        return xml

