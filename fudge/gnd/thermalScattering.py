# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

__metaclass__ = type

from fudge.core.math.xData import XYs
from fudge.legacy.converting import endfFormats
import itertools

coherentElasticToken = 'coherentElastic'
incoherentElasticToken = 'incoherentElastic'
incoherentInelasticToken = 'incoherentInelastic'

class thermalScattering:
    def __init__(self, material=None, MAT=None, mass=None, documentation=None, coherentElastic=None,
            incoherentElastic=None, incoherentInelastic=None):
        self.material = material
        self.MAT = MAT
        self.mass = mass
        self.documentation = documentation
        self.coherentElastic = coherentElastic
        self.incoherentElastic = incoherentElastic
        self.incoherentInelastic = incoherentInelastic

    def toXMLList( self, indent='' ):
        xml = ['<?xml version="1.0" encoding="UTF-8"?>']
        xml.append( indent+'<thermalScattering material="%s" MAT="%s" mass="%s">' % (self.material,self.MAT,self.mass) )
        if self.documentation is not None: xml += self.documentation.toXMLList(indent+'  ')
        for ts in (coherentElasticToken,incoherentElasticToken,incoherentInelasticToken):
            if getattr(self,ts) is not None:
                xml += getattr(self,ts).toXMLList(indent+'  ')
        xml.append(indent+'</thermalScattering>')
        return xml

    def saveToFile( self, filename ):
        with open(filename,'w') as fout:
            fout.write( '\n'.join( self.toXMLList() ) )

    def toENDF6( self, flags = {}, verbosityIndent = '' ):
        endfMFList = { 1 : { 451 : [] }, 7 : {} }
        targetInfo = {'ZA':self.MAT + 100, 'mass':self.mass}
        MAT = self.MAT
        NSUB, NVER = 12, 7  # 12: thermal scattering sub-library, 7: ENDF/B-VII

        for subsection in (coherentElasticToken,incoherentElasticToken,incoherentInelasticToken):
            if getattr(self, subsection) is not None:
                getattr(self, subsection).toENDF6( endfMFList, flags, targetInfo, verbosityIndent )

        endfDoc = self.documentation.getLines()
        docHeader = [
                endfFormats.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], -1, 0, 0, 0 ),
                endfFormats.endfHeadLine( 0, 0, 0, 0, 0, 6 ),    # ENDF-6
                endfFormats.endfHeadLine( 1.0, 0, 1, 0, NSUB, NVER ),
                endfFormats.endfHeadLine( 0.0, 0.0, 0, 0, len(endfDoc), 0 ) ]
        endfMFList[1][451] = docHeader + endfDoc

        return endfFormats.endfMFListToFinalFile( endfMFList, MAT, lineNumbers=True )


class coherentElastic:
    def __init__(self, S_table):
        self.S_table = S_table

    def toXMLList( self, indent='' ):
        xml = [indent+'<%s>' % coherentElasticToken]
        xml += self.S_table.toXMLList( indent=indent+'  ' )
        xml[-1] += '</%s>' % coherentElasticToken
        return xml

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        LTHR = 1
        endf = [endfFormats.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
        endf += self.S_table.toENDF6( flags, targetInfo, verbosityIndent = '' )
        endfMFList[7][2] = endf + [99999]

class S_table:
    """ For elastic coherent, cumulative structure factor 'S' is given as function of incident energy and temp. """ 
    def __init__(self, axes, table):
        self.axes = axes
        self.table = table

    def toXMLList( self, indent='' ):
        xml = [indent+'<S_table>']
        xml += self.axes.toXMLList( indent+'  ' )
        xml += self.table.toXMLList( indent+'  ' )
        xml[-1] += '</S_table>'
        return xml

    def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):

        # first temperature also has energy list:
        data = zip(*self.table.data)
        Tlist = [float( temp['value'].getValueAs('K') ) for temp in self.table.columns[1:] ]

        LT = len(Tlist)-1
        # first temperature includes the energy list:
        endf = [endfFormats.endfHeadLine( Tlist[0], 0, LT, 0, 1, len(data[0]) ) ]
        independentInterp = endfFormats.twoAxesToENDFInterpolation( self.axes, 0 )
        dependentInterp = endfFormats.twoAxesToENDFInterpolation( self.axes, 1 )
        #endf += endfFormats.endfInterpolationList( (len(data[0])/2, independentInterp ) )
        endf += ['%11i%11i%44s' % (len(data[0]), independentInterp, '' )]    # no trailing zeros
        endf += endfFormats.endfDataList( list( itertools.chain( *zip(data[0], data[1]) ) ) )

        # remaining temperatures:
        for T, datList in zip( Tlist[1:], data[2:] ):
            endf += [endfFormats.endfHeadLine( T, 0, dependentInterp, 0, len(datList), 0 )]
            endf += endfFormats.endfDataList( datList )
        return endf


class incoherentElastic:
    def __init__(self, characteristicCrossSection, DebyeWaller):
        self.characteristicCrossSection = characteristicCrossSection
        self.DebyeWaller = DebyeWaller

    def toXMLList( self, indent='' ):
        xml = [indent+'<%s characteristicCrossSection="%s">' % (incoherentElasticToken, 
            self.characteristicCrossSection) ]
        xml += self.DebyeWaller.toXMLList( indent=indent+'  ' )
        xml[-1] += '</%s>' % incoherentElasticToken
        return xml

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        LTHR = 2
        endf = [endfFormats.endfHeadLine( ZAM, AWT, LTHR, 0, 0, 0 )]
        targetInfo['characteristicCrossSection'] = self.characteristicCrossSection.getValueAs('b')
        endf += self.DebyeWaller.toENDF6( flags, targetInfo, verbosityIndent = '' )
        del targetInfo['characteristicCrossSection']
        endfMFList[7][2] = endf + [99999]

class DebyeWaller( XYs.XYs ):
    """For incoherent elastic sections, all we need is a characteristic cross section and a
    temperature-dependent list of the Debye-Waller integral """

    mutableYUnit = False

    def __init__(self, *args, **kwargs):
        XYs.XYs.__init__(self, *args, **kwargs)
        self.moniker = "DebyeWaller"

    def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
        NR = 1; NP = len(self)
        endf = [endfFormats.endfHeadLine( targetInfo['characteristicCrossSection'], 0, 0, 0, NR, NP )]
        interp = endfFormats.twoAxesToENDFInterpolation( self.axes, 0 )
        endf += ['%11i%11i%44s' % (len(self), interp, '')]
        endf += endfFormats.endfDataList( list( itertools.chain( *self.copyDataToXYs() ) ) )
        return endf


class incoherentInelastic:
    def __init__(self, S_tables, calculatedAtThermal = False, asymmetric = False, atoms = []):
        self.scatteringAtoms = atoms
        self.S_tables = S_tables
        self.calculatedAtThermal = calculatedAtThermal
        self.asymmetric = asymmetric

    def toXMLList( self, indent='' ):
        xml = [indent+'<%s' % incoherentInelasticToken]
        if self.calculatedAtThermal: xml[0] += ' calculatedAtThermal="true"'
        if self.asymmetric: xml[0] += ' asymmetric="true"'
        xml[0] += '>'
        xml += [indent+'  <scatteringAtoms>']
        for atom in self.scatteringAtoms:
            xml += atom.toXMLList( indent+'    ' )
        xml[-1] += '</scatteringAtoms>'
        xml += self.S_tables.toXMLList(indent+'  ')
        xml[-1] += '</%s>' % incoherentInelasticToken
        return xml

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ):
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        LTHR, LAT, LASYM = 0, self.calculatedAtThermal, self.asymmetric
        endf = [endfFormats.endfHeadLine( ZAM, AWT, LTHR, LAT, LASYM, 0 )]
        # describe scattering atoms:
        LLN, NS = 0, len( self.scatteringAtoms ) - 1
        NI = 6*(NS+1)
        endf += [endfFormats.endfHeadLine( 0, 0, LLN, 0, NI, NS )]
        # principal scattering atom:
        atom = self.scatteringAtoms[0]
        endf += [endfFormats.endfDataLine( [atom.freeAtomCrossSection.getValueAs('b'),
            atom.e_critical, atom.mass, atom.e_max.getValueAs('eV'), 0, atom.numberPerMolecule] )]
        for atom in self.scatteringAtoms[1:]:
            a1 = {'SCT':0.0, 'free_gas':1.0, 'diffusive_motion':2.0} [ atom.functionalForm ]
            endf += [endfFormats.endfDataLine( [a1, atom.freeAtomCrossSection.getValueAs('b'),
                atom.mass, 0, 0, atom.numberPerMolecule] )]

        # convert data form: sort first by beta, then E, then T
        betas = [beta['value'] for beta in self.S_tables[0][1].columns[1:] ]
        NR = 1; NB = len(betas)
        endf += [endfFormats.endfHeadLine( 0, 0, 0, 0, NR, NB )]
        #endf += endfFormats.endfInterpolationList( (NB, 4) )
        endf += ['%11i%11i%44s' % ( NB, 4, '' )]

        Tlist = [temp[0].getValueAs('K') for temp in self.S_tables]
        elist = self.S_tables[0][1] [:,0]
        LT = len(Tlist)-1

        if LT:
            T_interp = endfFormats.twoAxesToENDFInterpolation( self.S_tables.axes, 0 )
            alpha_interp = endfFormats.twoAxesToENDFInterpolation( self.S_tables.axes, 1 )
            beta_interp = endfFormats.twoAxesToENDFInterpolation( self.S_tables.axes, 2 )
        else:
            T_interp = None
            alpha_interp = endfFormats.twoAxesToENDFInterpolation( self.S_tables.axes, 0 )
            beta_interp = endfFormats.twoAxesToENDFInterpolation( self.S_tables.axes, 1 )


        for index, beta in enumerate( betas ):
            data = [dat[1][:,index+1] for dat in self.S_tables]

            endf += [endfFormats.endfHeadLine( Tlist[0], beta, LT, 0, 1, len(data[0]) )]
            #endf += endfFormats.endfInterpolationList( (len(data[0])/2, alpha_interp ) )
            endf += ['%11i%11i%44s' % (len(data[0]), alpha_interp, '' )]    # no trailing zeros
            # For each beta, the first temperature needs to include the energy list:
            endf += endfFormats.endfDataList( list( itertools.chain( *zip(elist, data[0]) ) ) )

            # remaining temperatures:
            for T, datList in zip( Tlist[1:], data[1:] ):
                endf += [endfFormats.endfHeadLine( T, beta, T_interp, 0, len(datList), 0 )]
                endf += endfFormats.endfDataList( datList )

        for atom in self.scatteringAtoms:
            if atom.effectiveTemperature is not None:
                endf += atom.effectiveTemperature.toENDF6( flags, targetInfo, verbosityIndent = '' )
        endfMFList[7][4] = endf + [99999]

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

    def toXMLList( self, indent='' ):
        xml = [indent+'<atom mass="%s" numberPerMolecule="%i" freeAtomCrossSection="%s"' % (self.mass,
            self.numberPerMolecule, self.freeAtomCrossSection)]
        for attribute in ('e_critical','e_max','functionalForm'):
            if getattr(self, attribute) is not None:
                xml[-1] += ' %s="%s"' % (attribute, getattr(self,attribute))
        if self.effectiveTemperature is not None:
            xml[-1] += '>'
            xml += self.effectiveTemperature.toXMLList( indent=indent+'  ' )
            xml[-1] += '</atom>'
        else:
            xml[-1] += '/>'
        return xml

class T_effective( XYs.XYs ):
    """ In incoherent inelastic sections, each scattering atom using the short collision time (SCT)
    approximation needs an effective temperature. """

    mutableYUnit = False

    def __init__(self, *args, **kwargs):
        XYs.XYs.__init__(self, *args, **kwargs)
        self.moniker = "T_effective"

    def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
        NR = 1; NP = len(self)
        endf = [endfFormats.endfHeadLine( 0, 0, 0, 0, NR, NP )]
        interp = endfFormats.twoAxesToENDFInterpolation( self.axes, 0 )
        endf += ['%11i%11i%44s' % (len(self), interp, '')]
        endf += endfFormats.endfDataList( list( itertools.chain( *self.copyDataToXYs() ) ) )
        return endf

class S_alpha_beta:
    """ Inelastic incoherent section contains one S_ab table per temperature """
    def __init__(self, axes, temperatures):
        self.axes = axes
        self.temperatures = temperatures

    def __getitem__(self, index): return self.temperatures[index]

    def toXMLList( self, indent='' ):
        xml = [indent+'<S_alpha_beta nTemperatures="%i">' % len(self.temperatures)]
        xml += self.axes.toXMLList( indent+'  ' )
        for temp, s_table in self.temperatures:
            xml += [indent+'  <temperature value="%s">' % temp]
            xml += s_table.toXMLList( indent+'    ' )
            xml[-1] += '</temperature>'
        xml[-1] += '</S_alpha_beta>'
        return xml

