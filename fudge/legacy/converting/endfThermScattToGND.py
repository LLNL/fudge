#!/usr/bin/env python

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

""" Translate thermal scattering files from ENDF-6 into GND.

Basic usage:
    python readThermScatt.py <path_to_endf_file>
results:
    test.endf6.xml  (translated into GND/XML)
    test.endf6      (GND/XML translated back into ENDF-6 for comparison

Alternate usage within python:
    >>> from fudge.legacy.converting import readThermScatt
    >>> TS = readThermScatt.parseThermalScattering( "path_to_endf_file" )
    >>> TS.saveToFile( "test.endf6.xml" )
"""

from fudge.gnd import documentation, thermalScattering as TS
from fudge.legacy.converting import endfFileToGNDMisc
from fudge.core.math import table
from fudge.core.math.xData import axes
from pqu import PQU


def readSTable( line, mf7 ):
    line, dat = endfFileToGNDMisc.getTAB1( line, mf7 )
    temps, betas = [PQU( dat['C1'], 'K' )], [dat['C2']]

    # energy interpolation should be 'flat':
    e_interp = dat['interpolationInfo']
    if len(e_interp) > 1: raise Exception("only one interpolation region allowed for S table")
    e_interp = endfFileToGNDMisc.ENDFInterpolationToGND2d( e_interp[0][1] )
    data = zip( *dat['data'] )

    # higher temperatures:
    ntemps = int( dat['L1'] )
    t_interp = []
    for j in range( ntemps ):
        line, dat = endfFileToGNDMisc.getList( line, mf7 )
        temps.append( PQU( dat['C1'], 'K' ) )
        betas.append( dat['C2'] )
        data.append( dat['data'] )
        t_interp.append( endfFileToGNDMisc.ENDFInterpolationToGND2d( dat['L1'] ) )

    # sanity checks: beta and temperature interpolation should be identical for all T
    if len(set(betas)) != 1: raise Exception("inconsistent beta values encountered!")
    beta = betas[0]
    if len(set(t_interp)) > 1: raise Exception("inconsistent temperature interpolations encountered!")
    if t_interp: t_interp = t_interp[0]
    else: t_interp = None

    return line, temps, beta, data, e_interp, t_interp

def readT_effective( line, mf7 ):
    line, dat = endfFileToGNDMisc.getTAB1( line, mf7 )
    e_interp = dat['interpolationInfo']
    if len(e_interp) > 1: raise Exception("only one interpolation region allowed for T_eff data")
    e_interp = endfFileToGNDMisc.ENDFInterpolationToGND2d( e_interp[0][1] )

    t_axes = axes.defaultAxes( labelsUnits = {0: ('temperature','K'), 1: ('t_effective','K')} )
    t_axes[0].interpolation.independent, t_axes[0].interpolation.dependent = e_interp
    return line, TS.T_effective( t_axes, dat['data'], accuracy=0.001 )

def readMF1( mf1 ):
    za, mass, lrp, lfi, nlib, nmod = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( mf1[0] )
    exc, sta, lis, liso, dum, nfor = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( mf1[1] )
    prjmass, dum, lrel, dum, nsub, nver = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( mf1[2] )
    temp, dum, ldrz, dum, nwd, nxc = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( mf1[3] )
    material = mf1[4][:11].strip()
    nwd = int(nwd)
    docs = documentation.documentation( 'endfDoc', '\n'.join( mf1[4:4+nwd] ) )
    return material, za, mass, docs

def readMF7( mt, mf7 ):

    ZA,AWR,LTHR,LAT,LASYM,dum = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( mf7[0], range(2,6) )

    if LTHR==1: # coherent elastic scattering

        line, temps, dummy, data, e_interp, t_interp = readSTable( 1, mf7 )
        Stab = table.table()
        Stab.addColumn( table.columnHeader(0, 'energy_in', units='eV') )
        for i in range(len(temps)):
            Stab.addColumn( table.columnHeader(i+1, 'temp%i' % i, value=temps[i]) )
        for row in zip(*data): Stab.addRow( row )

        if e_interp != ('linear','flat'): raise Exception("unexpected energy interpolation encountered")

        axes_ = axes.defaultAxes( 3,
                labelsUnits={0:('energy_in','eV'), 1:('temperature','K'),2:('S_cumulative','eV*b')} )
        axes_[0].interpolation.independent, axes_[0].interpolation.dependent = e_interp
        axes_[1].interpolation.independent, axes_[1].interpolation.dependent = t_interp

        section = TS.coherentElastic( TS.S_table( axes_, Stab ) )

    elif LTHR == 2: # incoherent elastic
        line, dat  = endfFileToGNDMisc.getTAB1( 1, mf7 )
        SB = PQU( float( dat['C1'] ), 'b' )

        e_interp = dat['interpolationInfo']
        if len(e_interp) > 1: raise Exception("only one interpolation region allowed for T_eff data")
        e_interp = endfFileToGNDMisc.ENDFInterpolationToGND2d( e_interp[0][1] )

        t_axes = axes.defaultAxes( labelsUnits = {0: ('temperature','K'), 1: ('DebyeWallerIntegral','1/eV')} )
        t_axes[0].interpolation.independent, t_axes[0].interpolation.dependent = e_interp
        if len(dat['data'])==2 and dat['data'][0] == dat['data'][1]:
            dat['data'] = dat['data'][0:1]
        DbW = TS.DebyeWaller( t_axes, dat['data'], accuracy=0.001 )
        section = TS.incoherentElastic( SB, DbW )

    elif LTHR == 0: # incoherent inelastic
        line, listy = endfFileToGNDMisc.getList( 1, mf7 )
        LLN, NI, NS, b_n = listy['L1'], listy['NPL'], listy['N2'], listy['data']
        if LLN:
            raise Exception("LLN != 0 not yet handled!")

        # b_n array contains information about principal / secondary scattering atoms
        atoms = [ TS.scatteringAtom( mass = b_n[2], numberPerMolecule = int(b_n[5]), freeAtomCrossSection = PQU(b_n[0],'b'),
            e_critical = b_n[1], e_max = PQU(b_n[3],'eV') ) ]
        for i in range(1,NS+1):
            functionalForm = {0.0: 'SCT', 1.0: 'free_gas', 2.0: 'diffusive_motion'} [b_n[6*i]]
            atoms.append( TS.scatteringAtom( mass = b_n[6*i+2], numberPerMolecule = b_n[6*i+5],
                freeAtomCrossSection = PQU(b_n[6*i+1],'b'), functionalForm=functionalForm ) )

        line, t2header = endfFileToGNDMisc.getTAB2Header( line, mf7 )
        n_betas = int( t2header['NZ'] )

        # read first beta:
        line, temps, beta, data, e_interp, t_interp = readSTable( line, mf7 )
        elist = data[0]
        betas = [[beta, data[1:]]]

        # read in remaining betas:
        for i in range(n_betas-1):
            line, temps_, beta, data, e_interp_, t_interp_ = readSTable( line, mf7 )
            if temps_ != temps or e_interp_ != e_interp or t_interp_ != t_interp: raise Exception("inconsistent values!")
            if data[0] != elist: raise Exception("inconsistent energy list!")
            betas.append( [beta, data[1:]] )

        blist = [b[0] for b in betas]
        T_tables = []
        for i in range(len(temps)):
            Stab = table.table()
            Stab.addColumn( table.columnHeader(0, 'energy_in', units='eV') )
            for j in range(len(blist)):
                Stab.addColumn( table.columnHeader(j+1, 'beta%i' % j, value=blist[j]) )
            data = zip(*([elist]+[beta[1][i] for beta in betas]))
            for row in data: Stab.addRow( row )
            T_tables.append( (temps[i], Stab ) )

        if t_interp is not None:
            axes_ = axes.defaultAxes( 4,
                    labelsUnits={0:('temperature','K'), 1:('energy_in','eV'), 2:('beta',''), 3:('S_alpha_beta','eV*b')} )
            axes_[0].interpolation.independent, axes_[0].interpolation.dependent = t_interp
            axes_[1].interpolation.independent, axes_[1].interpolation.dependent = e_interp
            axes_[2].interpolation.independent, axes_[2].interpolation.dependent = ('linear','flat')
        else:
            axes_ = axes.defaultAxes( 3,
                    labelsUnits={0:('energy_in','eV'), 1:('beta',''), 2:('S_alpha_beta','eV*b')} )
            axes_[0].interpolation.independent, axes_[0].interpolation.dependent = e_interp
            axes_[1].interpolation.independent, axes_[1].interpolation.dependent = ('linear','flat')

        S_tables = TS.S_alpha_beta( axes_, T_tables )

        # last read T_eff tables. There must be at least one (for the principal scattering atom), plus more for
        # any other atoms that specify the short-collision-time approximation:
        line, t_eff = readT_effective( line, mf7 )
        atoms[0].effectiveTemperature = t_eff
        for atom in atoms[1:]:
            if atom.functionalForm == 'SCT':
                line, t_eff = readT_effective( line, mf7 )
                atom.effectiveTemperature = t_eff

        section = TS.incoherentInelastic( S_tables, calculatedAtThermal = LAT, asymmetric = LASYM,
                atoms = atoms )

    if line != len(mf7): raise Exception("Trailing data left in MT%i MF7!" % mt)
    return LTHR, section

def parseThermalScattering( filename ):
    header, mat, mtmfs = endfFileToGNDMisc.parseENDFByMT_MF( filename )
    material, za, mass, doc = readMF1( mtmfs[451][1] )

    thermScat = TS.thermalScattering( material, MAT=mat, mass=mass, documentation=doc )

    if 2 in mtmfs:  # elastic
        LTHR, section = readMF7( 2, mtmfs[2][7] )
        if LTHR==1:
            thermScat.coherentElastic = section
        elif LTHR==2:
            thermScat.incoherentElastic = section

    if 4 in mtmfs:  # inelastic
        LTHR, section = readMF7( 4, mtmfs[4][7] )
        thermScat.incoherentInelastic = section

    return thermScat

if __name__ == '__main__':
    import sys
    thermScat = parseThermalScattering( sys.argv[1] )
    thermScat.saveToFile( "test.endf6.xml" )
    with open('test.endf6','w') as fout:
        fout.write( thermScat.toENDF6() )
