#!/usr/bin/env python

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
ENDF ITYPE=2 files store thermal neutron scattering data.
This module is for translating ITYPE=2 from ENDF-6 into GNDS.

These methods are meant to be called from endfFileToGNDS.
"""

import numpy
from fudge.gnds import documentation, thermalScattering as TS
from fudge.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc
from xData import axes as axesModule
from xData import array as arrayModule
from xData import gridded as griddedModule
from xData import values as valuesModule
from xData import standards as standardsModule

def ITYPE_2( MTDatas, info, verbose ):
    errors = []
    material = info.materialName
    mass = TS.mass( info.targetMass, 'amu' )
    emax = TS.cutoffEnergy( info.EMAX, 'eV' )
    doc = info.documentation

    thermScat = TS.thermalScattering( cutoffEnergy=emax, material=material, MAT=info.MAT, mass=mass, documentation=doc )

    if 2 in MTDatas:  # elastic
        LTHR, section = readMF7( 2, MTDatas[2][7] )
        if LTHR==1:
            thermScat.coherentElastic = section
        elif LTHR==2:
            thermScat.incoherentElastic = section

    if 4 in MTDatas:  # inelastic
        LTHR, section = readMF7( 4, MTDatas[4][7] )
        thermScat.incoherentInelastic = section

    return { 'reactionSuite': thermScat, 'errors': errors, 'covarianceSuite': None }

def readSTable( line, mf7 ):
    line, dat = endfFileToGNDSMisc.getTAB1( line, mf7 )
    temps, betas = [dat['C1']], [dat['C2']]

    # energy interpolation should be 'flat':
    e_interp = dat['interpolationInfo']
    if len(e_interp) > 1: raise ValueError("only one interpolation region allowed for S table")
    e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d( e_interp[0][1] )
    energies, data = map(list, zip( *dat['data'] ) )

    # higher temperatures:
    ntemps = int( dat['L1'] )
    t_interp = []
    for j in range( ntemps ):
        line, dat = endfFileToGNDSMisc.getList( line, mf7 )
        temps.append( dat['C1'] )
        betas.append( dat['C2'] )
        data.extend( dat['data'] )
        t_interp.append( endfFileToGNDSMisc.ENDFInterpolationToGNDS1d( dat['L1'] ) )

    # sanity checks: beta and temperature interpolation should be identical for all T
    if len(set(betas)) != 1: raise ValueError("inconsistent beta values encountered!")
    beta = betas[0]
    if len(set(t_interp)) > 1: raise ValueError("inconsistent temperature interpolations encountered!")
    if t_interp: t_interp = t_interp[0]
    else: t_interp = None

    return line, temps, energies, beta, data, e_interp, t_interp

def readT_effective( line, mf7 ):
    line, dat = endfFileToGNDSMisc.getTAB1( line, mf7 )
    e_interp = dat['interpolationInfo']
    if len(e_interp) > 1: raise ValueError("only one interpolation region allowed for T_eff data")
    e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d( e_interp[0][1] )

    t_axes = axesModule.axes( labelsUnits = {1: ('temperature','K'), 0: ('t_effective','K')} )
    return line, TS.T_effective( dat['data'], axes=t_axes, interpolation=e_interp )

def readMF7( mt, mf7 ):

    ZA,AWR,LTHR,LAT,LASYM,dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( mf7[0], range(2,6) )

    if LTHR==1: # coherent elastic scattering

        line, temps, energies, dummy, data, e_interp, t_interp = readSTable( 1, mf7 )
        temps, energies = map( valuesModule.values, (temps, energies) )

        if e_interp != 'flat': raise ValueError("unexpected energy interpolation encountered")

        axes = axesModule.axes( rank=3 )
        axes[2] = axesModule.grid( 'temperature', 2, 'K', axesModule.pointsGridToken, values=temps, interpolation=t_interp )
        axes[1] = axesModule.grid( 'energy_in', 1, 'eV', axesModule.pointsGridToken, values=energies, interpolation=e_interp )
        axes[0] = axesModule.axis( 'S_cumulative', 0, 'eV*b' )
        array = arrayModule.full( shape=(len(temps),len(energies)), data=data )
        Stab = griddedModule.gridded2d( axes, array, label="eval" )

        section = TS.coherentElastic( TS.S_table( Stab ) )

    elif LTHR == 2: # incoherent elastic
        line, dat  = endfFileToGNDSMisc.getTAB1( 1, mf7 )
        SB = TS.characteristicCrossSection( float( dat['C1'] ), 'b' )

        e_interp = dat['interpolationInfo']
        if len(e_interp) > 1: raise ValueError("only one interpolation region allowed for T_eff data")
        e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d( e_interp[0][1] )

        t_axes = axesModule.axes( labelsUnits = {1: ('temperature','K'), 0: ('DebyeWallerIntegral','1/eV')} )
        if len(dat['data'])==2 and dat['data'][0] == dat['data'][1]:
            dat['data'] = dat['data'][0:1]
        DbW = TS.DebyeWaller( dat['data'], axes=t_axes, interpolation=e_interp )
        section = TS.incoherentElastic( SB, DbW )

    elif LTHR == 0: # incoherent inelastic
        line, listy = endfFileToGNDSMisc.getList( 1, mf7 )
        LLN, NI, NS, b_n = listy['L1'], listy['NPL'], listy['N2'], listy['data']
        if LLN:
            raise NotImplementedError("LLN != 0 not yet handled!")

        # b_n array contains information about principal / secondary scattering atoms
        atoms = [ TS.scatteringAtom( label = "0", numberPerMolecule = int(b_n[5]),
            mass = TS.mass( b_n[2], 'amu' ),
            freeAtomCrossSection = TS.freeAtomCrossSection( b_n[0],'b'),
            e_critical = TS.e_critical( b_n[1], 'eV' ),
            e_max = TS.e_max( b_n[3],'eV') ) ]
        for idx in range(1,NS+1):
            functionalForm = {0.0: 'SCT', 1.0: 'free_gas', 2.0: 'diffusive_motion'} [b_n[6*idx]]
            atoms.append( TS.scatteringAtom( label = str(idx), numberPerMolecule = b_n[6*idx+5],
                mass = TS.mass( b_n[6*idx+2], 'amu' ),
                freeAtomCrossSection = TS.freeAtomCrossSection( b_n[6*idx+1],'b' ),
                functionalForm=functionalForm ) )

        line, t2header = endfFileToGNDSMisc.getTAB2Header( line, mf7 )
        n_betas = int( t2header['NZ'] )

        # read first beta:
        line, temps, alphas, beta, data, a_interp, t_interp = readSTable( line, mf7 )
        betas = [beta]

        # read in remaining betas:
        for idx in range(n_betas-1):
            line, temps_, alphas_, beta, data_, a_interp_, t_interp_ = readSTable( line, mf7 )
            if (temps_ != temps or alphas_ != alphas or a_interp_ != a_interp or t_interp_ != t_interp):
                raise ValueError("inconsistent values!")
            betas.append( beta )
            data.extend( data_ )

        temps, betas, alphas = map( valuesModule.values, (temps, betas, alphas) )

        axes = axesModule.axes( rank=4 )
        axes[3] = axesModule.grid('temperature', 3, 'K', axesModule.pointsGridToken, values=temps, interpolation=t_interp )
        axes[2] = axesModule.grid('beta', 2, '', axesModule.pointsGridToken, values=betas,
            interpolation=standardsModule.interpolation.flatToken )
        axes[1] = axesModule.grid('alpha', 1, '', axesModule.pointsGridToken, values=alphas, interpolation=a_interp )
        axes[0] = axesModule.axis('S_alpha_beta', 0, 'eV*b' )

        arrayShape_orig = (len(betas),len(temps),len(alphas))
        data = numpy.array( data ).reshape( arrayShape_orig )
        # ENDF data are stored as 'beta,T,alpha'. Switch order to 'T,beta,alpha':
        data = numpy.transpose( data, axes=(1,0,2) )
        arrayShape_new = (len(temps),len(betas),len(alphas))
        array = arrayModule.full( shape=arrayShape_new, data=data.flatten() )
        Stab = griddedModule.gridded3d( axes, array, label="eval" )

        S_tables = TS.S_alpha_beta( Stab )

        # last read T_eff tables. There must be at least one (for the principal scattering atom), plus more for
        # any other atoms that specify the short-collision-time approximation:
        line, t_eff = readT_effective( line, mf7 )
        atoms[0].T_effective = t_eff
        for atom in atoms[1:]:
            if atom.functionalForm == 'SCT':
                line, t_eff = readT_effective( line, mf7 )
                atom.T_effective = t_eff

        section = TS.incoherentInelastic( S_tables, calculatedAtThermal = LAT, asymmetric = LASYM,
                atoms = atoms )

    if line != len(mf7): raise ValueError("Trailing data left in MT%i MF7!" % mt)
    return LTHR, section

