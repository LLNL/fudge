# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
ENDF ITYPE=2 files store thermal neutron scattering data.
This module is for translating ITYPE=2 from ENDF-6 into GNDS.

These methods are meant to be called from endfFileToGNDS.
"""

import os
import numpy

from pqu import PQU as PQUModule

from xData import axes as axesModule
from xData import xDataArray as arrayModule
from xData import gridded as griddedModule
from xData import values as valuesModule
from xData import XYs as XYs1dModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.families import unorthodox as unorthodoxModule

from fudge import physicalQuantity as physicalQuantityModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

from fudge.reactions import reaction as reactionModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import coherentElastic as coherentElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentElastic as incoherentElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentInelastic as incoherentInelasticModule
import fudge.productData.distributions.reference as referenceModule

from fudge import institution as institutionModule
from fudge import outputChannel as outputChannelModule
from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

from .. import endf_endl as endf_endlModule
from .. import toGNDSMisc as toGNDSMiscModule

def ITYPE_2( fileName, MAT, MTDatas, info, evaluation, verbose ) :
    """Parses ENDF Thermal neutron scattering law data."""

    errors = []
    EMin = 1e-5
    EMax = info.EMAX
    if( EMax <= EMin ) :
        print( '    WARNING: EMax = %s, setting to 5.0 eV' % EMax )
        EMax = 5.0
    if( EMax > 5 ) :
        print( '    WARNING: EMax = %s is unphysical, setting to 5.0 eV' % EMax )
        EMax = 5.0

    projectile = toGNDSMiscModule.getTypeNameGamma( info, 1 )
    neutronMass = projectile.getMass( 'amu' )

    if( evaluation is None ) : evaluation = "%s-%d.%d" % ( info.library, info.NVER, info.LREL )

    projectileDomain = stylesModule.projectileEnergyDomain( EMin, EMax, 'eV' )
    axes = crossSectionModule.defaultAxes( 'eV' )
    domainXYs1d = XYs1dModule.XYs1d( [ [ EMin, 1.0 ], [ EMax, 1.0 ] ], axes = axes )

    evaluatedStyle = stylesModule.evaluated( info.PoPsLabel, '',
            physicalQuantityModule.temperature( PQUModule.pqu_float.surmiseSignificantDigits( info.targetTemperature ), 'K' ),
            projectileDomain, info.library, info.libraryVersion, date = info.Date )

    fileBaseName = os.path.basename( fileName )[4:][:-5]
    targetID = None

    if( fileBaseName[0].isdigit( ) ) :
        Z, Symbol, A = fileBaseName.split( '_' )
        ZA = int( Z ) * 1000 + int( A )
        targetID = 'tnsl-%s%d' % ( Symbol, int( A ) )
    elif( ( 'para' in fileBaseName ) or ( 'ortho' in fileBaseName ) ) :
        ZA = 1001
        if( fileBaseName[-1] == 'D' ) : ZA = 1002
        targetID = fileBaseName
    elif( 'graphite' in fileBaseName ) :
        ZA = 6012
        targetID = fileBaseName
    elif( 'Be-metal' in fileBaseName ) :
        ZA = 4009
        targetID = fileBaseName
    elif( 'CH4' in fileBaseName ) :
        target = toGNDSMiscModule.getTypeNameGamma( info, 6012 )
        mass = massModule.double( info.PoPsLabel, neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        target.mass.add( mass )
        ZA = 1001
        targetID = fileBaseName
    elif( 'SiO2' in fileBaseName ) :
        target = toGNDSMiscModule.getTypeNameGamma( info, 14028 )
        mass = massModule.double( info.PoPsLabel, neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        target.mass.add( mass )
        ZA = 8016
        targetID = fileBaseName
    elif( 'benzene' in fileBaseName or 'benzine' in fileBaseName ) :    # ENDF-VII.1 had spelling error
        target = toGNDSMiscModule.getTypeNameGamma( info, 1001 )
        mass = massModule.double( info.PoPsLabel, neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        target.mass.add( mass )
        ZA = 6012
        targetID = fileBaseName
    elif( 'in' in fileBaseName ) :
        targetSymbol = fileBaseName.split( 'in' )[0]
        ZA = { 'H' : 1001, 'D' : 1002, 'Li' : 3006, 'Be' : 4009, 'C' : 6012, 'N' : 7014, 'O' : 8016, 'Al' : 13027, 'Si' : 14028, 'Y' : 39089, 'Zr' : 40090, 'U' : 92238 }[targetSymbol]
        targetID = fileBaseName
    elif( 'BeO' in fileBaseName ) :  # ENDF-VII.0 and earlier didn't break up into Be and O components
        target = toGNDSMiscModule.getTypeNameGamma( info, 8016 )
        mass = massModule.double( info.PoPsLabel, neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        target.mass.add( mass )
        ZA = 4009
        targetID = fileBaseName
    else :
        # FIXME try using the MAT number to get targetID? Unfortunately libraries use different conventions for TNSL MAT numbers...
        raise Exception("Cannot determine TNSL target ZA or id from filename %s" % fileName)

    target = toGNDSMiscModule.getTypeNameGamma( info, ZA )
    mass = massModule.double( info.PoPsLabel, neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) )
    target.mass.add( mass )
    if( targetID is None ) : targetID = "tnsl-" + target.id
    info.PoPs.add(unorthodoxModule.particle(targetID))

    interaction = reactionSuiteModule.Interaction.TNSL
    evaluatedStyle.documentation.endfCompatible.body = info.documentation
    reactionSuite = reactionSuiteModule.reactionSuite( info.projectile, targetID, evaluation, style = evaluatedStyle, interaction = interaction,
            formatVersion = info.formatVersion, MAT = MAT, PoPs = info.PoPs )

    for MT in [ 2, 4 ] :
        if( MT in MTDatas ) :                                # elastic
            forms = readMF7( info, MT, MTDatas[MT][7] )

            for form in forms:
                outputChannel = outputChannelModule.outputChannel( outputChannelModule.Genre.twoBody, process = form.process )
                outputChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, 0.0, domainXYs1d ) )

                product = toGNDSMiscModule.newGNDSParticle( info, projectile, domainXYs1d )
                outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )

                reaction = reactionModule.reaction( outputChannel.genre, ENDF_MT = MT )
                endf_endlModule.setReactionsOutputChannelFromOutputChannel( info, reaction, outputChannel )
                reaction.doubleDifferentialCrossSection.add( form )
                reactionSuite.reactions.add( reaction )

                reaction.crossSection.add( crossSectionModule.thermalNeutronScatteringLaw1d( link = form, label = info.style, relative = True ) )
                product.distribution.add( referenceModule.thermalNeutronScatteringLaw( link = form, label = info.style, relative = True ) )

    # add applicationData section for help converting back to ENDF-6
    info.ENDFconversionFlags.add(reactionSuite, 'MAT=%d,ZA=%d' % (MAT, ZA))
    LLNLdata = institutionModule.institution("LLNL")
    LLNLdata.append( info.ENDFconversionFlags )
    reactionSuite.applicationData.add(LLNLdata)

    if len(info.doRaise) > 0:
        info.logs.write( '\nRaising due to following errors:\n' )
        for err in info.doRaise : info.logs.write( '    ' + err + '\n' )
        raise Exception( 'len( info.doRaise ) > 0' )

    return { 'reactionSuite' : reactionSuite, 'errors' : errors, 'covarianceSuite' : None }

def readSTable( line, MF7 ) :

    line, dat = endfFileToGNDSMisc.getTAB1(line, MF7)
    temps, betas = [dat['C1']], [dat['C2']]

    # energy interpolation should be 'flat':
    e_interp = dat['interpolationInfo']
    if( len( e_interp ) > 1 ) : raise ValueError( "only one interpolation region allowed for S table" )
    e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(e_interp[0][1])
    energies, data = map(list, zip( *dat['data'] ) )

    # higher temperatures:
    ntemps = int( dat['L1'] )
    t_interp = []
    for j in range( ntemps ) :
        line, dat = endfFileToGNDSMisc.getList(line, MF7)
        temps.append( dat['C1'] )
        betas.append( dat['C2'] )
        data.extend( dat['data'] )
        t_interp.append(endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(dat['L1']))

    # sanity checks: beta and temperature interpolation should be identical for all T
    if( len( set( betas ) ) != 1 ) : raise ValueError( "inconsistent beta values encountered!" )
    beta = betas[0]
    if( len( set( t_interp ) ) > 1 ) : raise ValueError( "inconsistent temperature interpolations encountered!" )
    if( t_interp ) :
        t_interp = t_interp[0]
    else :
        t_interp = None

    return line, temps, energies, beta, data, e_interp, t_interp

def readT_effective( line, MF7 ) :

    line, data = endfFileToGNDSMisc.getTAB1(line, MF7)
    interpolation = data['interpolationInfo']
    if( len( interpolation ) > 1 ) : raise ValueError( "Only one interpolation region allowed for T_eff data" )
    interpolation = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(interpolation[0][1])

    t_axes = axesModule.axes( labelsUnits = { 1 : ( 'temperature', 'K' ), 0 : ( 't_effective', 'K' ) } )
    function1d = XYs1dModule.XYs1d( data = data['data'], axes = t_axes, interpolation = interpolation )
    return( line, incoherentInelasticModule.T_effective( function1d ) )

def readMF7( info, MT, MF7 ) :

    ZA, AWR, LTHR, LAT, LASYM, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF7[0], range(2, 6))
    if LTHR not in range(0,4):
        raise NotImplementedError("MF7 with LTHR=%d" % LTHR)

    line = 1
    forms = []
    if( LTHR in (1, 3) ) :                                           # coherent elastic scattering

        line, temps, energies, dummy, data, e_interp, t_interp = readSTable( line, MF7 )
        temps, energies = map( valuesModule.values, ( temps, energies ) )

        if( e_interp != 'flat' ) : raise ValueError( "Unexpected energy interpolation encountered" )

        axes = axesModule.axes( rank = 3 )
        axes[2] = axesModule.grid( 'temperature', 2, 'K', axesModule.pointsGridToken, values = temps, interpolation = t_interp )
        axes[1] = axesModule.grid( 'energy_in', 1, 'eV', axesModule.pointsGridToken, values=energies, interpolation = e_interp )
        axes[0] = axesModule.axis( 'S_cumulative', 0, 'eV*b' )
        array = arrayModule.full( shape = ( len( temps ), len( energies ) ), data = data )
        Stab = griddedModule.gridded2d( axes, array )

        forms.append( coherentElasticModule.form( info.style, coherentElasticModule.S_table( Stab ) ) )

    if( LTHR in (2, 3) ) :                                      # incoherent elastic
        line, dat  = endfFileToGNDSMisc.getTAB1(line, MF7)
        SB = incoherentElasticModule.characteristicCrossSection( float( dat['C1'] ), 'b' )

        e_interp = dat['interpolationInfo']
        if( len( e_interp ) > 1 ) : raise ValueError( "only one interpolation region allowed for Debye/Waller data" )
        e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(e_interp[0][1])

        t_axes = axesModule.axes( labelsUnits = { 1 : ( 'temperature','K' ), 0 : ( 'DebyeWallerIntegral', '1/eV' ) } )
        if( len( dat['data'] ) == 2 and dat['data'][0] == dat['data'][1] ) : dat['data'] = dat['data'][0:1]

        function1d = XYs1dModule.XYs1d( data = dat['data'], axes = t_axes, interpolation = e_interp )
        DbW = incoherentElasticModule.DebyeWaller( function1d )

        forms.append( incoherentElasticModule.form( info.style, SB, DbW ) )

    if( LTHR == 0 ) :                                           # incoherent inelastic
        line, listy = endfFileToGNDSMisc.getList(line, MF7)
        LLN, NI, NS, b_n = listy['L1'], listy['NPL'], listy['N2'], listy['data']
        if( LLN != 0 ) : raise NotImplementedError( "LLN != 0 not yet handled!" )

            # b_n array contains information about principal / secondary scattering atoms
        numberPerMolecule = int( b_n[5] )
        freeAtomCrossSection = b_n[0] / numberPerMolecule

        neutronMassAMU = info.PoPs['n'].mass.float('amu')
        atom = incoherentInelasticModule.scatteringAtom( label = "0", numberPerMolecule = numberPerMolecule,
                    mass = incoherentInelasticModule.mass( b_n[2] * neutronMassAMU, 'amu' ),
                    freeAtomCrossSection = incoherentInelasticModule.freeAtomCrossSection( freeAtomCrossSection, 'b' ),
                    e_critical = incoherentInelasticModule.e_critical( b_n[1], 'eV' ), e_max = incoherentInelasticModule.e_max( b_n[3], 'eV' ) )
        atoms = [ atom ]

        for index in range( 1, NS + 1 ) :
            functionalForm = { 0.0 : 'SCT', 1.0 : 'free_gas', 2.0: 'diffusive_motion' }[b_n[6*index]]
            numberPerMolecule = int( b_n[6*index+5] )
            freeAtomCrossSection = b_n[6*index+1] / numberPerMolecule
            atoms.append( incoherentInelasticModule.scatteringAtom( label = str( index ), numberPerMolecule = numberPerMolecule,
                    mass = incoherentInelasticModule.mass( b_n[6*index+2] * neutronMassAMU, 'amu' ),
                    freeAtomCrossSection = incoherentInelasticModule.freeAtomCrossSection( freeAtomCrossSection,'b' ), functionalForm = functionalForm ) )

        line, t2header = endfFileToGNDSMisc.getTAB2Header(line, MF7)
        b_interp = t2header['interpolationInfo']
        if len(b_interp) > 1: raise ValueError( "only one interpolation region allowed for S table" )
        b_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(b_interp[0][1])
        n_betas = int( t2header['NZ'] )

        line, temps, alphas, beta, data, a_interp, t_interp = readSTable( line, MF7 )   # Read first beta.
        betas = [ beta ]

        for index in range( n_betas - 1 ) :                     # Read in remaining betas.
            nextLine, temps2, alphas2, beta, data2, a_interp2, t_interp2 = readSTable( line, MF7 )
            warnings = []
            if( temps != temps2 ) :
                warnings.append( "inconsistent temperatures for beta index %d starting on line %d of MF7 MT4" % (index, line) )
            if( alphas != alphas2 ) :
                warnings.append( "inconsistent alpha grids for beta index %d starting on line %d of MF7 MT4" % (index, line) )
            if( a_interp != a_interp2 or t_interp != t_interp2 ) :
                warnings.append( "inconsistent interpolation flags for beta index %d starting on line %d of MF7 MT4" % (index, line) )
            if warnings:
                warnings = ["WARNING: " + w for w in warnings]
                print("\n".join(warnings))
                info.doRaise.extend(warnings)
            betas.append( beta )
            data.extend( data2 )
            line = nextLine

        temps, betas, alphas = map( valuesModule.values, ( temps, betas, alphas ) )

        axes = axesModule.axes( rank = 4 )
        axes[3] = axesModule.grid( 'temperature', 3, 'K', axesModule.pointsGridToken, values = temps, interpolation = t_interp )
        axes[2] = axesModule.grid( 'beta', 2, '', axesModule.pointsGridToken, values = betas, interpolation = b_interp )
        axes[1] = axesModule.grid( 'alpha', 1, '', axesModule.pointsGridToken, values = alphas, interpolation = a_interp )
        axes[0] = axesModule.axis( 'S_alpha_beta', 0, '' )

        arrayShape_orig = ( len( betas ), len( temps ), len( alphas ) )
        data = numpy.array( data )
        zeroFraction = sum(data==0) / len(data)
        data = data.reshape( arrayShape_orig )
        data = numpy.transpose( data, axes = ( 1, 0, 2 ) )              # ENDF data are stored as 'beta,T,alpha'. Switch order to 'T,beta,alpha'.
        arrayShape_new = ( len( temps ), len( betas ), len( alphas ) )
        if zeroFraction > 0.1:
            array = arrayModule.flattened.fromNumpyArray(data)
        else:
            array = arrayModule.full( shape = arrayShape_new, data = data.flatten( ) )
        Stab = griddedModule.gridded3d( axes, array )

        S_tables = incoherentInelasticModule.S_alpha_beta( Stab )

        # last read T_eff tables. There must be at least one (for the principal scattering atom), plus more for
        # any other atoms that specify the short-collision-time approximation:
        line, t_eff = readT_effective( line, MF7 )
        atoms[0].T_effective = t_eff
        for atom in atoms[1:] :
            if( atom.functionalForm == 'SCT' ) :
                line, t_eff = readT_effective( line, MF7 )
                atom.T_effective = t_eff

        options = incoherentInelasticModule.options( calculatedAtThermal = LAT != 0, asymmetric = LASYM != 0 )

        forms.append( incoherentInelasticModule.form( info.style, options, S_tables, atoms = atoms ) )

    if( line != len( MF7 ) ) : raise ValueError( "Trailing data left in MT%i MF7!" % MT )

    return( forms )
