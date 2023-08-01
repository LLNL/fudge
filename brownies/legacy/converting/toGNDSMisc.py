# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import fractions

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import alias as PoPsAliasModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nuclide as nuclideModule
from PoPs.families import unorthodox as unorthodoxModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule

from fudge.core.utilities import brb
from fudge.core.math import fudgemath

from fudge import product as productModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.outputChannelData import Q as QModule

from fudge.productData.distributions import unspecified as unspecifiedModule

from . import massTracker as massTrackerModule

def ZAToName( ZA ) :

    if( ZA == 99120 ) : return( 'FissionProductENDL99120' )
    if( ZA == 99125 ) : return( 'FissionProductENDL99125' )
    return( chemicalElementMiscPoPsModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementMiscPoPsModule.symbolFromZ[ZA//1000], ZA % 1000 ) )

def addParticleData( particle, info, massValue, spinValue, parityValue, chargeValue, halflifeValue ) :

    if( massValue is not None ) :
        mass = massModule.Double( info.PoPsLabel, massValue, quantityModule.stringToPhysicalUnit( 'amu' ) )
        particle.mass.add( mass )

    if( spinValue is not None ) :
        spin = spinModule.Fraction( info.PoPsLabel, fractions.Fraction( spinValue ), spinModule.baseUnit )
        particle.spin.add( spin )

    if( parityValue is not None ) :
        parity = parityModule.Integer( info.PoPsLabel, parityValue, parityModule.baseUnit )
        particle.parity.add( parity )

    if( chargeValue is not None ) :
        charge = chargeModule.Integer( info.PoPsLabel, chargeValue, chargeModule.baseUnit )
        particle.charge.add( charge )

    if( halflifeValue is not None ) :
        if( isinstance( halflifeValue, str ) ) :
            halflife = halflifeModule.String( info.PoPsLabel, halflifeValue, halflifeModule.baseUnit )
        else :
            halflife = halflifeModule.Double( info.PoPsLabel, halflifeValue, halflifeModule.baseUnit )
        particle.halflife.add( halflife )

def addNucleusInfoForLightNuclei( ZA, nucleus, info ) :

    spinValue, parityValue, halflifeValue = None, None, None
    if( ZA == 1001 ) :
        massValue, spinValue, parityValue, chargeValue, halflifeValue = 1.00727646662, '1/2', 1, 1, 'stable'
    elif( ZA == 1002 ) :
        massValue, spinValue, parityValue, chargeValue, halflifeValue = 2.01355321275, '1', 1, 1, 'stable'
    elif( ZA == 1003 ) :
        massValue, spinValue, parityValue, chargeValue, halflifeValue = 3.01550071621, '1/2', 1, 1, 3.88789e+8
    elif( ZA == 2003 ) :
        massValue, spinValue, parityValue, chargeValue, halflifeValue = 3.01493224717, '1/2', 1, 2, 'stable'
    elif( ZA == 2004 ) :
        massValue, spinValue, parityValue, chargeValue, halflifeValue = 4.0015061791, '0', 1, 2, 'stable'
    else :
        return

    addParticleData( nucleus, info, massValue, spinValue, parityValue, chargeValue, halflifeValue )

def getPoPsParticle( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :

    if( levelIndex is not None ) :
        levelIndex = int( levelIndex )
        if( levelIndex < 0 ) : raise Exception( 'levelIndex = %d must be >= 0' % levelIndex )

    particle = None

    if( ZA in [ 0, 17 ] ) : ZA = 7
    if( name is not None ) :
        pass
    elif( ZA == 1 ) :
        particle = baryonModule.Particle(IDsPoPsModule.neutron)
        addParticleData( particle, info, 1.00866491574, "1/2", 1, 0, 881.5 )
        name = particle.id
    elif( ZA == 7 ) :
        particle = gaugeBosonModule.Particle( IDsPoPsModule.photon )
        addParticleData( particle, info, 0.0, "1", 1, 0, halflifeModule.STABLE )
        name = particle.id
    elif( ZA == 8 ) :
        particle = leptonModule.Particle( 'e-_anti', generation = leptonModule.Generation.electronic )
        addParticleData( particle, info, 5.485799090e-4, "1/2", -1, 1, halflifeModule.STABLE )
        name = particle.id
    elif( ZA == 9 ) :
        particle = leptonModule.Particle( 'e-', generation = leptonModule.Generation.electronic )
        addParticleData( particle, info, 5.485799090e-4, "1/2", 1, -1, halflifeModule.STABLE )
        name = particle.id
    elif( ZA in [ 99120, 99125 ] ) :
        name = ZAToName( ZA )
        particle = unorthodoxModule.Particle( name )
        mass = massModule.Double( info.PoPsLabel, 117.5, quantityModule.stringToPhysicalUnit( 'amu' ) )
        particle.mass.add( mass )
        charge = chargeModule.Integer( info.PoPsLabel, 46, chargeModule.baseUnit )
        particle.charge.add( charge )
    else :
        name = ZAToName( ZA )
        if( levelIndex is None ) :
            levelIndex = 0
            level = 0.
        if( level < 0. ) : raise Exception( 'Negative value = %s for continuum is not allowed' % level )
        level = PQUModule.PQU( PQUModule.PQU_float.surmiseSignificantDigits( level ), levelUnit )
        name = chemicalElementMiscPoPsModule.nuclideIDFromIsotopeSymbolAndIndex( name, levelIndex )

    if( ( particle is None ) and ( name not in info.PoPs ) ) :

        if( level is not None ) :                                   # Add a nuclide.
            baseName = name.split('_')[0]                           # Always need to add ground state before excited level.
            if( 'FissionProduct' in name ) : baseName = name.split('_')[0]

            particle = nuclideModule.Particle( name )
            charge = chargeModule.Integer( info.PoPsLabel, 0, chargeModule.baseUnit )
            particle.charge.add( charge )

            energy = nuclearEnergyLevelModule.Double( info.PoPsLabel, float( level ), level.unit )
            particle.nucleus.energy.add( energy )

            charge = chargeModule.Integer( info.PoPsLabel, ZA // 1000, chargeModule.baseUnit )
            particle.nucleus.charge.add( charge )
            addNucleusInfoForLightNuclei( ZA, particle.nucleus, info )
        else :
            if( particle is None ) : raise Exception( 'FIX ME' )
    else :
        if( particle is None ) : particle = info.PoPs[name]

    if name not in info.PoPs:
        info.PoPs.add(particle)
        if name != IDsPoPsModule.photon:
            familiarID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name, specialNuclearParticleIDPoPsModule.Mode.familiar)
            if familiarID != name:
                nucleusID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name, specialNuclearParticleIDPoPsModule.Mode.nucleus)
                if familiarID not in info.PoPs:
                    if familiarID == IDsPoPsModule.proton:
                        proton = baryonModule.Particle(IDsPoPsModule.proton)
                        addNucleusInfoForLightNuclei(1001, proton, info)
                        info.PoPs.add(proton)
                    else:
                        info.PoPs.add(PoPsAliasModule.Alias(familiarID, nucleusID))

    return particle

def getTypeName( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :
    """
    Returns the name for this ZA and level if present. Returned name is of the form Am242, Am242_m1, Am242_e2.
    levelIndex must be None or an integer > 0.
    """

    return( getPoPsParticle( info, ZA, name, levelIndex, level, levelUnit ) )

def newGNDSParticle(info, particle, crossSection, multiplicity=1, outputChannel=None):
    """
    The arguments to this function must match those of product.product excluding self and label.
    """

    name = particle
    if not isinstance(particle, str):
        name = particle.id
    if name in info.PoPs.aliases:
        name = info.PoPs.aliases[name].pid
    if name != IDsPoPsModule.photon:
        name = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name, info.specialNuclearParticleID)

    product = productModule.Product(name, label=name, outputChannel=outputChannel)

    if isinstance(multiplicity, XYs1dModule.XYs1d):
        if multiplicity.rangeMin == multiplicity.rangeMax == 1.0:
            multiplicity = multiplicity.rangeMin
    if isinstance(multiplicity, (int, float)):
        axes = multiplicityModule.defaultAxes( crossSection.domainUnit )
        multiplicity = multiplicityModule.Constant1d(multiplicity, crossSection.domainMin, crossSection.domainMax, axes=axes, label=info.style)
    product.multiplicity.add(multiplicity)

    return product

def getTypeNameGamma( info, ZA, level = None, levelIndex = None ) :

    if( level is not None ) :
        if( fudgemath.isNumber( level ) ) :
            if( level < 0 ) :
                if( ( level > -100 ) and ( levelIndex == 0 ) ) :
                    level = 0
                else :
                    message = 'Negative excitation level = %s for ZA = %s and levelIndex = %s is not allowed' % ( level, ZA, levelIndex )
                    sys.stderr.write( message+'\n' )
                    info.doRaise.append( message )
                    level = 0
        elif( type( level ) != type( '' ) ) :
            raise Exception( 'for ZA %d, level (type ="%s") must be a number or a string' % ( ZA, brb.getType( level ) ) )
    if( ZA == 0 ) : ZA = 7                             # Special case for gammas which are yo = 7 for ENDL.
    p = getTypeName( info, ZA, level = level, levelIndex = levelIndex, levelUnit = 'eV' )
    return( p )

def getTypeNameENDF( info, ZA, undefinedLevelInfo ) :

    levelIndex, level = None, 0.
    if( undefinedLevelInfo is not None ) :
        if( undefinedLevelInfo['ZA'] == ZA ) :
            if undefinedLevelInfo['level'] > 0. and ZA in [ 1001, 1002, 1003, 2003, 2004 ] :
                info.logs.write("  Excited state of ZA %d encountered" % ZA)
            undefinedLevelInfo['count'] += 1
            if( undefinedLevelInfo['count'] > 1 ) :
                if( undefinedLevelInfo['level'] > 0. ) :
                    raise Exception("undefinedLevelInfo['count'] > 1 for ZA = %s: %s" % (ZA, undefinedLevelInfo))
            else:
                levelIndex, level = undefinedLevelInfo['levelIndex'], undefinedLevelInfo['level']
    return( getTypeNameGamma( info, ZA, level = level, levelIndex = levelIndex ) )

class Infos:

    def __init__( self, formatVersion, style ) :

        self.formatVersion = formatVersion
        self.style = style
        self.massTracker = massTrackerModule.MassTracker()
        self.__logs = None
        self.resonanceReconstructionFailed = ''
        self.ZA_massLineInfo = massTrackerModule.ZA_massLineInfo()

    @property
    def logs(self): return self.__logs

    @logs.setter
    def logs(self, logs): self.__logs = logs

    def addMassAWR( self, ZA, AWR, asTarget=True ):
        warning = self.massTracker.addMassAWR( ZA, AWR, asTarget=asTarget )
        if warning: self.logs.write( warning )

def returnConstantQ( style, Q, crossSection ) :

    axes = QModule.defaultAxes( crossSection.domainUnit )
    return( QModule.Constant1d( Q, crossSection.domainMin, crossSection.domainMax, axes = axes, label = style ) )

def addUnspecifiedDistributions(info, outputChannel, frame=xDataEnumsModule.Frame.none) :
    """
    For products with no distribution, assign an unspecified distribution.
    Interim products before breakup get an 'implicitProduct' conversion flag (so they aren't written back to ENDF-6),
    as are heavy residuals when the first product is a neutron.
    Others are explicit, to ensure that the product (plus mass) appears in resulting ENDF file.
    """

    if outputChannel is None:
        return
    firstProductId = outputChannel[0].pid
    for product in outputChannel:
        if len(product.distribution) == 0:
            product.distribution.add(unspecifiedModule.Form(info.style, productFrame=frame))
            if product.outputChannel or firstProductId == IDsPoPsModule.neutron:
                info.ENDFconversionFlags.add( product, 'implicitProduct' )
        else :
            if frame == xDataEnumsModule.Frame.none:
                frame = product.distribution[info.style].productFrame
        addUnspecifiedDistributions(info, product.outputChannel, frame=xDataEnumsModule.Frame.lab)
