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

import sys
import fractions

from pqu import PQU as PQUModule
from xData import standards as standardsModule

from PoPs import IDs as IDsPoPsModule
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
from PoPs.families import nucleus as nucleusModule
from PoPs.families import unorthodox as unorthodoxModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule

from fudge.core.utilities import brb
from fudge.core.math import fudgemath

from fudge.gnds import product as productModule
from fudge.gnds.productData import multiplicity as multiplicityModule
from fudge.gnds.channelData import Q as QModule

from fudge.gnds.productData.distributions import unspecified as unspecifiedModule

import massTracker as massTrackerModule

def ZAToName( ZA ) :

    if( ZA == 99120 ) : return( 'FissionProductENDL99120' )
    if( ZA == 99125 ) : return( 'FissionProductENDL99125' )
    return( chemicalElementMiscPoPsModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementMiscPoPsModule.symbolFromZ[ZA/1000], ZA % 1000 ) )

def addParticleData( particle, info, massValue, spinValue, parityValue, chargeValue, halflifeValue ) :

    if( massValue is not None ) :
        mass = massModule.double( info.PoPsLabel, massValue, quantityModule.stringToPhysicalUnit( 'amu' ) )
        particle.mass.add( mass )

    if( spinValue is not None ) :
        spin = spinModule.fraction( info.PoPsLabel, fractions.Fraction( spinValue ), spinModule.baseUnit )
        particle.spin.add( spin )

    if( parityValue is not None ) :
        parity = parityModule.integer( info.PoPsLabel, parityValue, parityModule.baseUnit )
        particle.parity.add( parity )

    if( chargeValue is not None ) :
        charge = chargeModule.integer( info.PoPsLabel, chargeValue, chargeModule.baseUnit )
        particle.charge.add( charge )

    if( halflifeValue is not None ) :
        if( isinstance( halflifeValue, str ) ) :
            halflife = halflifeModule.string( info.PoPsLabel, halflifeValue, halflifeModule.baseUnit )
        else :
            halflife = halflifeModule.double( info.PoPsLabel, halflifeValue, halflifeModule.baseUnit )
        particle.halflife.add( halflife )

def addNucleusInfoForLightNuclei( ZA, nucleus, info ) :

    spinValue, parityValue, halflifeValue = None, None, None
    if( ZA == 1001 ) :
        spinValue, parityValue, halflifeValue = '1/2', 1, 'stable'
    elif( ZA == 1002 ) :
        spinValue, parityValue, halflifeValue = '1', 1, 'stable'
    elif( ZA == 1003 ) :
        spinValue, parityValue, halflifeValue = '1/2', 1, 3.88789e+8
    elif( ZA == 2003 ) :
        spinValue, parityValue, halflifeValue = '1/2', 1, 'stable'
    elif( ZA == 2004 ) :
        spinValue, parityValue, halflifeValue = '0', 1, 'stable'
    else :
        return

    addParticleData( nucleus, info, None, spinValue, parityValue, None, halflifeValue )

def getPoPsParticle( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :

    if( levelIndex is not None ) :
        levelIndex = int( levelIndex )
        if( levelIndex < 0 ) : raise Exception( 'levelIndex = %d must be >= 0' % levelIndex )

    particle = None

    if( ZA in[ 0, 17 ] ) : ZA = 7
    if( name is not None ) :
        pass
    elif( ZA == 1 ) :
        particle = baryonModule.particle( 'n' )
        addParticleData( particle, info, 1.00866491574, "1/2", 1, 0, 881.5 )
        name = particle.id
    elif( ZA == 7 ) :
        particle = gaugeBosonModule.particle( IDsPoPsModule.photon )
        addParticleData( particle, info, 0.0, "1", 1, 0, halflifeModule.STABLE )
        name = particle.id
    elif( ZA == 8 ) :
        particle = leptonModule.particle( 'e+', generation = 'electronic' )
# BRB need to make it 'e-_anti' and alias to 'e+'.
#        particle = leptonModule.particle( 'e-_anti', generation = 'electronic' )
        addParticleData( particle, info, 5.485799090e-4, "1/2", -1, 1, halflifeModule.STABLE )
        name = particle.id
    elif( ZA == 9 ) :
        particle = leptonModule.particle( 'e-', generation = 'electronic' )
        addParticleData( particle, info, 5.485799090e-4, "1/2", 1, -1, halflifeModule.STABLE )
        name = particle.id
    elif( ZA in [ 99120, 99125 ] ) :
        name = ZAToName( ZA )
        particle = unorthodoxModule.particle( name )
        mass = massModule.double( info.PoPsLabel, 117.5, quantityModule.stringToPhysicalUnit( 'amu' ) )
        particle.mass.add( mass )
    else :
        name = ZAToName( ZA )
        if( levelIndex is None ) :
            levelIndex = 0
            level = 0.
        if( level < 0. ) : raise Exception( 'Negative value = %s for continuum is not allowed' % level )
        level = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( level ), levelUnit )
        name = chemicalElementMiscPoPsModule.nuclideIDFromIsotopeSymbolAndIndex( name, levelIndex )

    if( ( particle is None ) and ( name not in info.PoPs ) ) :

        if( level is not None ) :                                   # Add a nuclide.
            baseName = name.split('_')[0]                           # Always need to add ground state before excited level.
            if( 'FissionProduct' in name ) : baseName = name.split('_')[0]

            particle = nuclideModule.particle( name )
            charge = chargeModule.integer( info.PoPsLabel, 0, chargeModule.baseUnit )
            particle.charge.add( charge )

            energy = nuclearEnergyLevelModule.double( info.PoPsLabel, float( level ), level.unit )
            particle.nucleus.energy.add( energy )

            charge = chargeModule.integer( info.PoPsLabel, ZA / 1000, chargeModule.baseUnit )
            particle.nucleus.charge.add( charge )
            addNucleusInfoForLightNuclei( ZA, particle.nucleus, info )
        else :
            if( particle is None ) : raise Exception( 'FIX ME' )
    else :
        if( particle is None ) : particle = info.PoPs[name]

    if( name not in info.PoPs ) : info.PoPs.add( particle )
    return( particle )

def getTypeName( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :
    """
    Returns the name for this ZA and level if present. Returned name is of the form Am242, Am242_m1, Am242_e2.
    levelIndex must be None or an integer > 0.
    """

    return( getPoPsParticle( info, ZA, name, levelIndex, level, levelUnit ) )

def newGNDSParticle( info, particle, crossSection, attributes = None, multiplicity = 1, outputChannel = None ) :
    """
    The arguments to this function must match those of product.product excluding self and label.
    """

    if attributes is None: attributes = {}
    name = particle
# BRBBRB
    if( not( isinstance( particle, str ) ) ) : name = particle.id
    product = productModule.product( id = name, label = name, attributes = attributes, outputChannel = outputChannel )
    if( isinstance( multiplicity, ( int, float ) ) ) :
        axes = multiplicityModule.defaultAxes( crossSection.domainUnit )
        multiplicity = multiplicityModule.constant1d( multiplicity, crossSection.domainMin,
                crossSection.domainMax, axes = axes, label = info.style )
    product.multiplicity.add( multiplicity )
    return( product )

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
                    raise Exception( "undefinedLevelInfo['count'] > 1 for ZA = %s" % ZA )
            else:
                levelIndex, level = undefinedLevelInfo['levelIndex'], undefinedLevelInfo['level']
    return( getTypeNameGamma( info, ZA, level = level, levelIndex = levelIndex ) )

class infos :

    def __init__( self, style ) :

        self.style = style
        self.massTracker = massTrackerModule.massTracker()
        self.__logs = None

    @property
    def logs(self): return self.__logs

    @logs.setter
    def logs(self, logs): self.__logs = logs

    def addMassAWR( self, ZA, AWR, asTarget=True ):
        warning = self.massTracker.addMassAWR( ZA, AWR, asTarget=asTarget )
        if warning: self.logs.write( warning )

def returnConstantQ( style, Q, crossSection ) :

    axes = QModule.defaultAxes( crossSection.domainUnit )
    return( QModule.constant1d( Q, crossSection.domainMin, crossSection.domainMax, axes = axes, label = style ) )

def addUnspecifiedDistributions( info, outputChannel, frame = None ) :
    """
    For products with no distribution, assign an unspecified distribution.
    Interim products before breakup get an 'implicitProduct' conversion flag (so they aren't written back to ENDF-6),
    as are heavy residuals when the first product is a neutron.
    Others are explicit, to ensure that the product (plus mass) appears in resulting ENDF file.
    """

    if( outputChannel is None ) : return
    firstProductId = outputChannel[0].pid
    for product in outputChannel :
        if( len( product.distribution ) == 0 ) :
            product.distribution.add(unspecifiedModule.form(info.style, productFrame=frame))
            if product.outputChannel or firstProductId=='n':
                info.ENDFconversionFlags.add( product, 'implicitProduct' )
        else :
            if( frame is None ) : frame = product.distribution[info.style].productFrame
        addUnspecifiedDistributions( info, product.outputChannel, frame = standardsModule.frames.labToken )

