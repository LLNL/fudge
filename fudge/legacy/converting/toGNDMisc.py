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

from pqu import PQU
from xData import standards as standardsModule

from fudge import gnd
from fudge.structure import masses as massModule
from fudge.legacy.endl import endl2
from fudge.core.utilities import brb
from fudge.core.math import fudgemath

from fudge.gnd.productData.distributions import unspecified as unspecifiedModule

def getTypeName( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :
    """Returns the name for this ZA and level if present. Returned name is of the form Am_242, Am_242_m1, Am_242_e2.
    levelIndex must be None, 'c', 's' or an integer > 0."""

    if( levelIndex not in [ 'c', 's', None ] ) :
        levelIndex = int( levelIndex )
        if( levelIndex < 0 ) : raise Exception( 'levelIndex = %d must be >= 0' % levelIndex )
    if( ZA == 17 ) : ZA = 7
    particleQualifiers = {}
    if( name is not None ) :
        pass
    elif( ZA == 7 ) :
        name = 'gamma'
    elif( ZA == 8 ) :
        name = 'e+'
    elif( ZA == 9 ) :
        name = 'e-'
    else :
        name = endl2.endlToGNDName( ZA )
        if( levelIndex not in [ None, 0 ] ) :
            if( level < 0. ) : raise Exception( 'Negative value = %s for continuum is not allowed' % level )
            if( levelIndex in [ 'c', 's' ] ) :        # to continuum or sum of all levels.
                particleQualifiers['level'] = gnd.xParticle.undefinedLevel( PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( level ), levelUnit ) )
                name += "_%s" % levelIndex
            else :
                particleQualifiers['level'] = PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( level ), levelUnit )
                if( levelIndex != 0 ) : name += '_e%d' % levelIndex

    if( not( info.particleList.hasParticle( name ) ) ) :
        if( ZA in info.ZAMasses ) :
            if( ( info.ZAMasses[ZA] is not None ) and ( info.ZAMasses[ZA] < 0 ) ) : info.ZAMasses[ZA] *= -1
            mass = info.ZAMasses[ZA]
        else :
            mass = info.masses.getMassFromZA( ZA )
            if( mass is None ) :                                        # If mass is not found, let's attempt to make something up that is close.
                A = ZA % 1000
                for idx in xrange( 1, A ) :
                    massLower = info.masses.getMassFromZA( ZA - idx )
                    massHigher = info.masses.getMassFromZA( ZA + idx )
                    if( massLower is not None ) :
                        mass = massLower + idx * info.masses.getMassFromZA( 1 )
                        break
                    elif( massHigher is not None ) :
                        mass = massHigher - idx * info.masses.getMassFromZA( 1 )
                        break
                if( mass is None ) : raise Exception( 'Could not obtain mass for ZA = %s' % ZA )
                mass = round(mass,1)    # since we're extrapolating
        if( mass is None ) :
            massUnit = None
        else :
            massUnit = PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( mass ), 'amu' )
        if( name in info.transportables ) : particleQualifiers['transportable'] = 'true'
        if( ( levelIndex is not None ) and ( levelIndex != 0 ) ) : particleQualifiers['levelIndex'] = levelIndex
        if( name in [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4', 'gamma' ] ) : particleQualifiers['transportable'] = 'true'

        if( particleQualifiers.get( 'level' ) is not None ) :                   # Add new particle or level.
                # careful, always need to add ground state before excited level:
            baseName = name.split('_')[0]
            if( 'natural' in name ) :
                baseName = '_'.join( name.split('_')[:2] )
            elif( 'FissionProduct' in name ) :
                baseName = name.split('_')[0]
            if( baseName not in info.particleList ) :
                info.particleList.addParticle( gnd.xParticle.isotope( baseName, mass = massUnit, attributes={} ) )

                # now ground state is present, add excited state
            if( 'levelIndex' in particleQualifiers ) :
                index = particleQualifiers.pop('levelIndex')
            elif( particleQualifiers['level'].getValueAs( 'MeV' ) == '0.' ) :
                index = 0
            else :
                raise Exception( 'Index could not be determined for name = "%s", level = "%s", levelIndex = "%s"' % ( name, level, levelIndex ) )
            energy = particleQualifiers.pop( 'level' )
            particleOrLevel = gnd.xParticle.nuclearLevel( name, energy, index, attributes = particleQualifiers )
        else:
            if( name == 'e+' ) :
                particleOrLevel = gnd.xParticle.lepton( name, 'electronic', massUnit, +1, anti = True, attributes = particleQualifiers )
            elif( name == 'e-' ) :
                particleOrLevel = gnd.xParticle.lepton( name, 'electronic', massUnit, -1, attributes = particleQualifiers )
            elif( name == 'gamma' ) :
                particleOrLevel = gnd.xParticle.photon( attributes = particleQualifiers )
            elif( 'FissionProduct' in name ) :
                particleOrLevel = gnd.xParticle.FissionProduct( name, mass = massUnit, attributes = particleQualifiers )
            elif( 'TNSL' in name ) :
                particleOrLevel = gnd.xParticle.thermalNeutronScatteringLawIsotope( name, mass = massUnit, attributes = particleQualifiers )
            else :
                particleOrLevel = gnd.xParticle.isotope( name, mass = massUnit, attributes = particleQualifiers )
        info.particleList.addParticle( particleOrLevel )

    return( info.particleList.getParticle( name ) )

def newGNDParticle( info, particle, attributes = None, multiplicity = 1, outputChannel = None ) :
    """The arguments to this function must match those of product.product excluding self and label."""

    if attributes is None: attributes = {}
    product = gnd.product.product( particle, label = particle.name, attributes = attributes, outputChannel = outputChannel )
    if( isinstance( multiplicity, ( int, float ) ) ) : multiplicity = gnd.productData.multiplicity.constant( info.style, multiplicity )
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
        if( undefinedLevelInfo['ZA'] is not None ) :
            if( ( undefinedLevelInfo['ZA'] == ZA ) and ( ZA not in [ 1001, 1002, 1003, 2003, 2004 ] ) ) :
                undefinedLevelInfo['count'] += 1
                if( undefinedLevelInfo['count'] > 1 ) : raise Exception( "undefinedLevelInfo['count'] > 1 for ZA = %s" % ZA )
                levelIndex, level = undefinedLevelInfo['levelIndex'], undefinedLevelInfo['level']
    return( getTypeNameGamma( info, ZA, level = level, levelIndex = levelIndex ) )

class infos :

    def __init__( self, style, xenslIsotopes = None, transportables = None ) :

        self.style = style
        self.masses = massModule
        self.xenslIsotopes = xenslIsotopes
        if( xenslIsotopes is None ) : self.xenslIsotopes = {}
        self.transportables = transportables or []
        self.particleList = gnd.xParticleList.xParticleList( )
        self.ZAMasses = {}

    def setReactionSuite( self, reactionSuite )  :

        self.reactionSuite = reactionSuite

def returnConstantQ( style, Q, unit = 'eV' ) :

    return( gnd.channelData.Q.constant( style, PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( Q ), unit ) ) )

def addUnspecifiedDistributions( info, outputChannel, frame = None ) :

    if( outputChannel is None ) : return
    for product in outputChannel :
        if( len( product.distribution ) == 0 ) :
            product.distribution.add( unspecifiedModule.form( info.style, productFrame = frame ) )
        else :
            if( frame is None ) : frame = product.distribution[info.style].productFrame
        addUnspecifiedDistributions( info, product.outputChannel, frame = standardsModule.frames.labToken )
