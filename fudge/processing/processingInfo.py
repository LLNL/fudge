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

conserveParticle = 'conserveParticle'
conserveEnergy = 'conserveEnergy'
conserveParticleAndEnergy = 'conserveParticleAndEnergy'

from PoPs import IDs as IDsPoPsModule

from fudge.gnd import styles as stylesModule
from fudge.gnd import reactionSuite as reactionSuiteModule

class tempInfo :

    def __init__( self ) :

        self.dict = {}

    def __contains__( self, key ) :

        return( key in self.dict )

    def __delitem__( self, key ) :

        del self.dict[key]

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) : 

        self.dict[key] = value

    def get( self, key, default=None ):

        return self.dict.get( key, default )

    def keys( self ) :

        return( self.dict.keys( ) )

class processInfo2 :

    def __init__( self, style, reactionSuite, flux = None, logFile = None, verbosity = 0, verbosityIndentStep = '  ',
            energyAccuracy = 1e-6, momentumAccuracy = 1e-3 ) :

        if( not( isinstance( style, stylesModule.style ) ) ) : raise TypeError( 'invalid style' )
        self.__style = style

        if( not( isinstance( reactionSuite, reactionSuiteModule.reactionSuite ) ) ) : raise TypeError( 'invalid reactionSuite' )
        self.__reactionSuite = reactionSuite

        self.__flux = flux
        self.__logFile = logFile
        self.__verbosity = int( verbosity )

        if( not( isinstance( verbosityIndentStep, str ) ) ) : raise TypeError( 'invalid verbosityIndentStep' )
        if( len( verbosityIndentStep ) != verbosityIndentStep.count( ' ' ) ) :
            raise ValueError( 'verbosityIndentStep must only contain the space character' )
        self.__verbosityIndentStep = verbosityIndentStep

        self.__energyAccuracy = energyAccuracy
        self.__momentumAccuracy = momentumAccuracy

    @property
    def energyAccuracy( self ) :

        return( self.__energyAccuracy )

    @property
    def momentumAccuracy( self ) :

        return( self.__momentumAccuracy )

    @property
    def reactionSuite( self ) :

        return( self.__reactionSuite )

    @property
    def style( self ) :

        return( self.__style )

    @property
    def verbosityIndentStep( self ) :

        return( self.__verbosityIndentStep )

    @property
    def flux( self ) :

        return( self.__flux )

    @property
    def logFile( self ) :

        return( self.__logFile )

    @property
    def verbosity( self ) :

        return( self.__verbosity )

    def checkStyle( self, _class ) :

        if( not( issubclass( _class, stylesModule.style ) ) ) : raise TypeError( 'class not a sub-class of style' )
        if( isinstance( self.style, _class ) ) : return
        raise TypeError( 'style %s not of style %s' % ( _class.moniker, self.style.moniker ) )

class processInfo :

    def __init__( self, target, particles = None, flux = None, logFile = None, verbosity = 0 ) :

        self.target = target
        self.dict = {}
        self['particles'] = {}
        if( not( particles is None ) ) :
            for name in particles : self.addParticle( particles[name] )
        if( not flux is None ) : self['flux'] = flux
        self['verbosityIndent'] = ''
        self['verbosity'] = verbosity
        self['logFile'] = logFile

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) :

        self.dict[key] = value

    def addParticle( self, particle ) :

        self['particles'][particle.id]  = particle

    def getParticleGroups( self, name ) :

        return( self['particles'][name].groups )

    def getParticle_lMax( self, name ) :

        return( self['particles'][name].lMax )

    def getProjectileName( self ) :

        return( self.target.projectile.name )

    def getTargetName( self ) :

        return( self.target.target.name )

    def isProcessParticle( self, name ) :

        return( name in self['particles'] )

    def process( self, verbosity = None, verbosityIndent = None ) :

        verbositySave = self['verbosity']
        if( not ( verbosity is None ) ) : self['verbosity'] = verbosity
        if( verbosityIndent is None ) : verbosityIndent = self['verbosityIndent']
        doRaise = False
        try :
            self.target.process( self, verbosityIndent )
        except :
            doRaise = True
        self['verbosity'] = verbositySave
        if( doRaise ) : raise

class processInfoParticle :

    def __init__( self, name, groups, lMax, conservationFlag = conserveParticle ) :

        self.name = name
        self.groups = groups
        self.lMax = lMax
        self.conservationFlag = conservationFlag

    def __repr__( self ) :

        s = '\n%s %s %s\n' % ( self.name, self.lMax, self.conservationFlag )
        s += `self.groups`
        return( s )

class processInfoLLNL( processInfo ) :

    def __init__( self, target, groups = None, flux = None, LLNL_Pn = True, lMax = 3, logFile = None, verbosity = 0 ) :

        particles = {}
        for particle in groups :
            lMax_ = 0
            if( particle == target.projectile.name ) : lMax_ = lMax
            conservationFlag = conserveParticleAndEnergy
            if( particle == 'n' ) : conservationFlag = conserveParticle
            if( particle == IDsPoPsModule.photon ) : conservationFlag = conserveEnergy
            particles[particle] = processInfoParticle( particle, groups[particle], lMax_, conservationFlag )
        processInfo.__init__( self, target, particles, flux = flux, logFile = logFile, verbosity = verbosity )
        self['workDir'] = 'xndfgen.work'
        styles = []
        if( LLNL_Pn ) : styles.append( 'LLNL_Pn' )
        self['styles'] = styles
