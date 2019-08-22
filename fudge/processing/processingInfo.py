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

conserveParticle = 'conserveParticle'
conserveEnergy = 'conserveEnergy'
conserveParticleAndEnergy = 'conserveParticleAndEnergy'

class tempInfo :

    def __init__( self ) :

        self.dict = {}

    def __contains__( self, key ) :

        return( key in self.dict )

    def __getitem__( self, key ) :

        return( self.dict[key] )

    def __setitem__( self, key, value ) : 

        self.dict[key] = value

    def get( self, key, default=None ):

        return self.dict.get( key, default )

    def keys( self ) :

        return( self.dict.keys( ) )

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

        self['particles'][particle.name]  = particle

    def getParticleGroups( self, name ) :

        return( self['particles'][name].groups )

    def getParticle_lMax( self, name ) :

        return( self['particles'][name].lMax )

    def getProjectileName( self ) :

        return( self.target.projectile.getName( ) )

    def getTargetName( self ) :

        return( self.target.target.getName( ) )

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

    def __init__( self, target, groups = None, flux = None, LLNL_MC = True, LLNL_Pn = True, lMax = 3, logFile = None, verbosity = 0 ) :

        particles = {}
        for particle in groups :
            lMax_ = 0
            if( particle == target.projectile.getName( ) ) : lMax_ = lMax
            conservationFlag = conserveParticleAndEnergy
            if( particle == 'n' ) : conservationFlag = conserveParticle
            if( particle == 'gamma' ) : conservationFlag = conserveEnergy
            particles[particle] = processInfoParticle( particle, groups[particle], lMax_, conservationFlag )
        processInfo.__init__( self, target, particles, flux = flux, logFile = logFile, verbosity = verbosity )
        self['workDir'] = 'xndfgen.work'
        styles = []
        if( LLNL_MC ) : styles.append( 'LLNL_MC' )
        if( LLNL_Pn ) : styles.append( 'LLNL_Pn' )
        self['styles'] = styles
