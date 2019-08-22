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

from . import base, productBuilder
from fudge import particles
from fudge.gnd import xParticle, xParticleList
from pqu import PQU # for physical quantities

# particleLibrary keeps track of all particles that have been defined,
# so that an evaluation contains only one copy of each particle.
class particleLibrary:
    def __init__(self, database=None):
        '''Creates particle library to keep track of all particles used in the evaluation.
        Once an external particle database is ready, this library could point to that database.'''
        if database:
            raise Exception("No support yet for using external particle database!")
        self.particleList = xParticleList.xParticleList()

    def getParticle(self, name, mass=None, levelIndex=None, levelEnergy=None):
        """Each particle should be defined only once. This checks whether the particle is already defined,
        and either returns it or defines and returns the result.

        Currently we must manually add the energy level for an excited state.
        for example:
            >prod = getParticle('Fe56')
            >lev1 = getParticle('Fe56_e1', levelIndex=1, levelEnergy=PQU('...') )
            >lev2 = getParticle('Fe56', levelIndex=2, levelEnergy=PQU('...') )

        In the future we'll get that information from the particle database instead."""

        if type(name) != type(''):
            raise Exception("Particle name must be a string, not %s" % type(name))

        if levelIndex is not None and '_' not in name:
            if type(levelIndex) == type(0): name += '_e%i' % levelIndex
            else: name += '_%s' % levelIndex

        if type(mass) is str: mass = PQU(mass)

        if self.particleList.hasParticle( name ): return self.particleList.getParticle( name )

        if name == "gamma":
            particle = xParticle.xParticle( name, genre="photon", mass = PQU(0,"amu") )
        elif levelIndex is not None:
            particle = xParticle.nuclearLevel( name, energy = PQU(0,'eV'), label = levelIndex )
            if levelEnergy is not None:
                if isinstance( levelEnergy, PQU ): particle.energy = levelEnergy
                else: raise Exception( "Level energy must be given as a physical quantity" )
            else: particle.needsEnergy = True
        else:
            particle = xParticle.xParticle( name, xParticle.particleType_Nuclear, mass = PQU(0,"amu") )
            if mass is not None: particle.setMass( mass )
            else: particle.needsMass = True
        self.particleList.addParticle( particle )
        return particle

    def newProduct(self, name, mass=None, levelIndex=None, levelEnergy=None, **kwargs):
        particle = self.getParticle( name, mass=mass, levelIndex=levelIndex, levelEnergy=levelEnergy )
        return productBuilder.productBuilder( particle, **kwargs )

    def setMass(self, name, mass, massUnit):
        """
        Set the mass for particle named <name>.
        """
        particle = self.particleList.getParticle( name )
        if not hasattr( particle, 'needsMass' ) and particle.getMass(massUnit) != mass:
            raise Exception( "Encountered multiple values for %s mass: %s vs %s %s" % (particle,
                particle.getMass(massUnit), mass, massUnit) )
        particle.setMass( PQU( mass, massUnit ) )
        if hasattr( particle, 'needsMass' ): del particle.needsMass

    def setExcitationEnergy(self, name, energy, energyUnit):
        """
        Set the excitation for level named <name>
        """
        level = self.particleList.getParticle( name )
        if not hasattr( level, 'needsEnergy' ) and level.getLevelAsFloat(energyUnit) != energy:
            raise Exception( "Encountered multiple values for %s energy: %s vs %s %s" % (name,
                level.getLevelAsFloat(energyUnit), energy, energyUnit) )
        level.energy = PQU(energy, energyUnit)
        if hasattr( level, 'needsEnergy' ): del level.needsEnergy
    
    def setGammaBranchings( self, name, gammas=[] ):
        '''
        Sets the gamma branchings out of a level named <name>
        
        gammas must be a list of xParticle.nuclearLevelGamma instances
        '''
        level = self.particleList.getParticle( name )
        for gamma in gammas: level.addGamma( gamma )

    def checkForUnknownMassOrExcitation(self):
        unknownMass, unknownExcitation = [], []
        for particle in self.particleList.values():
            if hasattr(particle,'needsMass'): unknownMass.append( particle )
            for level in particle.levels.values():
                if hasattr(level,'needsEnergy'): unknownExcitation.append( level )
        return unknownMass, unknownExcitation

