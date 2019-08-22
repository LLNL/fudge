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
            >lev1 = getParticle('Fe56_e1', levelIndex=1, levelEnergy=PQU.PQU('...') )
            >lev2 = getParticle('Fe56', levelIndex=2, levelEnergy=PQU.PQU('...') )

        In the future we'll get that information from the particle database instead."""

        if type(name) != type(''):
            raise Exception("Particle name must be a string, not %s" % type(name))

        if levelIndex is not None and '_' not in name:
            if type(levelIndex) == type(0): name += '_e%i' % levelIndex
            else: name += '_%s' % levelIndex

        if type(mass) is str: mass = PQU.PQU(mass)

        if self.particleList.hasParticle( name ): return self.particleList.getParticle( name )

        if name == "gamma":
            particle = xParticle.photon( )
        elif levelIndex is not None:
            particle = xParticle.nuclearLevel( name, energy = PQU(0,'eV'), label = levelIndex )
            if levelEnergy is not None:
                if isinstance( levelEnergy, PQU.PQU ): particle.energy = levelEnergy
                else: raise Exception( "Level energy must be given as a physical quantity" )
            else: particle.needsEnergy = True
        else:
            particle = xParticle.isotope( name, mass = PQU.PQU(0,"amu") )
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
# BRB, help
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

    def __contains__(self, name): return self.particleList.hasParticle( name )
