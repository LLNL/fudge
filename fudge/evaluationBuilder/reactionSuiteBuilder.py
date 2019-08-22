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

from . import reactionBuilder, productBuilder
import fudge
from fudge import gnd
from pqu import PQU

class reactionSuiteBuilder:
    def __init__(self, particleLibrary, projectile, target, library = "ENDF/B", version = None,
            documentation = None, transportables = ['n']):
        if version is None:
            import datetime
            version = datetime.datetime.today().strftime("%Y")
        self.rs = gnd.reactionSuite.reactionSuite( projectile, target,
                style = gnd.miscellaneous.style("evaluated", {'library':library,'version':version} ) )
        self.rs.addDocumentation( gnd.documentation.documentation(library,
                documentation or "skeleton for documentation" ) )
        self.rs.setAttribute( 'temperature', PQU(0,'K') )

        self.rs.projectile.attributes['transportable'] = True
        self.particleLibrary = particleLibrary
        self.transportables = transportables

        self.units = None

    def setUnits(self, energy, crossSection, mass, time, **kwargs):
        self.units = { 'energy': energy, 'crossSection': crossSection, 'mass': mass, 'time': time }
        self.units.update( kwargs )

    def addReaction( self, reaction ):
        if isinstance( reaction, fudge.gnd.reactions.reaction.reaction ):
            self.rs.addReaction( reaction )
        elif isinstance( reaction, reactionBuilder.reactionBuilder ):
            self.rs.addReaction( reaction.finalize( self.rs, self.particleLibrary ) )
        else:
            raise Exception("addReaction() argument must be gnd.reaction or evaluator.reactionBuilder instance")

    def addSummedReaction( self, reaction ):
        raise NotImplementedError()

    def finalize( self ):
        """
        Fix up some data, and do some testing before returning the reactionSuite instance
        """

        badMass, badEnergy = self.particleLibrary.checkForUnknownMassOrExcitation()
        if badMass:
            raise Exception( "Masses of the following products are still undefined: %s" %
                    [p.getName() for p in badMass] )
        if badEnergy:
            raise Exception( "Energies of the following levels are still undefined: %s" %
                    [p.getName() for p in badEnergy] )
        self.rs.particles = self.particleLibrary.particleList

        for index, reaction in enumerate( self.rs.reactions ):
            reaction.setLabel( index )
            mt = self.__calculateMT( reaction )
            reaction.attributes['ENDF_MT'] = mt
        
        for transportable in self.transportables:
            particle = self.rs.getParticle( transportable )
            particle.attributes['transportable'] = True
        #warnings = self.rs.check()
        return self.rs

    def __calculateMT( self, reaction ):
        """
        For internal use only: get the ENDF MT number from the list of products.
        """
        pdict = dict.fromkeys( ('n','H1','H2','H3','He3','He4'), 0 )
        for prod in reaction.outputChannel:
            particle = prod.getToken()
            if particle == 'gamma': continue
            elif particle in pdict: pdict[ particle ] += prod.multiplicity.getConstant()
            else:
                # should be the residual. Check if it's in an excited state:
                excitation = prod.particle.getLevelIndex()

        def getProductTuple( particleDictionary, excitation ):
            return [particleDictionary[key] for key in ('n','H1','H2','H3','He3','He4')] + [excitation]

        productTuple = getProductTuple( pdict, excitation )
        # some special cases:
        if not any( productTuple ): return 102    # capture
        if ( sorted( [prod.particle for prod in reaction.outputChannel] ) ==
                sorted( [self.rs.projectile, self.rs.target] ) ): return 2    # elastic

        from fudge.legacy.converting.endf_endl import endfMTtoC_ProductLists
        mts = [mt for mt in endfMTtoC_ProductLists if getProductTuple(
            endfMTtoC_ProductLists[ mt ].productCounts, endfMTtoC_ProductLists[mt].residualLevel
            ) == productTuple]

        if len(mts) == 1: return mts[0]
        if not mts: return -1
        raise Exception("Got multiple MTs: %s" % mts)

