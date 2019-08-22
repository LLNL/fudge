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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
        self.rs.setAttribute( 'temperature', PQU.PQU(0,'K') )

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
                    [p.name for p in badMass] )
        if badEnergy:
            raise Exception( "Energies of the following levels are still undefined: %s" %
                    [p.name for p in badEnergy] )
        self.rs.particles = self.particleLibrary.particleList

        for index, reaction in enumerate( self.rs.reactions ):
            reaction.setLabel( index )
            mt = self.__calculateMT( reaction )
            reaction.ENDF_MT = mt
        
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
        excitation = 0
        for prod in reaction.outputChannel:
            particle = prod.name
            if particle == 'gamma':
                continue
            elif particle in pdict :
                pdict[ particle ] += prod.multiplicity.getConstant()
            else:               # Should be the residual. Check if it's in an excited state:
                if( hasattr( prod.particle, 'getLevelIndex' ) ) : excitation = prod.particle.getLevelIndex()

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

