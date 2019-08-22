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

import xData.ancestry as ancestryModule
from fudge.gnd import xParticle as xParticleModule

class xParticleList( dict, ancestryModule.ancestry ) : 

    moniker = 'particles'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

    def __iter__( self ) :

        names = sorted( self.keys( ) )
        for name in names :
            yield self[name]

    def addParticle( self, particle ) :
        """Add either a particle or an excited level to xParticleList """

        if( isinstance( particle, ( xParticleModule.isotope, xParticleModule.lepton, xParticleModule.photon ) ) ) :
            self[particle.name] = particle
            particle.setAncestor( self, 'name' )
        elif( isinstance( particle, xParticleModule.nuclearLevel ) ) :
            baseName = particle.name.split('_')[0]
            if( 'natural' in particle.name ) :
                baseName = '_'.join(particle.name.split('_')[:2])
            if baseName in self :
                self[ baseName ].addLevel( particle )
            else:
                raise Exception( "Can't add excited level %s before ground state!" % particle.name )
        else :
            raise Exception( 'Object is not a particle: type = "%s"' % type( particle ) )

    def hasID( self, ID ) :

        for particle in self :
            if( particle.name == ID ) : return( True )
            try :
                if( particle.hasID( ID ) ) : return( True )
            except :
                pass
        return( False )

    def hasParticle( self, name ) :

        if( '_e' in name ) :
            baseName, levelId = name.split('_e')
            return( ( baseName in self ) and ( int( levelId ) in self[baseName].levels ) )
        elif( name.endswith( '_c' ) or name.endswith( '_s' ) ):
            baseName, levelId = name[:-2], name[-1]
            return( ( baseName in self ) and ( levelId in self[baseName].levels ) )
        return( name in self )

    def getParticle( self, name ) :

        try:
            if( '_e' in name ) :
                baseName, levelId = name.split('_e')
                return( self[baseName][ int(levelId) ] )
            elif( name.endswith('_c') or name.endswith('_s') ) :
                baseName, levelId = name[:-2], name[-1]
                return( self[baseName][ levelId ] )
            return( self[name] )
        except KeyError:
            raise KeyError, ("Particle '%s' is missing from the particle list!" % name)

    def names( self ) :

        return( self.values( ) )

    def check( self, info ) :
        from fudge.gnd import warning
        warnings = []

        for particle in self:
            particleWarnings = particle.check( info )
            if particleWarnings:
                warnings.append( warning.context("%s" % particle, particleWarnings) )
        return warnings

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xml = [ '%s<%s>' % ( indent, self.moniker ) ]
        particles = sorted( self.values( ) )
        for particle in particles: xml += particle.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return( xml )

def parseXMLNode( particlesElement ):
    """ starting with the xml '<particles>' element, translate XML into xParticleList """

    from pqu import PQU

    def toIntIfPossible( label ):
        try: return int(label)
        except ValueError: return label

    def fq( val ):

        if( val.startswith( 'u:' ) ) : return( xParticleModule.undefinedLevel( PQU.PQU( val[2:] ) ) )
        return PQU.PQU( val )

    def getAttrs( element, required=() ):
        conversionTable = {'mass':fq, 'energy':fq, 'probability':float, 'label':toIntIfPossible,
                'spin':xParticleModule.spin, 'parity':xParticleModule.parity}
        attrs = dict( element.items() )
        retD = {}
        for key in attrs.keys():
            if key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
            if key in required: retD[key] = attrs.pop(key)
        if attrs: retD['attributes'] = attrs
        return retD

    particles = xParticleList()
    gammas = []
    for p in particlesElement:
        if p.tag=='photon':
            particle = xParticleModule.photon( **getAttrs(p, required=('name')) )
        elif p.tag == 'lepton':
            particle = xParticleModule.lepton( **getAttrs(p, required=('name','generation','mass','charge')) )
        elif p.tag in ('isotope','FissionProduct'):
            class_ = {'isotope': xParticleModule.isotope, 'FissionProduct': xParticleModule.FissionProduct} [p.tag]
            particle = class_( **getAttrs(p, required=('name','mass')) )
            for l in p:
                level = xParticleModule.nuclearLevel( groundState=particle,
                        **getAttrs(l, required=('name','label','energy')) )
                for g in l:
                    gamma = xParticleModule.nuclearLevelGamma( angularDistribution=None,
                            **getAttrs(g, required=('finalLevel','probability')) )
                    level.addGamma( gamma )
                    gammas.append( gamma )
                particle.addLevel( level )
        else:
            raise Exception("Encountered unknown particle type '%s'" % p.tag)
        particles.addParticle( particle )
    for gamma in gammas:    # link to final levels
        gamma.finalLevel = particles.getParticle( gamma.finalLevel )
    return particles
