# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry

monikerParticles = 'particles'

class xParticleList( dict, ancestry ) : 

    def __init__( self ) :

        ancestry.__init__(self, monikerParticles, None)

    def addParticle( self, particle ) :
        """Add either a particle or an excited level to xParticleList """

        if( 'mass' in particle.attributes ) :               # it's a basic particle
            self[ particle.name ] = particle
            particle.setParent( self )
        else:
            baseName = particle.name.split('_')[0]
            if( 'natural' in particle.name ) :
                baseName = '_'.join(particle.name.split('_')[:2])
            elif( 'FissionProduct' in particle.name ) :
                baseName = '_'.join(particle.name.split('_')[:2])
            if baseName in self :
                self[ baseName ].addLevel( particle )
            else:
                raise Exception( "Can't add excited level %s before ground state!" % particle.name )

    def hasParticle( self, name ) :

        if( '_e' in name ) :
            baseName, levelId = name.split('_e')
            return( ( baseName in self ) and ( int( levelId ) in self[baseName].levels ) )
        elif( name.endswith( '_c' ) or name.endswith( '_s' ) ):
            baseName, levelId = name[:-2], name[-1]
            return( ( baseName in self ) and ( levelId in self[baseName].levels ) )
        elif( '_m' in name ):
            baseName, m_index = name.split('_m')
            return( ( baseName in self ) and ( bool( [ lev for lev in self[baseName].levels.values()
                if( lev.attributes.get('metaStableIndex' ) == int( m_index ) ) ] ) ) )
        return( name in self )

    def getParticle( self, name ) :

        try:
            if( '_e' in name ) :
                baseName, levelId = name.split('_e')
                return( self[baseName][ int(levelId) ] )
            elif( name.endswith('_c') or name.endswith('_s') ) :
                baseName, levelId = name[:-2], name[-1]
                return( self[baseName][ levelId ] )
            elif( '_m' in name ) :
                baseName, m_index = name.split('_m')
                for lev in self[baseName].levels.values() :
                    if( lev.getMetaStableIndex() == int(m_index) ) : return( lev )
            return( self[name] )
        except KeyError:
            raise KeyError, ("Particle '%s' is missing from the particle list!" % name)

    def names( self ) :

        return( self.values( ) )

    def check( self, info ) :
        from fudge.gnd import warning
        warnings = []

        for particle in self:
            particleWarnings = self.getParticle(particle).check( info )
            if particleWarnings:
                warnings.append( warning.context("%s" % particle, particleWarnings) )
        return warnings

def parseXMLNode( particlesElement ):
    """ starting with the xml '<particles>' element, translate XML into xParticleList """
    from pqu import physicalQuantityWithUncertainty
    from . import xParticle
    def toIntIfPossible( label ):
        try: return int(label)
        except ValueError: return label
    def fq( val ):
        if val.startswith('u:'):
            from fudge.gnd import miscellaneous
            return miscellaneous.undefinedLevel( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( val[2:] ) )
        return physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( val )
    def getAttrs(element, required=None):
        conversionTable = {'mass':fq, 'energy':fq, 'probability':float, 'label':toIntIfPossible,
                'spin':xParticle.spin, 'parity':xParticle.parity}
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
        particle = xParticle.xParticle( **getAttrs(p, required=('name','genre','mass')) )
        for l in p:
            level = xParticle.nuclearLevel( groundState=particle,
                    **getAttrs(l, required=('name','label','energy')) )
            for g in l:
                gamma = xParticle.nuclearLevelGamma( angularDistribution=None,
                        **getAttrs(g, required=('finalLevel','probability')) )
                level.addGamma( gamma )
                gammas.append( gamma )
            particle.addLevel( level )
        particles.addParticle( particle )
    for gamma in gammas:    # link to final levels
        gamma.finalLevel = particles.getParticle( gamma.finalLevel )
    return particles
