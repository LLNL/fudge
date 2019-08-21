# <<BEGIN-copyright>>
# <<END-copyright>>

"""
containers for a particle and a reaction, plus utility to parse in string
representations to containers:
>>>parseReaction("n + Fe56 -> n[multiplicity:'2'] + (Fe55_s -> gamma)")
"""

from pyparsing import *

__metaclass__ = type

class particleString:
    optionList = ('multiplicity','emissionMode','decayRate')
    specialNames = ('gamma','electron') # treat differently
    def __init__(self, symbol, A, excitation=0, opts=None, decaysTo=None):
        self.symbol = symbol
        self.A = A
        self.excitation = excitation
        self.opts = opts or {}
        self.decaysTo = decaysTo or []
    
    def __str__(self):
        ret = ''
        if self.decaysTo: ret += '('
        if self.symbol in particleString.specialNames:
            ret += self.symbol
        else:
            ret += '%s%i' % (self.symbol, self.A)
        if self.excitation == 's': ret += '_s'
        elif self.excitation: ret += '_e%i' % self.excitation
        opts = ["%s:'%s'" % (key,self.opts[key]) for key in particleString.optionList
                if self.opts.get(key) is not None]
        if opts:
            ret += "[" + ', '.join(opts) + "]"
        if self.decaysTo:
            ret += ' -> ' + ' + '.join([str(a) for a in self.decaysTo]) + ')'
        return ret

class reactionString:
    def __init__(self, projectile, target, products, info=None):
        self.projectile = projectile
        self.target = target
        self.products = products
        self.info = info
    
    def __str__(self):
        ret = '%s + %s -> ' % (self.projectile,self.target)
        ret += ' + '.join( [str(p) for p in self.products] )
        if self.info: ret += ' [%s]' % self.info
        return ret

def particleParser():
    """ 
    parse string of form "Pu239_e2[option1:'value', option2:'value' ... ]" into particle class
    """
    integer = Word(nums).setParseAction( lambda t: int(t[0]) )
    optionDict = (Suppress('[') + delimitedList( Word(alphas)+Suppress(':') +
            Suppress("'")+Word(alphanums+'+-.')+Suppress("'") , delim=",") + 
            Suppress(']') )
    optionDict.setParseAction( lambda t: dict(zip( t[::2],t[1::2] )) )
    
    # now define particle: name, excitation, [list of options]
    particle = (
            Word(alphas)+integer +
            Optional( Suppress('_') + ((Suppress('e')+integer) | Word('s')), default=0 ) +
            Optional( optionDict, default={} )
            )
    # special cases: gamma and electron have no A or excitation:
    specialParticle = (
        (Literal('gamma')|Literal('electron'))+Optional(optionDict, default={})
        ).setParseAction( lambda t: [t[0],0,0,t[1]] )
    particleParser = particle | specialParticle
    #particleParser.setParseAction( lambda t: particleString(*t) )
    return particleParser

def reactionParser():
    """
    parse reaction of form "n + Fe56 -> n[options...] + (Fe55_u -> gamma)" into reaction class
    """
    # reaction string may include one or more decays,
    # of form (Th232 -> He4 + (Ra228 -> He4 + Rd224)).
    # decays may contain subsequent decays, so they are defined recursively:
    decay = Forward()
    particle = particleParser()
    atom = Group(particle) | Group( Suppress('(') + decay + Suppress(')') )
    decay << Group(particle) + Suppress('->') + Group(delimitedList(atom,delim='+'))
    
    inputChannel = Group(particle)+Suppress('+')+Group(particle)
    outputChannel = delimitedList(atom, delim='+')
    # extra channel information, "[total fission]" for example:
    info = Optional(Suppress('[') + Word(alphas+' ') + Suppress(']'), default='')
    return inputChannel+Literal('->').suppress() + Group(outputChannel) + info

def parseReaction( str ):
    reaction = reactionParser()
    proj, targ, outgoing, info = reaction.parseString(str)
    proj, targ = [particleString(*a) for a in (proj,targ)]
    products = []
    def readProduct(element):
        if len(element)==4:
            return particleString( *element )
        elif len(element)==2:   # decay
            parent = particleString(*element[0])
            for daughter in element[1]:
                parent.decaysTo.append( readProduct( daughter ) )
            return parent
        else:
            raise Exception, "parsing problem!"
    for prod in outgoing:
        products.append( readProduct(prod) )
    
    return reactionString(proj, targ, products, info)

