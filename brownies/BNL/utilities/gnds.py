from PoPs.groups import misc as chemicalElementMiscPoPsModule
from fudge.reactionData.crossSection import component


# Crack the isotope name to get the ZA
def getZAFromGNDSName(name):
    sym, A, m = getSymAFromGNDSName(name)
    Z = chemicalElementMiscPoPsModule.symbolFromZ[sym]
    return Z * 1000 + int(A)


# Crack the isotope name to get the A & symbol
def getSymAFromGNDSName(name):
    sym = ''
    A = ''
    if '_' in name:
        m = name.split('_')[1]
    else:
        m = None
    for c in name.split('_')[0]:
        if c.isalpha():
            sym += c
        else:
            A += c
    if sym == 'n':
        return sym, 1, None
    if sym == 'g':
        return sym, 0, None
    if m == 'natural':
        A = '0'
    return sym, A, m


def check_is_cross_section(x):
    if not isinstance(x, component):
        raise TypeError("Not instance of fudge.reactionData.crossSection.component")
