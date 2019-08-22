from PoPs.groups import misc as chemicalElementMiscPoPsModule

# Crack the isotope name to get the ZA
def getZAFromGNDSName( name ) :

    sym, A, m = getSymAFromGNDSName( name )
    Z = chemicalElementMiscPoPsModule.symbolFromZ[sym]
    return Z*1000+int(A)

# Crack the isotope name to get the A & symbol
def getSymAFromGNDSName( name ):
    sym = ''
    A = ''
    if '_' in name: m = name.split('_')[1]
    else: m = None
    for c in name.split('_')[0]:
        if c.isalpha(): sym += c
        else: A += c
    if sym == 'n': return sym, 1, None
    if sym == 'g': return sym, 0, None
    if m =='natural': A = '0'
    return sym, A, m
