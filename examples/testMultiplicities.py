#!/usr/bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

""" demonstrates checking data in the new GND format.

python testMultiplicities.py *.xml  # to check all translated gnd files in a directory """

import os
import sys
from xml.etree import cElementTree

def getData( data ):
    try:
        if 'piecewise' in data.tag.lower():
            dat = [region[1].text for region in data.findall(".//region")]
        elif data.tag=="pointwise":
            dat = [data.find("data").text]
        else:
            #print ("unknown data type: %s" % data.tag)
            return []
        dat = [map(float, d.split()) for d in dat]
        return [(d[::2],d[1::2]) for d in dat]
    except Exception, e:
        print("error trying to extract data")
        return []

def checkReaction( reac ):
    """ look for negative cross sections or multiplicities """
    warnings = []
    
    xsc = reac.find("crossSection")
    data = xsc.find( xsc.get("nativeData") )
    for x,y in getData( data ):
        if min(y)<0: warnings.append("negative cross section")
        if len(x)!=len(set(x)): warnings.append("duplicate cross section energies")

    for prod in reac.findall(".//product"):
        mult = prod.get("multiplicity") # attribute
        if mult=="energyDependent":
            mult = prod.find("multiplicity") # element
            data = mult.find( mult.get("nativeData") )
            for x,y in getData( data ):
                if min(y)<0: warnings.append("negative energy-dependent multiplicity")
                if len(x)!=len(set(x)): warnings.append("duplicate multiplicity energies")
        elif mult=="unknown":
            warnings.append( "unknown multiplicity" )
        else:
            mult = float(mult)
            if mult<=0: warnings.append( "negative multiplicity" )
    return warnings

if __name__ == '__main__':
    files = sys.argv[1:]
    for filename in files:
        if '-covar' in filename: continue   # ignore covariances for now

        xdoc = cElementTree.parse(filename)
        for reac in xdoc.findall( "reaction" ):
            for warning in checkReaction( reac ):
                print( "%s in %s, reaction %s (label %s)" % (warning, 
                    os.path.split(filename)[-1], reac.get("outputChannel"), reac.get("label") ) )

