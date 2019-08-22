#!/usr/bin/env python

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

