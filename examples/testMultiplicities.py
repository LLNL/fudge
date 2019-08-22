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

