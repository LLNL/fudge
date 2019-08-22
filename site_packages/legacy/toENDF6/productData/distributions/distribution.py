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

import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import fudge.gnd.productData.distributions.distribution as distributionModule
import fudge.gnd.productData.distributions.reference as referenceModule
import fudge.gnd.productData.distributions.unspecified as unspecifiedModule

def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    if( len( self ) == 0 ) :
        raise Exception( 'This needs to be tested' )     # In particular, the method endfMultiplicityList has been removed.
        if( MT in [ 452, 455, 456 ] ) : return
        targetInfo['MF6LCTs'].append( None )
        particleID = targetInfo[ 'zapID' ]    # get ZAP for the outgoing particle
        if( particleID == 'gamma' ) :
            ZAP = 0
        else :
            Z, A, suffix, ZAP = nuclear.getZ_A_suffix_andZAFromName( particleID )
        nPoints, multiplicityList = targetInfo['multiplicity'].endfMultiplicityList( targetInfo )
        endfMFList[6][MT] += [ endfFormatsModule.endfContLine( ZAP, targetInfo['particleMass'] / targetInfo['neutronMass'], 0, 0, 1, nPoints ) ]
        endfMFList[6][MT] += multiplicityList
        return
    form = self[targetInfo['style']]
    if( hasattr( form, 'toENDF6' ) ) :
        form.toENDF6( MT, endfMFList, flags, targetInfo )
    elif( isinstance( form, referenceModule.form ) ) :
        pass
    elif( isinstance( form, unspecifiedModule.form ) ) :
        pass
    else :
        print 'WARNING: Distribution, no toENDF6 for class = %s' % form.__class__

distributionModule.component.toENDF6 = toENDF6
