# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import fudge.processing.resonances.getScatteringMatrices, unittest, numpy

def python_getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities):
    NE = len(E)
    dim = len(widths)
    dE = Eres - E[:,numpy.newaxis]
    DEN = 1/(dE**2 + captureWidth**2)
    captOverDEN = captureWidth*DEN
    dEoverDEN = dE*DEN
    del dE, DEN

    R = numpy.zeros((NE,dim,dim))
    S = numpy.zeros((NE,dim,dim))
    for i in range(dim):
        for j in range(dim):
            width_ij = widths[i]*widths[j]
            p_ij = penetrabilities[i]*penetrabilities[j]
            R[:,i,j] = p_ij.T * numpy.sum( width_ij * captOverDEN, axis=1)
            S[:,i,j] = p_ij.T * numpy.sum( width_ij * dEoverDEN, axis=1)
            # symmetrize:
            R[:,j,i] = R[:,i,j]
            S[:,j,i] = S[:,i,j]
    return R,S

if __name__ == '__main__':

    class Test_Driver( unittest.TestCase ):
    
        def test_getScatteringMatrices( self ):
            E = numpy.array([1e-5, 0.025, 0.1, 1.0, 10.0, 20.0])
            Eres = numpy.array([0.35, 2.4, 7.8])
            captureWidth = numpy.array([0.08, 1.2, 1.6])
            widths = numpy.array([[0.01,0.3,0.7],
                                [0,0,1.2e-4]])
            penetrabilities = numpy.array([[0.1,0.2,0.3,0.4,0.5,0.55],
                                        [0.05,0.15,0.25,0.35,0.45,0.55]])

            R,S = fudge.processing.resonances.getScatteringMatrices.getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            R1,S1 = python_getScatteringMatrices(E,Eres,captureWidth,widths,penetrabilities)
            self.assertTrue( numpy.allclose( R,R1 ) and numpy.allclose( S,S1 ) )

    unittest.main()

