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

import _getCoulombWavefunctions
import numpy

def getCoulombWavefunctions(rho, eta, L):
    # .... Check arguments:
    NE = len(rho)
    for array in (rho,eta):
        if type(array) != numpy.ndarray: raise TypeError("rho and eta must be numpy arrays!")
        if array.dtype != numpy.float64: raise TypeError("rho and eta must be float arrays!")
        if array.shape not in ((NE,), (NE,1)):
            raise TypeError("rho and eta must be 1-d arrays, of the same length!")
    if type(L) is not int: raise TypeError("L must be an integer!")

    F,G = _getCoulombWavefunctions.getCoulombWavefunctions( rho, eta, L )
    F.shape = G.shape = rho.shape
    return (F,G)

#==== Tests ====

def test_getCoulombWavefunctions():
    rho = numpy.arange(1.0,6.0)
    for l in range(0,21):
        print (" L =          %2i"%l)
        for i in range(1,21):
            eta = numpy.array( [0.5*i] * len(rho) )
            F,G = getCoulombWavefunctions(rho,eta,l)
            print ("%25.17E"*len(rho)) % tuple(G)   # print F or G
        print ' '

if __name__ == '__main__':
    test_getCoulombWavefunctions()
