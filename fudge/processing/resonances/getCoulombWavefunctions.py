#!/usr/bin/env python

# <<BEGIN-copyright>>
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
