# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Functions for dealing with polynomials:
    -degree( poly ) -> returns the highest-order term
    -poly_div(N, D) -> returns (quotient, remainder) of polynomial division
    -poly_derivative( poly ) -> return derivative of a polynomial
    -sturm_sequence( poly ) -> returns (F1,F2,F3,F4,...): original function, derivative, -rem(F1/F2), -rem(F2/F3)
        The Sturm sequence is an efficient way to locate roots of the polynomial

    polynomials are handled as lists of coefficients: [1,0.5,3] == 1 + 0.5*x + 3*x**2
"""

import math
from itertools import izip

""" degree() and poly_div() taken from http://rosettacode.org/wiki/Polynomial_long_division
Small change: don't take absolute value in poly_div, need to preserve original sign! """
def degree(poly):
    while poly and poly[-1] == 0:
        poly.pop()   # normalize
    return len(poly)-1
 
def poly_div(N, D):
    dD = degree(D)
    dN = degree(N)
    if dD < 0: raise ZeroDivisionError
    if dN >= dD:
        q = [0] * dN
        while dN >= dD:
            d = [0]*(dN - dD) + D
            mult = q[dN - dD] = N[-1] / float(d[-1])
            d = [coeff*mult for coeff in d]
            #N = [math.fabs ( coeffN - coeffd ) for coeffN, coeffd in izip(N, d)]
            N = [( coeffN - coeffd ) for coeffN, coeffd in izip(N, d)]
            dN = degree(N)
        r = N
    else:
        q = [0]
        r = N
    return q, r

def poly_derivative( poly ):
    result = []
    for i in range(1,len(poly)):
        result.append( i*poly[i] )
    return result

def sturm_sequence( poly ):
    """Get the Sturm sequence for a polynomial.
    This sequence is an easy way to detect polynomial roots in a given interval by evaluating each function
    at the two endpoints. Better than Newton's method, which only detects odd number of roots in the interval.

    First and second terms are the original polynomial F1 and its derivative F2.
    F3 = -(F1 modulo F2), F4 = -(F2 modulo F3) and so on (seq. for polynomial of order N has N+1 terms)

    For example, the polynomial x^3 + 3x^2 - 1 has two roots in the interval -1->1, which we can detect by
    evaluating the Sturm sequence:
        endpoint    F1  F2  F3  F4  cumulative sign changes:
           -1.0     1  -3  -1   9/4         2
            1.0     3   9   3   9/4         0
    The difference in cumulative sign changes between the two endpoints gives the number of roots in the interval"""

    sequence = [ poly, poly_derivative( poly ) ]
    for i in range( degree(poly)-1 ):
        # add next term: negative remainder from ratio of previous two terms
        nextTerm = [-1*val for val in poly_div( sequence[-2], sequence[-1] )[1]]
        if nextTerm == []: break
        sequence.append( nextTerm )

    return [generate_function( P ) for P in sequence]

def generate_function( poly ):
    """Convert list of polynomial coefficients into a function. """
    def func( x ):
        result = 0
        for i in range(len(poly)):
            result += poly[i] * math.pow(x,i)
        return result
    return func

if __name__ == '__main__':
    print "POLYNOMIAL LONG DIVISION"
    N = [-42, 0, -12, 1]
    D = [-3, 1, 0, 0]
    print "  %s / %s =" % (N,D),
    print " %s remainder %s" % poly_div(N, D)

