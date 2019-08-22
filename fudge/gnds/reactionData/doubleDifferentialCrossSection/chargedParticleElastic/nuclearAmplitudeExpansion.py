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
"""
The nuclearAmplitudeExpansion represents charged-particle elastic scattering using four terms:
 - pure Coulomb scattering  (i.e. RutherfordScattering),
 - pure nuclear scattering
 - complex interference between Coulomb and nuclear terms, broken into real and imaginary components
"""

from __future__ import division

import math
import numpy

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData.series1d import Legendre, LegendreSeries

from fudge.gnds.productData.distributions import angular as angularModule
from .misc import chargedParticleElasticTerm

__metaclass__ = type

class nuclearTerm( chargedParticleElasticTerm):

    moniker = 'nuclearTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class realInterferenceTerm( chargedParticleElasticTerm):

    moniker = 'realInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class imaginaryInterferenceTerm( chargedParticleElasticTerm):

    moniker = 'imaginaryInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)


class nuclearAmplitudeExpansion( ancestryModule.ancestry):

    moniker = 'nuclearAmplitudeExpansion'

    def __init__( self, nuclearTerm = None, realInterference = None, imaginaryInterference = None ):

        ancestryModule.ancestry.__init__( self )
        self.nuclearTerm = nuclearTerm
        self.realInterferenceTerm = realInterference
        self.imaginaryInterferenceTerm = imaginaryInterference

    @property
    def etaCoefficient( self ) :

        return( self.ancestor.etaCoefficient )

    @property
    def nuclearTerm( self ):
        return self.__nuclear

    @nuclearTerm.setter
    def nuclearTerm(self, value):
        if not isinstance(value, nuclearTerm):
            raise TypeError( '%s nuclear term must be instance of %s' % (self.moniker, nuclearTerm.moniker))
        value.setAncestor( self )
        self.__nuclear = value

    @property
    def realInterferenceTerm( self ):
        return self.__realInterference

    @realInterferenceTerm.setter
    def realInterferenceTerm(self, value):
        if not isinstance(value, realInterferenceTerm):
            raise TypeError( '%s real interference term must be instance of %s' %
                             (self.moniker, realInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__realInterference = value

    @property
    def imaginaryInterferenceTerm( self ):
        return self.__imaginaryInterference

    @imaginaryInterferenceTerm.setter
    def imaginaryInterferenceTerm(self, value):
        if not isinstance(value, imaginaryInterferenceTerm):
            raise TypeError( '%s imaginary interference term must be instance of %s' %
                             (self.moniker, imaginaryInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__imaginaryInterference = value

    @property
    def domainMin(self) :

        return self.nuclearTerm.data.domainMin

    @property
    def domainMax(self) :

        return self.nuclearTerm.data.domainMax

    @property
    def domainUnit(self) :

        return self.nuclearTerm.data.domainUnit

    def check( self, info ):
        """
        Check that incident energy domain of all three terms are the same
        FIXME: I'm not sure whether each term is meant to be normalized  (I see L=0 terms != 0)
        There may be more checks we can perform depending on whether particles are identical
        :return:
        """
        return []
        raise NotImplementedError()

    def convertUnits( self, unitMap ):

        self.nuclearTerm.convertUnits( unitMap )
        self.realInterferenceTerm.convertUnits( unitMap )
        self.imaginaryInterferenceTerm.convertUnits( unitMap )

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
        :param E: incident energy
        :param mu: scattering angle cosine
        :return: differential cross section at E,mu (integrated over phi) in b/sr
        """

        if E < self.nuclearTerm.data.domainMin: return 0
        if E > self.nuclearTerm.data.domainMax:
            raise ValueError( "Attempted to evaluate at %s, outside of domain limit %s %s" % ( E, self.nuclearTerm.data.domainMax, self.nuclearTerm.data.domainUnit ) )


        eta = self.etaCoefficient / numpy.sqrt( E )

        # get Legendre series at desired energy, interpolating coefficients if necessary
        coefficients = {}
        for term in ( 'nuclearTerm', 'realInterferenceTerm', 'imaginaryInterferenceTerm' ) :
            data = getattr(self, term).data
            position, function1, function2, frac, interpolation, interpolationQualifier = data.getBoundingSubFunctions( E )
            if position == '<':
                coefficients[term] = [0]
            elif position == '=':
                coefficients[term] = function1.coefficients
            elif position == '>':
                raise ValueError("Attempted to evaluate outside of domain")
            else :                                  # Ignore interpolation as it may be "charged-particle" which does not make sense.
                coef1 = function1.coefficients[:]
                coef2 = function2.coefficients[:]
                coef1 += max( 0, len( coef2 ) - len( coef1 ) ) * [ 0 ]
                coef2 += max( 0, len( coef1 ) - len( coef2 ) ) * [ 0 ]
                coefficients[term] = [ ( 1.0 - frac ) * coef1[idx] + frac * coef2[idx] for idx in range( len( coef1 ) ) ]

        if self.findAttributeInAncestry('identicalParticles'):
            NL = len(coefficients['nuclearTerm'])
            NL_half = len(coefficients['realInterferenceTerm'])

            P_lAtMu = [Legendre(l, mu) for l in range(NL)]

            term1 = (1+mu) * numpy.exp(1j * eta * numpy.log((1-mu)/2))
            term2 = (1-mu) * numpy.exp(1j * eta * numpy.log((1+mu)/2))
            interferenceCoefficients = [ term1 + (-1)**ell * term2 for ell in range(NL) ]
            complexAs = [ real + 1j * imag for (real,imag) in zip(
                coefficients['realInterferenceTerm'], coefficients['imaginaryInterferenceTerm']) ]
            interferenceSum = numpy.sum( [interferenceCoefficients[ell] * (2*ell+1)/2 * complexAs[ell] * P_lAtMu[ell]
                                          for ell in range(NL_half)])
            interference = -2*eta / (1 - mu**2) * numpy.real( interferenceSum )

            nuclear = numpy.sum( [(4*ell+1)/2 * coefficients['nuclearTerm'][ell] * P_lAtMu[ell]
                                  for ell in range(NL)] )

        else:   # target != projectile
            RI = LegendreSeries( coefficients['realInterferenceTerm'] ).evaluate( mu )
            II = LegendreSeries( coefficients['imaginaryInterferenceTerm'] ).evaluate( mu )

            interference = -2*eta / (1-mu) * numpy.real(numpy.exp(1j * eta * numpy.log((1-mu)/2)) * (RI + 1j*II))

            nuclear = LegendreSeries( coefficients['nuclearTerm'] ).evaluate( mu )

        return( nuclear + interference )

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s>' % (indent, self.moniker) ]
        xml += self.nuclearTerm.toXMLList( indent2, **kwargs )
        xml += self.realInterferenceTerm.toXMLList( indent2, **kwargs )
        xml += self.imaginaryInterferenceTerm.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag == nuclearTerm.moniker:
                nuclear = nuclearTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == realInterferenceTerm.moniker:
                real = realInterferenceTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == imaginaryInterferenceTerm.moniker:
                imag = imaginaryInterferenceTerm.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        CoulExp = nuclearAmplitudeExpansion( nuclearTerm=nuclear, realInterference=real, imaginaryInterference=imag)
        xPath.pop()
        return CoulExp
