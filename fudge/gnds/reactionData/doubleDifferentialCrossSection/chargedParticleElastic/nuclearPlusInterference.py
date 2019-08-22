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
nuclearPlusInterference separates charged-particle elastic scattering into two terms:
 - (cross section) - (pure Rutherford cross section), integrated from mu = -1 up to a cut-off angle mu_cutoff
 - the distribution (P(mu|E) - P(mu|E)_Rutherford) / cross section), defined for mu = -1 to mu_cutoff
"""

import math

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule

from fudge.gnds.reactionData import crossSection as crossSectionModule
from fudge.gnds.productData.distributions import angular as angularModule

from .misc import chargedParticleElasticTerm


__metaclass__ = type

class crossSection( chargedParticleElasticTerm):

    moniker = 'crossSection'

    allowedDataForms = (crossSectionModule.XYs1d, crossSectionModule.regions1d)

class distribution( chargedParticleElasticTerm):

    moniker = 'distribution'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)


class nuclearPlusInterference( ancestryModule.ancestry ):

    moniker = "nuclearPlusInterference"

    def __init__( self, muCutoff, crossSection = None, distribution = None ):

        ancestryModule.ancestry.__init__( self )
        self.__muCutoff = muCutoff
        self.crossSection = crossSection
        self.distribution = distribution

    @property
    def muCutoff(self): return self.__muCutoff

    @property
    def crossSection( self ):
        return self.__crossSection

    @crossSection.setter
    def crossSection(self, value):
        if not isinstance(value, crossSection):
            raise TypeError( '%s cross section must be instance of %s' % (self.moniker, crossSection.moniker))
        value.setAncestor( self )
        self.__crossSection = value

    @property
    def distribution( self ):
        return self.__distribution

    @distribution.setter
    def distribution(self, value):
        if not isinstance(value, distribution):
            raise TypeError( '%s distribution must be instance of %s' % (self.moniker, distribution.moniker))
        value.setAncestor( self )
        self.__distribution = value

    @property
    def domainMin(self): return self.crossSection.data.domainMin

    @property
    def domainMax(self): return self.crossSection.data.domainMax

    @property
    def domainUnit(self): return self.crossSection.data.domainUnit

    def check( self, info ):
        """
        Things to check: mu_cutoff should be between 0.99 and 1.0,
        domain of effective xsc / distribution should match,
        if identical particles distribution should only be for 0->muCutoff, otherwise -1->muCutoff
        distribution should be normalized
        If the cross section goes negative, also check that the pure Coulomb portion is large enough to 'win'
        :return:
        """
        return []
        raise NotImplementedError

    def convertUnits( self, unitMap ):

        self.crossSection.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
        :param E: incident energy
        :param mu: scattering angle cosine
        :return: differential cross section at E,mu (integrated over phi) in b/sr
        """

        if mu > self.muCutoff: return 0
        if E < self.crossSection.data.domainMin: return 0
        if E > self.crossSection.data.domainMax:
            raise ValueError( "Attempted to evaluate at %s, outside of domain limit %s %s" % ( E, self.crossSection.data.domainMax, self.crossSection.data.domainUnit ) )
        return self.crossSection.data.evaluate( E ) * self.distribution.data.evaluate( E ).evaluate( mu ) / (2*math.pi)

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s muCutoff="%s">' % (indent, self.moniker, self.muCutoff) ]
        xml += self.crossSection.toXMLList( indent2, **kwargs )
        xml += self.distribution.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag == crossSection.moniker:
                Xsc = crossSection.parseXMLNode( child, xPath, linkData)
            elif child.tag == distribution.moniker:
                Dist = distribution.parseXMLNode( child, xPath, linkData)
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        NPCI = nuclearPlusInterference(
                float(element.get('muCutoff')), crossSection=Xsc, distribution=Dist )

        xPath.pop()
        return NPCI
