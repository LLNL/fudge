# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import math

from LUPY import ancestry as ancestryModule

from xData import axes as axesModule


def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

    axes = axesModule.Axes(3)
    axes[0] = axesModule.Axis( 'dSigma_dMu', 0, crossSectionUnit )
    axes[1] = axesModule.Axis( 'mu', 1, '' )
    axes[2] = axesModule.Axis( 'energy_in', 2, energyUnit )

    return( axes )

class ChargedParticleElasticTerm(ancestryModule.AncestryIO, abc.ABC):
    """
    In addition to the Rutherford scattering term, charged-particle elastic sections
    may have additional terms: 'crossSection', 'distribution',
    'nuclearTerm', etc. Each is just a context layer that contains the actual data.
    This serves as the base class for all terms.
    """

    @property
    @abc.abstractmethod
    def allowedDataForms(self): pass

    def __init__(self, data):

        ancestryModule.AncestryIO.__init__(self)
        self.data = data

    @property
    def data(self): return self.__data

    @data.setter
    def data(self, value):
        if not isinstance(value, self.allowedDataForms):
            raise TypeError("'%s' cannot contain data of type '%s'" % (self.moniker, type(value)))
        value.setAncestor( self )
        self.__data = value

    def convertUnits( self, unitMap ):

        self.data.convertUnits( unitMap )

    def toXML_strList(self, indent='', **kwargs):

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXML_strList(indent+'  ', **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        child = element[0]
        parseClass = None
        for allowedForm in cls.allowedDataForms:
            if child.tag == allowedForm.moniker:
                parseClass = allowedForm
                break
        if parseClass is None:
            raise TypeError("Unexpected element '%s' encountered inside %s" % (child.tag, cls.moniker))
        data = parseClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        term_ = cls( data )

        xPath.pop()

        return term_

def fixMuRange(data, epsilon=1e-6):
    '''This function ensures that the range of mu in *data* is [-1, 1] by adding 0.0's to the missing ranges. Note, *data* is modified.'''

    m1 = data[0][0]
    if m1 > -1.0:
        if abs(m1) == 0.0:
            m2 = -epsilon
        else:
            eps = math.copysign(epsilon, -m1)
            m2 = (1 + eps) * m1
            if m2 <= -1.0:
                data[0][0] = (1 - eps) * m1
                m2 = m1
        data.insert(0, [m2, 0.0])
        data.insert(0, [-1.0, 0.0])

    m1 = data[-1][0]
    if m1 < 1.0:
        if abs(m1) == 0.0:
            m2 = epsilon
        else:
            eps = math.copysign(epsilon, m1)
            m2 = (1 + eps) * m1
            if m2 >= 1.0:
                data[-1][0] = (1 - eps) * m1
                m2 = m1
        data.append([m2, 0.0])
        data.append([1.0, 0.0])
