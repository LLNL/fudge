# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc

from xData import ancestry as ancestryModule
from xData import axes as axesModule

__metaclass__ = type

def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

    axes = axesModule.axes( rank = 3 )
    axes[0] = axesModule.axis( 'dSigma_dMu', 0, crossSectionUnit )
    axes[1] = axesModule.axis( 'mu', 1, '' )
    axes[2] = axesModule.axis( 'energy_in', 2, energyUnit )

    return( axes )

class chargedParticleElasticTerm( ancestryModule.ancestry):
    """
    In addition to the Rutherford scattering term, charged-particle elastic sections
    may have additional terms: 'crossSection', 'distribution',
    'nuclearTerm', etc. Each is just a context layer that contains the actual data.
    This serves as the base class for all terms.
    """
    __metaclass__ = abc.ABCMeta

    @property
    @abc.abstractmethod
    def allowedDataForms(self): pass

    def __init__(self, data):
        ancestryModule.ancestry.__init__(self)
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

    def toXMLList(self, indent='', **kwargs):

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXMLList(indent+'  ', **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        child = element[0]
        parseClass = None
        for allowedForm in cls.allowedDataForms:
            if child.tag == allowedForm.moniker:
                parseClass = allowedForm
                break
        if parseClass is None:
            raise TypeError("Unexpected element '%s' encountered inside %s" % (child.tag, cls.moniker))
        data = parseClass.parseXMLNode( child, xPath, linkData )
        term_ = cls( data )
        xPath.pop()
        return term_
