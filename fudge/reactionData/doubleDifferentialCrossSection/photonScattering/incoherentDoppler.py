# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Doppler broadened photon incoherent double difference cross section form and its supporting classes.
Data only covers Bound->Free reactions that ionize the electron.
"""

import math

from pqu import PQU as PQUModule
from numericalFunctions import integration as nf_integrationModule

from LUPY import ancestry as ancestryModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import suites as suitesModule

from xData import enums as xDataEnumsModule
from xData import xs_pdf_cdf as xs_pdf_cdfModule

from fudge.core.math import fudgemath as fudgemathModule
from ... import crossSection as crossSectionModule
from .. import base as baseModule

class XYs1d(baseModule.XYs1d):

    pass

class Xs_pdf_cdf1d(xs_pdf_cdfModule.Xs_pdf_cdf1d):

    def toLinearXYsClass( self ) :
        """
        This method returns the class for representing linear point-wise 1-d data of *self*.
        """

        return XYs1d

class ComptonProfile(ancestryModule.AncestryIO):

    moniker = 'ComptonProfile'
    ENDF_MT = 1504

    def __init__(self, data):

        ancestryModule.AncestryIO.__init__(self)

        if not isinstance(data, (XYs1d, Xs_pdf_cdf1d)):
            raise TypeError( "Invalid data type: got %s." % type(data))

        self.__data = data
        self.__data.setAncestor(self)

    @property
    def data(self):

        return self.__data

    def convertUnits(self, unitMap):

        self.__data.convertUnits(unitMap)

    def check(self, info):

        return []

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        xmlStringList += self.__data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append( node.tag )

        child = node[0]
        if( child.tag == XYs1d.moniker ) :
            data = XYs1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        elif( child.tag == Xs_pdf_cdf1d.moniker ) :
            data = Xs_pdf_cdf1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else :
            raise TypeError('Invalid child "%s" for "%s"' % (child.tag, cls.tag))

        comptonProfile = cls(data)

        xPath.pop()

        return comptonProfile

class Form( baseModule.Form ) :

    moniker = 'incoherentBoundToFreePhotonScattering'
    keyName = 'label'

    subformAttributes = ('comptonProfile', )

    def __init__(self, pid, label, productFrame, comptonData):

        if not isinstance(comptonData, (ComptonProfile)):
            raise Exception('Instance is class "%s" and not ComptonProfile' % comptonData.__class__)

        baseModule.Form.__init__(self, pid, label, productFrame, (comptonData,))

    def check(self, info):

        return []

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        aveEnergy = []
        aveMomenta = []
        aveEnergy.append([ 1e-6, 1e-6 ])
        aveEnergy.append([ 1e5, 1e5 ])
        aveMomenta.append([ 1e-6, 1e-6 ])
        aveMomenta.append([ 1e5, 1e5 ])
        return( [ aveEnergy ], [ aveMomenta ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        from fudge.productData.distributions import photonScattering as photonScatteringModule

        data = self.comptonProfile.data
        if not isinstance(data, XYs1d):
            raise TypeError('Processing of type %s not supported.' % type(data))

        xs_opdf_cdf1d = Xs_pdf_cdf1d.fromXYs(data, thinEpsilon=1e-14)       # fromXYs will normalize first.
        comptonProfile = ComptonProfile(xs_opdf_cdf1d)

        doubleDifferentialForm = Form(self.pid, style.label, self.productFrame, comptonProfile)
        self.ancestor.add(doubleDifferentialForm)

        return photonScatteringModule.IncoherentBoundToFreePhotonScattering.Form(label=style.label, link=doubleDifferentialForm, relative=True)

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( None )

    '''
    def crossSection( self, domainMin, domainMax, domainUnit, tolerance=1e-3, interpolation=xDataEnumsModule.Interpolation.linlin):
        """
        Integrates the double differential cross section to get the cross section :math:`\\sigma(E)` where :math:`E` is the projectile's energy.
        """

        class Parameters :

            pass

        class Tester :

            def __init__( self, parameters, relativeTolerance, absoluteTolerance ) :

                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                return( integratedCrossSection( energy, self.parameters ) )

        def integrand( mu, parameters ) :

            oneMinusMu = 1.0 - mu
            _x = parameters.energy_hc * math.sqrt( 0.5 * oneMinusMu )
            kPrime_k = 1.0 / ( 1.0 + parameters.k * oneMinusMu )

            integral = kPrime_k * kPrime_k * ( 1.0 + mu * mu + kPrime_k * parameters.k**2 * oneMinusMu**2 )
            integral *= parameters.scatteringFactor.evaluate( _x )

            return( integral )

        def integratedCrossSection( energy, parameters ) :

            parameters.energy_hc = energy * parameters.factor_E2x
            parameters.k = energy * factor_electronMass
            return( parameters.constant * nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0] )

        scatteringFactor = self.scatteringFactor

        factor_electronMass = 1.0 / PQUModule.PQU( 1.0, 'me * c**2' ).getValueAs( domainUnit )
        factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % domainUnit ).getValueAs( scatteringFactor.axes[1].unit )

        parameters = Parameters( )
        tolerance = min( 0.1, max( tolerance, 1e-6 ) )
        parameters.tolerance = 1e-2 * tolerance
        parameters.factor_E2x = factor_E2x
        parameters.scatteringFactor = scatteringFactor
        electronRadius = PQUModule.PQU( 1.0, 'e**2 / ( 4 * pi * eps0 * me * c**2 )' ).getValueAs( '1e-12 * cm' )
        parameters.constant = math.pi * electronRadius**2

        domainMax = min( domainMax, scatteringFactor.domainMax / factor_E2x )
        energies = [ domainMin, math.sqrt( domainMin * domainMax ), domainMax ]
        energy_crossSection = []
        for energy in energies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters ) ] )
        _tester = Tester( parameters, tolerance, tolerance * energy_crossSection[-1][1] )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 20, interpolation = interpolation )

        return( crossSectionModule.XYs1d( data = energy_crossSection, axes = crossSectionModule.defaultAxes( domainUnit ),
                interpolation = interpolation ) )
    '''

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)

        # This implementation is not robust
        child = node[0]
        if child.tag == ComptonProfile.moniker:
            comptonProfile = ComptonProfile.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( child.tag, form.moniker ) )

        form = cls(node.get('pid'), node.get('label'), node.get('productFrame'), comptonProfile)

        xPath.pop()

        return form

class Reactions(suitesModule.Reactions):

    pass
