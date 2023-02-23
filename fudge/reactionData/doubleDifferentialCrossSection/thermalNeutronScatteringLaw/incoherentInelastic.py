# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Thermal neutron scattering law coherent elastic double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import physicalQuantity as physicalQuantityModule
from xData import XYs1d as XYs1dModule
from xData import gridded as griddedModule

from PoPs import IDs as IDsPoPsModule

from fudge import suites as suitesModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from . import base as baseModule
from . import incoherentInelasticMisc as incoherentInelasticMiscModule


def readBool(string):
    # convert xml 'false' / 'true' to False / True.
    # FIXME move to LUPY?
    return string.lower() == 'true'


class Mass( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'mass'


class BoundAtomCrossSection( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'boundAtomCrossSection'


class CoherentAtomCrossSection( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'coherentAtomCrossSection'


class E_critical( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'e_critical'


class E_max( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'e_max'

class T_effective( ancestryModule.AncestryIO ) :
    """
    In incoherent inelastic sections, each scattering atom using the short collision time (SCT)
    approximation needs an effective temperature.
    """

    moniker = 'T_effective'

    def __init__( self, function1d ) :

        ancestryModule.AncestryIO.__init__( self )

        assert isinstance(function1d, XYs1dModule.XYs1d)
        self.__function1d = function1d
        self.__function1d.setAncestor( self )

    @property
    def function1d( self ) :

        return self.__function1d

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__function1d.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.__function1d.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
    
        xPath.append( element.tag )

        xys1d = XYs1dModule.XYs1d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        T_effective1 = cls( xys1d )

        xPath.pop( )
        return( T_effective1 )


class PhononSpectrum(ancestryModule.AncestryIO):
    """
    Stores the phonon spectrum for Gaussian approximation scattering kernel
    """

    moniker = 'phononSpectrum'

    def __init__(self, function1d):
        ancestryModule.AncestryIO.__init__(self)

        assert isinstance(function1d, XYs1dModule.XYs1d)
        self.__xys1d = function1d
        self.__xys1d.setAncestor(self)

    @property
    def XYs1d(self):

        return self.__xys1d

    def convertUnits(self, unitMap):

        self.__xys1d.convertUnits( unitMap )

    def toXML_strList(self, indent='', **kwargs):

        xmlList = ['%s<%s>' % (indent, self.moniker)]
        xmlList += self.__xys1d.toXML_strList(indent + '  ', **kwargs)
        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        assert element[0].tag == XYs1dModule.XYs1d.moniker
        function1d = XYs1dModule.XYs1d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)

        instance = cls(function1d)
        xPath.pop( )
        return instance


class GaussianApproximation(ancestryModule.AncestryIO):
    """
    Used when the scattering kernel can be expressed as a Fourier transform of the material phonon spectrum
    """

    moniker = 'GaussianApproximation'

    def __init__(self, phononSpectrum):
        ancestryModule.AncestryIO.__init__(self)

        assert isinstance(phononSpectrum, PhononSpectrum)
        self.__phononSpectrum = phononSpectrum
        self.__phononSpectrum.setAncestor(self)

    @property
    def phononSpectrum(self):

        return self.__phononSpectrum

    def convertUnits(self, unitMap):

        self.__phononSpectrum.convertUnits( unitMap )

    def toXML_strList(self, indent='', **kwargs):

        xmlList = ['%s<%s>' % (indent, self.moniker)]
        xmlList += self.__phononSpectrum.toXML_strList(indent + '  ', **kwargs)
        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        assert element[0].tag == PhononSpectrum.moniker
        spectrum = PhononSpectrum.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)

        instance = cls(spectrum)
        xPath.pop( )
        return instance


class SCTApproximation(ancestryModule.AncestryIO):
    """
    Used when the scattering kernel can be computed entirely using the short-collision-time approximation.
    """

    moniker = 'SCTApproximation'

    def __init__(self):
        ancestryModule.AncestryIO.__init__(self)

    def convertUnits(self, unitMap):

        pass # nothing to do

    def toXML_strList(self, indent='', **kwargs):

        return [ '%s<%s/>' % ( indent, self.moniker ) ]

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs) :
    
        return cls()


class FreeGasApproximation(ancestryModule.AncestryIO):
    """
    Used when the scattering kernel should be computed using the free gas approximation
    """

    moniker = 'freeGasApproximation'

    def __init__(self):
        ancestryModule.AncestryIO.__init__(self)

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        pass

    def toXML_strList(self, indent='', **kwargs):

        return [ '%s<%s/>' % ( indent, self.moniker ) ]

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs) :

        return cls()


class SelfScatteringKernel(ancestryModule.AncestryIO):
    """
    The 'self' part of the scattering kernel, $S_{self}(\\alpha,\\beta,T)$
    """

    moniker = 'selfScatteringKernel'
    allowedKernels = (griddedModule.Gridded3d, GaussianApproximation, SCTApproximation, FreeGasApproximation)

    def __init__(self, kernel, symmetric=None):

        ancestryModule.AncestryIO.__init__(self)

        assert isinstance(kernel, SelfScatteringKernel.allowedKernels)
        self.__kernel = kernel
        self.__symmetric = symmetric
        self.__kernel.setAncestor(self)

    @property
    def kernel(self):

        return self.__kernel

    @property
    def symmetric(self):

        if self.__symmetric is None:
            if isinstance(self.kernel, (SCTApproximation, FreeGasApproximation)):
                self.__symmetric = True
            elif isinstance(self.kernel, griddedModule.Gridded3d):
                betaAxis, = [axis for axis in self.kernel.axes if axis.label == "beta"]
                self.__symmetric = min(betaAxis.values) >= 0
            else:
                raise NotImplementedError("symmetry in GaussianApproximation")
        return self.__symmetric

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__kernel.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attrs = ""
        if self.symmetric:
            attrs = ' symmetric="true"'
        xmlStringList = ['%s<%s%s>' % (indent, self.moniker, attrs)]
        xmlStringList += self.__kernel.toXML_strList(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        kernel = None
        symmetric = element.get('symmetric','') == 'true'
        for kernelClass in cls.allowedKernels:
            if kernelClass.moniker == element[0].tag:
                kernel = kernelClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        if kernel is None:
            raise NotImplementedError("Self scattering kernel of type %s" % element[0].tag)

        instance = cls(kernel, symmetric)
        xPath.pop( )
        return instance


class DistinctScatteringKernel(ancestryModule.AncestryIO):
    """
    The 'distinct' part of the scattering kernel, $S_{distinct}(\\alpha,\\beta,T)$
    """

    moniker = 'distinctScatteringKernel'
    allowedKernels = (griddedModule.Gridded3d, )

    def __init__(self, kernel):

        ancestryModule.AncestryIO.__init__(self)

        assert isinstance(kernel, DistinctScatteringKernel.allowedKernels)
        self.__kernel = kernel
        self.__kernel.setAncestor(self)

    @property
    def kernel(self):

        return self.__kernel

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__kernel.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.__kernel.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        kernel = None
        for kernelClass in DistinctScatteringKernel.allowedKernels:
            if kernelClass.moniker == element[0].tag:
                kernel = kernelClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        if kernel is None:
            raise NotImplementedError("Distinct scattering kernel of type %s" % element[0].tag)

        instance = cls(kernel)
        xPath.pop( )
        return instance


class ScatteringAtom(ancestryModule.AncestryIO):
    """
    Inelastic incoherent scattering requires a description of each atom that is scattered off.
    """

    moniker = 'scatteringAtom'
    attributes = {'pid': None, 'numberPerMolecule': None, 'primaryScatterer': False}

    def __init__(self, pid, numberPerMolecule, mass, e_max, boundAtomCrossSection, selfScatteringKernel, *,
                 primaryScatterer = False, e_critical = None, coherentAtomCrossSection = None,
                 distinctScatteringKernel = None, T_effective = None):

        ancestryModule.AncestryIO.__init__(self)
        self.pid = pid
        self.primaryScatterer = primaryScatterer
        self.numberPerMolecule = numberPerMolecule

        self.mass = mass
        self.mass.setAncestor(self)

        self.e_max = e_max
        self.e_max.setAncestor(self)

        self.boundAtomCrossSection = boundAtomCrossSection
        self.boundAtomCrossSection.setAncestor(self)

        self.selfScatteringKernel = selfScatteringKernel
        self.selfScatteringKernel.setAncestor(self)

        self.distinctScatteringKernel = distinctScatteringKernel
        if self.distinctScatteringKernel is not None: self.distinctScatteringKernel.setAncestor(self)

        self.e_critical = e_critical
        if self.e_critical is not None: self.e_critical.setAncestor(self)

        self.coherentAtomCrossSection = coherentAtomCrossSection
        if self.coherentAtomCrossSection is not None: self.coherentAtomCrossSection.setAncestor(self)

        self.T_effective = T_effective
        if self.T_effective is not None: self.T_effective.setAncestor(self)

    @property
    def label(self):
        # FIXME required if scatteringAtoms inherits from suite. Replace scatteringAtoms instead?

        return self.pid

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.mass.convertUnits( unitMap )
        self.e_max.convertUnits( unitMap )
        self.boundAtomCrossSection.convertUnits( unitMap )
        self.selfScatteringKernel.convertUnits( unitMap )
        if self.e_critical is not None: self.e_critical.convertUnits( unitMap )
        if self.e_max is not None: self.e_max.convertUnits( unitMap )
        if self.T_effective is not None: self.T_effective.convertUnits( unitMap )

    def shortCollisionTime( self, alpha, beta, T ):
        """
        Equation 7.8 in ENDF manual
        :param alpha:
        :param beta:
        :param T:
        :return:
        """

        Teff = self.T_effective.function1d.evaluate( T )
            # num = math.exp( -( alpha - abs( beta ) )**2 * T / ( 4 * alpha * Teff ) - ( beta + abs( beta ) ) / 2 )
            # reformulate using T/Teff = 1-f, should be more numerically stable
        f = 1 - T / Teff
        num = math.exp( ( f * ( alpha - abs( beta ) )**2 - ( alpha + beta )**2 ) / ( 4 * alpha ) )
        den = math.sqrt( 4 * math.pi * alpha * Teff / T )
        return( num / den )

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        formatVersion = kwargs.get('formatVersion')

        if formatVersion == GNDS_formatVersionModule.version_1_10:
            attrs = ' label="%s" numberPerMolecule="%d"' % (self.pid, self.numberPerMolecule)
            if not self.primaryScatterer:
                functionalForm = {
                    SCTApproximation: 'SCT',
                    FreeGasApproximation: 'free_gas'
                }[type(self.selfScatteringKernel.kernel)]
                attrs += ' functionalForm="%s"' % functionalForm
            xmlStringList = ['%s<%s%s>' % (indent, self.moniker, attrs)]
            xmlStringList += self.mass.toXML_strList(indent=indent2, **kwargs)

            neutronMass = self.rootAncestor.PoPs['n'].mass.float('amu')
            freeAtomCrossSection = self.boundAtomCrossSection.value * self.mass.value**2 / (self.mass.value + neutronMass)
            xmlStringList.append('%s<freeAtomCrossSection value="%f" unit="b"/>' % (indent2, freeAtomCrossSection))
            for child in ('e_critical', 'e_max', 'T_effective'):
                if getattr(self, child) is not None:
                    xmlStringList += getattr(self, child).toXML_strList(indent=indent2, **kwargs)

        else:
            attrs = ' pid="%s" numberPerMolecule="%d"' % (self.pid, self.numberPerMolecule)
            if self.primaryScatterer:
                attrs += ' primaryScatterer="true"'
            xmlStringList = ['%s<%s%s>' % (indent, self.moniker, attrs)]
            for child in ('mass', 'e_critical', 'e_max', 'boundAtomCrossSection', 'coherentAtomCrossSection',
                          'distinctScatteringKernel', 'selfScatteringKernel', 'T_effective') :
                if getattr(self, child) is not None:
                    xmlStringList += getattr(self, child).toXML_strList(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker
        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        formatVersion = kwargs.get('formatVersion')
        attrs = {
            'pid': element.get('pid'),
            'numberPerMolecule': int(element.get('numberPerMolecule')),
            'primaryScatterer': readBool(element.get('primaryScatterer', 'false'))
        }
        for child in element:
            _class = {  Mass.moniker : Mass,
                        E_critical.moniker : E_critical,
                        E_max.moniker : E_max,
                        BoundAtomCrossSection.moniker: BoundAtomCrossSection,
                        CoherentAtomCrossSection.moniker: CoherentAtomCrossSection,
                        SelfScatteringKernel.moniker: SelfScatteringKernel,
                        DistinctScatteringKernel.moniker: DistinctScatteringKernel,
                        T_effective.moniker : T_effective }.get( child.tag, None )
            if _class is None:
                if not (child.tag == 'freeAtomCrossSection' and formatVersion in [GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_4]):
                    print("Warning: encountered unexpected element '%s' in scatteringAtom!" % child.tag)
                continue
            attrs[child.tag] = _class.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        if formatVersion in [GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_4]:
            attrs['pid'] = element.get('label')
            if 'e_max' not in attrs:
                energyUnit = linkData['reactionSuite'].styles.getEvaluatedStyle().projectileEnergyDomain.unit
                attrs['e_max'] = E_max(0, energyUnit)
            if element.get('functionalForm'):
                kernel = {
                    'SCT': SCTApproximation,
                    'free_gas': FreeGasApproximation,
                }[element.get('functionalForm')]()
                attrs[SelfScatteringKernel.moniker] = SelfScatteringKernel(kernel)
            else:   # primary scatterer gets S_alpha_beta
                S_alpha_beta = kwargs.pop('S_alpha_beta')   # read from parent node
                attrs[SelfScatteringKernel.moniker] = SelfScatteringKernel(S_alpha_beta)
                attrs['primaryScatterer'] = True

            massAMU = attrs[Mass.moniker].getValueAs('amu')
            neutronMassAMU = linkData['reactionSuite'].PoPs['n'].mass.float('amu')
            freeAtomCrossSection = float(element.find('freeAtomCrossSection').get('value'))
            attrs[BoundAtomCrossSection.moniker] = BoundAtomCrossSection(
                freeAtomCrossSection * (massAMU + neutronMassAMU) / massAMU**2, 'b')

        SA = ScatteringAtom( **attrs )
        xPath.pop( )
        return SA


class ScatteringAtoms( suitesModule.Suite ) :

    moniker = 'scatteringAtoms'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, ( ScatteringAtom, ) )


class Form(baseModule.Form):

    moniker = 'thermalNeutronScatteringLaw_incoherentInelastic'
    keyName = 'label'

    process = 'thermalNeutronScatteringLaw incoherent-inelastic'
    subformAttributes = ( 'scatteringAtoms', )

    def __init__(self, label, primaryScatterer, pid=IDsPoPsModule.neutron, productFrame=xDataEnumsModule.Frame.lab,
            calculatedAtThermal=False, incoherentApproximation=True):

        baseModule.Form.__init__(self, pid, label, productFrame, (ScatteringAtoms(),))
        self.primaryScatterer = primaryScatterer
        self.calculatedAtThermal = calculatedAtThermal
        self.incoherentApproximation = incoherentApproximation

    def getPrimaryScatterer(self):
        for atom in self.scatteringAtoms:
            if atom.primaryScatterer:
                return atom

    @property
    def domainUnit(self):

        return str(self.getPrimaryScatterer().e_max.unit)

    def energyMaximum(self):

        return float(self.getPrimaryScatterer().e_max)

    @property
    def asymmetric(self):
        """ Is S_alpha_beta asymmetric WRT beta (i.e. tabulated for both positive and negative beta)? """

        if not hasattr(self, '_asymmetric'):
            kernel = self.getPrimaryScatterer().selfScatteringKernel.kernel
            if not isinstance(kernel, griddedModule.Gridded3d):
                return False
            betas = list(kernel.axes[2].values)
            self._asymmetric = min(betas) < 0

        return self._asymmetric

    def processThermalNeutronScatteringLaw(self, style, kwargs):

        from fudge import styles as stylesModule

        temperature = style.temperature.getValueAs(kwargs['temperatureUnit'])

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )

        energyDomain = style.findDerivedFromStyle(stylesModule.Evaluated).projectileEnergyDomain
        energyMax = PQUModule.PQU( energyDomain.max, energyDomain.unit ).getValueAs( kwargs['incidentEnergyUnit'] )

        return incoherentInelasticMiscModule.process(self, style, energyMin, energyMax, temperature, kwargs)

    def temperatures(self, _temperatures):

        primaryKernel = self.getPrimaryScatterer().selfScatteringKernel.kernel
        if isinstance(primaryKernel, griddedModule.Gridded3d):
            _temperatures['incoherent-inelastic'] = [
                primaryKernel.axes[3].unit, primaryKernel.axes[3].values.values ]

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        formatVersion = kwargs.get('formatVersion')

        attributeStr = ' label="%s" pid="%s" productFrame="%s"' % (self.label, self.pid, self.productFrame)
        if formatVersion != GNDS_formatVersionModule.version_1_10:
            attributeStr += ' primaryScatterer="%s"' % self.primaryScatterer
            if self.calculatedAtThermal:
                attributeStr += ' calculatedAtThermal="true"'
            if not self.incoherentApproximation:
                attributeStr += ' incoherentApproximation="false"'
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]

        if formatVersion == GNDS_formatVersionModule.version_1_10:
            # 1.10 had an '<options>' node:
            xmlString.append('%s<options calculatedAtThermal="%s" asymmetric="%s"/>' %
                    (indent2, str(self.calculatedAtThermal).lower(), str(self.asymmetric).lower()))

            xmlString += self.scatteringAtoms.toXML_strList(indent=indent2, **kwargs)

            # 1.10 also stored S_ab separate from the scattering atoms:
            xmlString += ["%s<S_alpha_beta>" % indent2]
            kernel = self.getPrimaryScatterer().selfScatteringKernel.kernel
            if not isinstance(kernel, griddedModule.Gridded3d):
                raise TypeError("GNDS-1.10 only supports tabulated S_alpha_beta, not %s" % type(kernel))
            indent3 = indent2 + kwargs.get('incrementalIndent', '  ')
            xmlString += kernel.toXML_strList(indent=indent3, **kwargs)
            xmlString[-1] += '</S_alpha_beta>'
        else:
            for subform in self.subforms : 
                if subform is not None: xmlString += subform.toXML_strList( indent2, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        formatVersion = kwargs.get('formatVersion')

        attrs = {
            'label': element.get('label'),
            'pid': element.get('pid'),
            'primaryScatterer': element.get('primaryScatterer'),
            'productFrame': element.get('productFrame'),
            'calculatedAtThermal': readBool(element.get('calculatedAtThermal', 'false')),
            'incoherentApproximation': readBool(element.get('incoherentApproximation', 'true')),
        }

        if formatVersion in [GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_4]:
            attrs['calculatedAtThermal'] = readBool(element.find('options').get('calculatedAtThermal'))

            # read S_ab array, to be added to primary scattering atom while parsing scatteringAtoms below
            kwargs['S_alpha_beta'] = griddedModule.Gridded3d.parseNodeUsingClass(
                    element.find('S_alpha_beta').find(griddedModule.Gridded3d.moniker), xPath, linkData, **kwargs)

        incoherentInelastic = Form(**attrs)
        incoherentInelastic.scatteringAtoms.parseNode(
            element.find(ScatteringAtoms.moniker), xPath, linkData, **kwargs
        )
        if formatVersion in [GNDS_formatVersionModule.version_1_10, GNDS_formatVersionModule.version_2_0_LLNL_4]:
            for atom in incoherentInelastic.scatteringAtoms:
                if atom.primaryScatterer:
                    incoherentInelastic.primaryScatterer = atom.pid
                    break

        xPath.pop( )
        return incoherentInelastic
