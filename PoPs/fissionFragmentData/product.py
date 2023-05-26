# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# FIXME CMM: we should merge this with the PoPs/decay/product.py module!
# Maybe as part of broader re-org of PoPs in future GNDS version.
# Also, this is currently unused although it should be used at least for Cf258 in the decay sub-library.

"""This module defines the PoPs.fissionFragmentData.Product class, used to store
spontaneous fission products including multiplicities and outgoing energy spectra."""

from LUPY import ancestry as ancestryModule

from .. import IDs as IDsPoPsModule

class Product(ancestryModule.AncestryIO):

    moniker = 'product'

    def __init__( self, pid, productFrame ) :
        """Construct a Product instance for the given particle and product frame.

        :param pid: string identifying the PoPs particle, e.g. 'n'
        :param productFrame: 
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__pid = pid
        self.__productFrame = productFrame

        self.__multiplicity = multiplicity.Suite()
        self.__multiplicity.setAncestor(self)

        self.__spectrum = spectrumModule.Suite()
        self.__spectrum.setAncestor(self)

    @property
    def pid(self):
        """Return the PoPs particle id or 'pid' for this product."""

        return self.__pid

    @property
    def productFrame(self):
        """Return the reference frame (typically 'lab') for the outgoing spectrum."""

        return self.__productFrame

    @property
    def multiplicity(self):
        """Return the product multiplicity (multiplicity.Suite instance)."""

        return self.__multiplicity

    @property
    def spectrum(self):
        """Return the product outgoing energy spectrum (spectrum.Suite instance)."""

        return self.__spectrum

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:    The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:    A keyword list.

        :return:          List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLStringList = [ '%s<%s pid="%s" productFrame="%s">' % (indent, self.moniker, self.__pid, self.__productFrame) ]
        XMLStringList += self.__multiplicity.toXML_strList(indent2, **kwargs)
        XMLStringList += self.__spectrum.toXML_strList(indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker
        return XMLStringList

    # FIXME need xml parsing logic! Not implemented since this class isn't used yet.
