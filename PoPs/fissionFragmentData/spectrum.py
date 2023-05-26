# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module defines the spectrum.Suite and spectrum.XYs1d classes, used to store outgoing energy
spectra for beta-delayed neutrons from spontaneous fission."""

from xData import XYs1d as XYs1dModule

from .. import suite as suiteModule

class XYs1d( XYs1dModule.XYs1d ) :
    """Stores an outgoing delayed neutron energy spectrum."""

    pass

class Suite( suiteModule.SortedSuite ) :
    """Stores one or more delayed neutron spectra, each corresponding to a unique data **style**."""

    moniker = 'spectrum'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( XYs1d, ) )
