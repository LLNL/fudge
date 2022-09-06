# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.distributions.energyAngular module.
"""

from xData import enums as xDataEnumsModule
from fudge.productData.distributions import KalbachMann as KalbachMannModule

def toACE( self, label, offset, weight, **kwargs ) :

    header = [ 0, 44, offset + len( weight ) + 4 ] + weight
    fSubformData = self.fSubform.data

    interpolation = fSubformData.interpolation
    if interpolation != xDataEnumsModule.Interpolation.linlin: raise ValueError( 'interpolation = "%s" not supported' % interpolation )

    NE = len( fSubformData )
    e_ins = []
    Ls = []
    epData = []
    offset += len( header ) + 1 + 1 + 2 * NE + 1        # header length plus NR, NE, Es, Ls, (1-based).
    for i1, POfEp in enumerate( fSubformData ) :
        e_ins.append( POfEp.outerDomainValue )
        Ls.append( offset + len( epData ) )

        eps = POfEp.xs.values.values
        pdf = POfEp.pdf.values.values
        cdf = POfEp.cdf.values.values

        function = self.rSubform.data[i1]
        Rs = [ function.evaluate( ep ) for ep in eps ]

        function = self.aSubform.data[i1]
        As = [ function.evaluate( ep ) for ep in eps ]

        interpolation = {xDataEnumsModule.Interpolation.flat: 1, xDataEnumsModule.Interpolation.linlin: 2}[POfEp.interpolation]
        epData += [ interpolation, len( eps ) ] + eps + pdf + cdf + Rs + As

    return( header + [ 0, NE ] + e_ins + Ls + epData )

KalbachMannModule.Form.toACE = toACE
