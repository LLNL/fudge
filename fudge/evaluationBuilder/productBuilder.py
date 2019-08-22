# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

from . import base
from fudge.gnd import product, productData

class productBuilder:
    def __init__(self, particle, multiplicity=None, distributionForm=None, distributionData=None):
        self.product = product.product( particle, "product", label = particle.name )
        if multiplicity is not None:
            self.addMultiplicity( multiplicity )
        if distributionForm is not None:
            self.addDistribution( distributionForm, distributionData )

    def addMultiplicity(self, mult):
        """Add product multiplicity. 'mult' may be an integer (for constant multiplicity)
        or a list of energy-dependent multiplicities of form [ [energies], [multiplicities] ].
        
        Currently, only one linear-linear interpolation region is supported for energy-dependent multiplicities."""

        if type(mult) is int:
            self.product.setMultiplicity( productData.multiplicity.constant(mult) )
        elif type(mult) is list:
            axes_ = productData.multiplicity.pointwise.defaultAxes()
            self.product.setMultiplicity( productData.multiplicity.pointwise( axes_, mult, base.defaultAccuracy,
                dataForm="XYs" ) )
        else:
            raise Exception("Can't understand multiplicity of form '%s'!" % mult)

    def addDistribution(self, dataForm, data=None, frame=None):
        from . import distributionBuilder
        db = distributionBuilder( dataForm, data, frame )
        self.product.distributions.addComponent( db.component )
        self.product.distributions.setNativeData( db.component.name )
