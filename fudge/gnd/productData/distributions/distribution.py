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

"""Distribution class."""

from fudge.core.utilities import brb

from fudge.gnd import abstractClasses as abstractClassesModule

from . import angular as angularModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import KalbachMann as KalbachMannModule
from . import angularEnergy as angularEnergyModule
from . import uncorrelated as uncorrelatedModule
from . import Legendre as LegendreModule
from . import photonScattering as photonScatteringModule
from . import reference as referenceModule
from . import multiGroup as multiGroupModule
from . import unspecified as unspecifiedModule
from . import unknown as unknownModule

# probably missing stuff from photonScattering.py.

__metaclass__ = type

class component( abstractClassesModule.component ) :

    moniker = 'distribution'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( angularModule.form, angularModule.twoBodyForm,
                angularModule.CoulombExpansionForm, angularModule.NuclearPlusCoulombInterferenceForm,
                energyModule.form, KalbachMannModule.form,
                energyAngularModule.form, angularEnergyModule.form, angularEnergyModule.form,
                angularEnergyModule.LLNLAngularEnergyForm,
                uncorrelatedModule.form, LegendreModule.form, referenceModule.form,
                photonScatteringModule.incoherent.form, photonScatteringModule.coherent.form, 
                multiGroupModule.form, 
                unspecifiedModule.form, unknownModule.form ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

        return( self.evaluated.getSpectrumAtEnergy( energy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        form = style.findFormMatchingDerivedStyles( self )
        if( form is None ) : raise Exception( 'No matching style' )
#        print form.moniker
        return( form.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. Need ability to find specific distribution component
        """
        if attribute is not None:
            for entity in self:
                if entity.moniker == entityName and getattr(entity,attribute) == value:
                    return entity
        else:
            for entity in self:
                if entity.moniker == entityName:
                    return entity
        raise KeyError( "Can't find entity '%s' in distributions!" % entityName )

    def hasData( self ) :
        """
        Returns False if self's only has unspecified or unknown form; otherwise, returns True.
        """

        for form in self :
            if( not( isinstance( form, ( unspecifiedModule.form, unknownModule.form ) ) ) ) : return( True )
        return( False )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.getEvaluated( ).toPointwise_withLinearXYs( accuracy, lowerEps, upperEps ) )
