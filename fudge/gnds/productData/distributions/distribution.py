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

from PoPs import IDs as IDsPoPsModule

from fudge.gnds import abstractClasses as abstractClassesModule
from xData import standards as standardsModule

from . import angular as angularModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import energyAngularMC as energyAngularMCModule
from . import angularEnergyMC as angularEnergyMCModule
from . import KalbachMann as KalbachMannModule
from . import angularEnergy as angularEnergyModule
from . import uncorrelated as uncorrelatedModule
from . import Legendre as LegendreModule
from . import photonScattering as photonScatteringModule
from . import reference as referenceModule
from . import multiGroup as multiGroupModule
from . import unspecified as unspecifiedModule
from . import branching3d as branching3dModule

# probably missing stuff from photonScattering.py.

__metaclass__ = type

class component( abstractClassesModule.component ) :

    moniker = 'distribution'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( angularModule.form, angularModule.twoBodyForm,
                KalbachMannModule.form,
                energyAngularModule.form, energyAngularMCModule.form,
                angularEnergyModule.form, angularEnergyMCModule.form,
                angularEnergyModule.LLNLAngularEnergyForm,
                uncorrelatedModule.form, LegendreModule.form, referenceModule.form,
                referenceModule.CoulombPlusNuclearElastic,
                photonScatteringModule.coherentPhotonScattering.form, photonScatteringModule.incoherentPhotonScattering.form, 
                multiGroupModule.form, unspecifiedModule.form, branching3dModule.form ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

        return( self.evaluated.getSpectrumAtEnergy( energy ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        form = style.findFormMatchingDerivedStyle( self )
        if( form is None ) : raise Exception( 'No matching style' )
        return( form.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def check( self, info ):
        """check all distribution forms"""
        from fudge.gnds import warning
        warnings = []

        for form in self:

            if info['isTwoBody']:
                if( form.productFrame != standardsModule.frames.centerOfMassToken ) :
                    warnings.append( warning.wrong2BodyFrame( form ) )

                if form.moniker not in (angularModule.twoBodyForm.moniker,
                                        referenceModule.form.moniker,
                                        referenceModule.CoulombPlusNuclearElastic.moniker,
                                        unspecifiedModule.form.moniker):
                    warnings.append( warning.wrongDistributionComponent( form.moniker, '2-body' ) )
            else:
                if form.moniker in (angularModule.twoBodyForm.moniker,
                                    angularModule.form.moniker,
                                    energyModule.form.moniker):
                    warnings.append( warning.wrongDistributionComponent( form.moniker, 'N-body' ) )

            def checkSubform( subform, contextMessage ):
                distributionErrors = []
                if hasattr(subform, 'domainMin') and (subform.domainMin, subform.domainMax) != info['crossSectionDomain']:
                    domain = (subform.domainMin, subform.domainMax)
                    # For gamma products, domainMin should be >= cross section start, upper bounds should match.
                    if( self.ancestor.id == IDsPoPsModule.photon ) :
                        startRatio = subform.domainMin / info['crossSectionDomain'][0]
                        endRatio = subform.domainMax / info['crossSectionDomain'][1]
                        if (startRatio < 1-standardsModule.floats.epsilon or endRatio < 1-standardsModule.floats.epsilon
                                or endRatio > 1+standardsModule.floats.epsilon):
                            distributionErrors.append( warning.domain_mismatch(
                                    *(domain + info['crossSectionDomain']), obj=subform ) )
                    # For all other products, check lower and upper edges: only warn if they disagree by > eps
                    else:
                        for e1,e2 in zip(domain, info['crossSectionDomain']):
                            ratio = e1 / e2
                            if (ratio < 1-standardsModule.floats.epsilon or ratio > 1+standardsModule.floats.epsilon):
                                distributionErrors.append( warning.domain_mismatch(
                                        *(domain + info['crossSectionDomain']), obj=subform ) )
                                break
                    if not hasattr(subform,'check'):
                        distributionErrors.append( warning.NotImplemented(subform.moniker, subform ) )
                        if info['failOnException']:
                            raise NotImplementedError("Checking distribution form '%s'" % subform.moniker)
                else:
                    distributionErrors += subform.check( info )
                if distributionErrors:
                    warnings.append( warning.context( contextMessage + " - %s:" % subform.moniker, distributionErrors) )

            if isinstance(form, uncorrelatedModule.form):
                for subformName in ('angularSubform','energySubform'):
                    subform = getattr(form, subformName ).data
                    checkSubform( subform, 'uncorrelated - ' + subformName.replace('Subform','') )

            elif isinstance(form, KalbachMannModule.form):
                checkSubform( form, form.moniker )

            else:
                for subform in form.subforms:
                    checkSubform( subform, form.moniker )

        return warnings


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
        return abstractClassesModule.component.findEntity( self, entityName, attribute, value )

    def hasData( self ) :
        """
        Returns False if self's only has unspecified form; otherwise, returns True.
        """

        for form in self :
            if( not( isinstance( form, unspecifiedModule.form ) ) ) : return( True )
        return( False )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.getEvaluated( ).toPointwise_withLinearXYs( **kwargs ) )
