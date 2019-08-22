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

"""This module contains the 'branching3d' form for distribution."""

from xData import standards as standardsModule

from PoPs import IDs as IDsPoPsModule

from fudge.processing.deterministic import transferMatrices as transferMatricesModule
from fudge.processing import group as groupModule

from .. import multiplicity as multiplicityModule

from .. import energyDeposition as energyDepositionModule
from .. import momentumDeposition as momentumDepositionModule

from . import base as baseModule

__metaclass__ = type

class form( baseModule.form ) :

    moniker = 'branching3d'
    subformAttributes = [ 'pids' ]
    ancestryMembers = ( 'pids', )

    def __init__( self, label, _pids, productFrame = standardsModule.frames.labToken ) :

        baseModule.form.__init__( self, label, productFrame, [] )

        if( not( isinstance( _pids, pids ) ) ) : raise TypeError( 'Invalid pids' )
        self.__pids = _pids

    @property
    def pids( self ) :

        return( self.__pids )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = kwargs['momentumDepositionUnit']
        domainMin = kwargs['EMin']
        domainMax = kwargs['EMax']

        Q = kwargs['outputChannel'].Q[0]
        QValue = Q.evaluate( domainMin )

        energyData = [ [ [ domainMin, QValue ], [ domainMax, QValue ] ] ]       # Gives all the energy to the gammas and non to the residual.
        momentumData = [ [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ] ]           # Assumes isotropic scattering.

        return( [ energyData, momentumData ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        def processDecayModesSetup( product, probability, gammas ) :

            initialEnergy = product.nucleus.energy[0].value
            for decayMode in product.decayData.decayModes :
                branchingRatio = decayMode.probability[0].value * probability

                decayPath = decayMode.decayPath[0]
                residual = decayPath.products[0].pid
                if( residual == IDsPoPsModule.photon ) : residual = decayPath.products[1].pid
                residual = PoPs[residual]

                gammaEnergy = initialEnergy - residual.nucleus.energy[0].value
                if( gammaEnergy not in gammas ) : gammas[gammaEnergy] = 0.0
                gammas[gammaEnergy] += branchingRatio

                processDecayModesSetup( residual, branchingRatio, gammas )

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']

        if( verbosity > 2 ) : print( '%sGrouping %s' % ( indent, self.moniker ) )

        crossSection = tempInfo['crossSection']
        domainMin = crossSection.domainMin
        domainMax = crossSection.domainMax

        PoPs = tempInfo['reactionSuite'].PoPs
        nuclide = PoPs[self.pids.initial]

        gammas = {}
        processDecayModesSetup( nuclide, 1.0, gammas )

        TM_1s = []
        TM_Es = []
        multipliticyAxes = multiplicityModule.defaultAxes( tempInfo['incidentEnergyUnit'] )
        for gammaEnergy in gammas :
            branchingRatio = gammas[gammaEnergy]
            multiplicity = multiplicityModule.XYs1d( data = [ [ domainMin, branchingRatio ], [ domainMax, branchingRatio ] ], axes = multipliticyAxes )

            TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( style, tempInfo, gammaEnergy, crossSection, None, 
                    multiplicity = multiplicity, comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            TM_1s.append( TM_1 )
            TM_Es.append( TM_E )
        TM_1 = transferMatricesModule.addTMs( TM_1s )
        TM_E = transferMatricesModule.addTMs( TM_Es )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame

        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        xmlString += self.pids.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker

        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        _pids = pids.parseXMLNode( element[0], xPath, linkData )
        _form = form( element.get( 'label' ), _pids, productFrame = element.get( 'productFrame' ) )

        xPath.pop( )
        return( _form )

class pids( baseModule.subform ) :

    moniker = 'pids'

    def __init__( self, initial, final ) :

        self.__initial = initial
        self.__final = final

    @property
    def initial( self ) :

        return( self.__initial )

    @property
    def final( self ) :

        return( self.__final )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [ '%s<%s initial="%s" final="%s"/>' % ( indent, self.moniker, self.initial, self.final ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        _pids = pids( element.get( 'initial' ), element.get( 'final' ) )
        xPath.pop( )

        return( _pids )
