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

"""
This module contains all the classes for the various types of ENDL I data.
"""

import copy

import endl2dmathClasses, endl2dmathmisc, endl3dmathClasses, endl3dmathmisc, endl4dmathClasses, endl4dmathmisc
from fudge.core.math import fudge2dGrouping
import endlmath
from fudge.core.utilities import fudgeExceptions
from fudge.vis.gnuplot import plotbase
import endlIClassesParameters
import endlNd
import endlmisc
import endlParameters
from endl2 import yoToZA

normCheckTolerance = 1e-5
fixThresholdMode_None = endlIClassesParameters.fixThresholdMode_None
fixThresholdMode_RaiseOnly = endlIClassesParameters.fixThresholdMode_RaiseOnly
fixThresholdMode_All = endlIClassesParameters.fixThresholdMode_All
fixThreshold_deltaFunctionEpsilon = 1e-4

def endlAddIObject( f, yo, C, I, S, h, points, bdflsFile = None ) :
    "For internal use only."

    if   ( I ==   0 ) : return endlI0( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==   7 ) : return endlI7( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==   9 ) : return endlI9( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  10 ) : return endlI10( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  11 ) : return endlI11( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  12 ) : return endlI12( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  13 ) : return endlI13( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  20 ) : return endlI20( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  80 ) : return endlI80( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  89 ) : return endlI89( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  90 ) : return endlI90( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  91 ) : return endlI91( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  92 ) : return endlI92( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I == 941 ) : return endlI941( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I == 942 ) : return endlI942( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==   1 ) : return endlI1( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  21 ) : return endlI21( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  22 ) : return endlI22( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  81 ) : return endlI81( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==  84 ) : return endlI84( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==   3 ) : return endlI3( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    elif ( I ==   4 ) : return endlI4( f, yo, C, I, S, h, points, bdflsFile = bdflsFile )
    else : raise Exception( "\nError in endlIClasses.endlAddIObject: I = %d not supported" % I )

class endlI0( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 0 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI0 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = 'crossSection'
        endlNd.endlNd.__init__( self, f, 0, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];cross_section[barn]'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = True, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, 0., EMin = EMin, fixThresholdMode = fixThresholdMode, 
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def heat( self, T, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False, doNotThin = True, EMin = 1e-11, 
        heatBelowThreshold = True, removeClosePoints = True, heatAllEDomain = True ) :
        """Calls crossSectionAdjustForHeatedTarget with self's data. Returns an endl2dmath
    instance. See module crossSectionAdjustForHeatedTarget for more details."""

        TData = self.getTemperature( )
        if( TData > T ) : raise Exception( "\nError in endlI0.heat: temperature  = %e must be greater than data's temperature = %e" % ( T, TData ) )
        unheated = self.copyData( )
        unheated.trim( )
        dT = T - TData
        if( ( len( self ) > 1 ) and ( dT > 0. ) ) :
            massRatio = self.bdflsFile.mass( self.ZA ) / self.bdflsFile.mass( self.yi )
            if( lowerlimit == None ) :
                lowerlimit = "oneOverV"
                if( self.C == 10 ) : lowerlimit = "constant"
                elif ( unheated.xMin( ) > 1e-8 ) : lowerlimit = "threshold"
            if( upperlimit == None ) : upperlimit = 'constant'
            if( not heatBelowThreshold ) : EMin = max( unheated.data[0][0], EMin )
            from crossSectionAdjustForHeatedTarget import heat as heatModule
            d  = heatModule.crossSectionAdjustForHeatedTarget( massRatio, dT, EMin, unheated.data, lowerlimit = lowerlimit, \
                upperlimit = upperlimit, interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                heatAllEDomain = heatAllEDomain )
            heated = endl2dmathClasses.endl2dmath( d )
            if( removeClosePoints ) : heated.removeClosePoints( )
        else :
            heated = unheated
        I0 = endlI0( None, self.yo, self.C, self.I, self.S, self.h, heated, bdflsFile = self.bdflsFile )
        if( dT < 0 ) : 
            I0.setTemperature( TData )
        else :
            I0.setTemperature( T )
        return( I0 )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "cross section (barn)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

    def toZAsFrame( self, newProjectileMass, newTargetMass, halflife, bdflsFile, ELevel = 0. ) :

        newI0, r = toZAsFrameMisc( endlI0, self, newProjectileMass, newTargetMass, ELevel, halflife, bdflsFile, True )
        return( [ newI0 ] )

class endlI7( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 7 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI7 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = 'multiplicity'
        endlNd.endlNd.__init__( self, f, 7, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];multiplicity'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = True, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode, 
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "nu_bar", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI9( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 9 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI9 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = 'multiplicity'
        endlNd.endlNd.__init__( self, f, 9, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];multiplicity'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = True, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "multiplicity", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )


    def toZAsFrame( self, newProjectileMass, newTargetMass, halflife, bdflsFile, ELevel = 0. ) :

        newI9, r = toZAsFrameMisc( endlI9, self, newProjectileMass, newTargetMass, ELevel, halflife, bdflsFile, True )
        return( [ newI9 ] )

class endlI10( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 10 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI10 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        from fudge import gnd

        self.name = gnd.productData.energyDeposition.component.moniker
        endlNd.endlNd.__init__( self, f, 10, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];energy_deposition[MeV]'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = True, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "E to yo (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI11( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 11 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI11 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        from fudge import gnd

        self.name = gnd.productData.energyDeposition.component.moniker
        endlNd.endlNd.__init__( self, f, 11, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];energy_deposition[MeV]'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = True, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "E to res. (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI12( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 12, Q(E), data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI12 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = 'energyDependentQ'
        endlNd.endlNd.__init__( self, f, 12, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];Q[MeV]'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "Q (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI13( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 13 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI13 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        from fudge import gnd

        self.name = gnd.productData.momentumDeposition.component.moniker
        endlNd.endlNd.__init__( self, f, 13, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];momentum_deposition[MeV/c]'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor2d to fix threshold."""

        fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, None, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        return( getThresholdsForChecker2d( self ) )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "momentum (MeV/c)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI80( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 80 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI80 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 80, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "T ave'd rates (barn * cm/sh)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI89( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 89 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI89 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 89, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "multiplicity tn", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI90( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 90 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI90 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 90, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "T ave'd E to yo (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI91( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 91 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI91 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 91, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "T ave'd E to res. (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI92( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 92 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI92 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 92, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "T (MeV)", yLabel = "T ave'd total E (MeV)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI941( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 941 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI941 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        from fudge import gnd

        self.name = 'pointwise'
        endlNd.endlNd.__init__( self, f, 941, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];form_factor'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = True, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "form factor", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

    def toInterpolation( self, interpolation, accuracy, diSectionMax = 3 ) :

        data = self
        if( ( interpolation == 0 ) and ( self.interpolation == 3 ) and ( len( self ) > 0 ) ) :
            x, y = self.data[0]
            if( ( x == 0 ) or ( y == 0 ) ) :
                data = endlI941( None, self.yo, self.C, self.I, self.S, self.h, self.data )
                data.data = data.data[1:]
        newData = endl2dmathClasses.endl2dmath.toInterpolation( data, interpolation, accuracy, diSectionMax = diSectionMax )
        if( data != self ) : newData.data.insert( 0, [ x, y ] )
        return( newData )

class endlI942( endlNd.endlNd, endl2dmathClasses.endl2dmath ) :
    """This class is for the I = 942 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI942 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        from fudge import gnd

        self.name = 'pointwise'
        endlNd.endlNd.__init__( self, f, 942, yo, C, I, S, h, points, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];scattering_function'

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = True, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "scattering function", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

    def toInterpolation( self, interpolation, accuracy, diSectionMax = 3 ) :

        data = self
        if( ( interpolation == 0 ) and ( self.interpolation == 3 ) and ( len( self ) > 0 ) ) :
            x, y = self.data[0]
            if( ( x == 0 ) or ( y == 0 ) ) :
                data = endlI942( None, self.yo, self.C, self.I, self.S, self.h, self.data )
                data.data = data.data[1:]
        newData = endl2dmathClasses.endl2dmath.toInterpolation( data, interpolation, accuracy, diSectionMax = diSectionMax )
        if( data != self ) : newData.data.insert( 0, [ x, y ] )
        return( newData )

class endlI951( endlNd.endlNd, endl2dmathClasses.endl2dmath ) : # Special I = 951 for NADS (i.e., sigma * v)
    """This class is for the I = 951 data. Special I value for NADS."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI951 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl2dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 951, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Calls the endl2dmathmisc.check2dData function and returns a list of endlCheckerObject instances.
        See endl2dmathmisc.check2dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        n, dummy, messages2d = endl2dmathmisc.check2dData( self.data, allowZeroX = allowZeroE, positiveY = False, printWarning = printWarning, \
            printErrors = printErrors, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) : ErrMsgs.append( endlmisc.endlCheckerObject( data = self, message = messages2d ) )
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "cross section * v (barn * cm/sec)", interpolation = 0 ) :

        endl2dmathClasses.endl2dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, interpolation = interpolation )

class endlI1( endlNd.endlNd, endl3dmathClasses.endl3dmath ) :
    """This class is for the I = 1 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI1 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        from fudge import gnd

        self.name = 'pointwise'
        endlNd.endlNd.__init__( self, f, 1, yo, C, I, S, h, points, i2 = 2, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];mu;P(mu|energy_in)'

    def check( self, normTolerance = 1e-5, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """
        Checks to see that the data is consistance with I = 1 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors.
        """

        normTolerance = max( normTolerance, normCheckTolerance )
        ErrMsgs = []
        messages = []
        d2 = endlmath.ZSum( self.data )
        for E, y in d2.data :
            if( abs( 1. - y ) > normTolerance ) : messages.append( 'endlI1.check: bad normalize = %.8e for E = %e' % ( y, E ) )
        messages += endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowZeroX = allowZeroE, allowNegativeY = True, allowZeroY = True, \
            positiveZ = False, printWarning = False, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        for E, muP in self.data :
            if( muP[0][0] < -1. ) : messages.append( 'endlI1.check: mu = %.8e < -1 for E = %e' % ( muP[0][0], E ) )
            if( muP[-1][0] > 1. ) : messages.append( 'endlI1.check: mu = %.8e > 1 for E = %e' % ( muP[-1][0], E ) )
        messages += self.checkMus( True )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def checkMus( self, internal = False ) :
        """
        This method checks that the mu values are reasonable. For two-body reactions, the mu domain must range from -1 to 1.
        For other data, the first mu value must starts at -1 if above threshold and must be greater than -1 if at (near) threshold.
        """

        import math

        messages, yo = [], self.yo
        if( yo > 9 ) : yo -= 10

        if( ( self.C == 10 ) or ( self.S == 1 ) ) :        # Two-body.
            for e_in, muPs in self.data :
                if( ( muPs[0][0] != -1 ) or ( muPs[-1][0] != 1 ) ) :
                    messages.append( 'endlI1.checkMus: for two-body, mu domain invalid: domain is %.5f for %.5f for E = %e' % ( muPs[0][0], muPs[-1][0], e_in ) )
        else :
            if( yo == 7 ) : return( [] )

            projectileMass, targetMass, productMass = self.bdflsFile.mass( self.yi ), self.bdflsFile.mass( self.ZA ), self.bdflsFile.mass( yo )
            compoundMass = projectileMass + targetMass
            Q = self.getQ( )

            for i1, e_muPs in enumerate( self.data ) :
                e_in, muPs = e_muPs
                uMax2 = 2 * ( e_in * targetMass / compoundMass + Q ) / productMass  # maximum speed of product squared in COM frame.
                vCOM2 = 2 * projectileMass * e_in / compoundMass**2                 # speed of COM squared.
                atThreshold = ( i1 == 0 ) and ( Q < 0 )
                if( atThreshold ) :
                    if( muPs[0][0] == -1 ) : messages.append( 'endlI1.checkMus: need delta function in mu, mu = -1 for E = %e' % ( e_in ) )
                else :
                    if( uMax2 >= vCOM2 ) :
                        if( muPs[0][0] > -1 ) : messages.append( 'endlI1.checkMus: mu range too short, mu = %.5e > -1 for E = %e' % ( muPs[0][0], e_in ) )
                    else :
                        if( muPs[0][0] == -1 ) : messages.append( 'endlI1.checkMus: need delta function, mu = %.8e == -1 for E = %e' % ( muPs[0][0], e_in ) )
        if( not( internal ) ) :
            if( len( messages ) > 0 ) : messages = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( messages )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor3d to fix threshold."""

        fixThresholdFor3d( self, thresholdCrossSectionIsZero, self.I, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        if( len( self ) > 1 ) : return [ self.data[0][0], self.data[1][0] ]
        if( len( self ) > 0 ) : return [ self.data[0][0] ]
        return []

    def getEData( self, E ) :
        """Returns an endl2dmath object for the data at incident energy E. If E is outside the
        domain of the data, then None is returned.  If the requested E is not a point in the data 
        then linear interpolation is performed. """

        EmuP = self.data
        if( EmuP == None ) : return( None )
        nE = len( EmuP )
        if( nE == 0 ) : return( None )
        if( ( E < EmuP[0][0] ) or ( E > EmuP[-1][0] ) ) : return( None )
        EPrior = None
        muPPrior = None
        for E_, muP_ in self.data :
            if( E_ >= E ) :
                muP = []
                if( E_ == E ) :
                    for mu, P in muP_ : muP.append( [ mu, P ] )
                    return( endl2dmathClasses.endl2dmath( muP ) )
                return( ( ( E - EPrior ) * endl2dmathClasses.endl2dmath( muP_ ) + ( E_ - E ) * endl2dmathClasses.endl2dmath( muPPrior ) ) / ( E_ - EPrior ) )
            EPrior = E_
            muPPrior = muP_

    def normalize( self ) :
        """Normalizes the data so that the integral P(E, mu ) dmu = 1."""

        d2 = endlmath.ZSum( self.data )
        i = 0
        for x_zy in self.data :
            s = d2.data[i][1]
            if( s != 0. ) :
                for yz in x_zy[1] : yz[1] = yz[1] / s
            i += 1

    def runningMuSum( self ) :
        """Does a running sum of the mu data, returning an endl3dmath object of the results."""

        d3 = endlmath.runningZSum( self.data, xLabel = self.xLabel, yLabel = self.yLabel, zLabel = "running Sum of P( E, mu ) vs mu" )
        return d3

    def setEData( self, E, muP ) :
        """This method adds muP at E. If E is in self than its distribution is over-written; otherwise, a new E is added."""

        i = 0
        muP_ = endl2dmathmisc.get2dmathData( muP, "setEData", "what ever" ) 
        muP = []
        for mu, P in muP_ : muP.append( [ mu, P ] )
        for E_, muP_ in self.data :
            if( E <= E_ ) : break
            i += 1
        if( i < len( self ) ) :
            if( self.data[i][0] == E ) :
                self.data[i][1] = muP
            else :
                self.data.insert( i, [ E, muP ] )
        else :
            self.data.append( [ E, muP ] )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "mu", zLabel = "pdf(E,mu) vs mu", interpolation = 0  ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel
            , interpolation = interpolation )

    def toZAsFrame( self, newProjectileMass, newTargetMass, halflife, bdflsFile, ELevel = 0. ) :

        newI1, r = toZAsFrameMisc( endlI1, self, newProjectileMass, newTargetMass, ELevel, halflife, bdflsFile, True )
        return( [ newI1 ] )

class endlI21( endlNd.endlNd, endl3dmathClasses.endl3dmath ) :
    """This class is for the I = 21 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI21 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        from fudge import gnd

        self.name = ''
        self.form = 'pointwise'
        endlNd.endlNd.__init__( self, f, 21, yo, C, I, S, h, points, i2 = 2, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 21 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        messages = endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowZeroX = allowZeroE, allowNegativeY = False, positiveZ = False, \
            printWarning = False, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor3d to fix threshold."""

        fixThresholdFor3d( self, thresholdCrossSectionIsZero, self.I, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "E' (MeV)", zLabel = "pdf(E,E') vs E' (1/MeV)", interpolation = 0 ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

class endlI22( endlNd.endlNd, endl3dmathClasses.endl3dmath ) :
    """This class is for the I = 22 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI22 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 22, yo, C, I, S, h, points, i2 = 2, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 22 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        messages = endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowZeroX = allowZeroE, allowNegativeY = False, positiveZ = False, \
            printWarning = False, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        for E, muP in self.data :
            if( muP[0][0] < 0. ) : messages.append( 'endlI22.check: 1 - mu = %.8e < 0 for E = %e' % ( muP[0][0], E ) )
            if( muP[-1][0] > 2. ) : messages.append( 'endlI22.check: 1 - mu = %.8e > 2 for E = %e' % ( muP[0][0], E ) )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor3d to fix threshold."""

        fixThresholdFor3d( self, thresholdCrossSectionIsZero, self.I, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.xMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "x = 1 - mu", zLabel = "pdf(E,x) vs x", interpolation = 0 ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel, 
            interpolation = interpolation )

class endlI81( endlNd.endlNd, endl3dmathClasses.endl3dmath ) :
    """This class is for the I = 81 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI81 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 81, yo, C, I, S, h, points, i2 = 2, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 81 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        messages = endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowZeroX = False, allowNegativeY = False, allowZeroY = allowZeroE, \
            positiveZ = False, printWarning = False, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "T (MeV)", yLabel = "E (MeV)", zLabel = "cross section (T,E)", interpolation = 0 ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

class endlI84( endlNd.endlNd, endl3dmathClasses.endl3dmath ) :
    """This class is for the I = 84 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI84 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 84, yo, C, I, S, h, points, i2 = 2, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 84 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        messages = endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowZeroX = False, allowNegativeY = False, positiveZ = False, \
            printWarning = False, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "E' (MeV)", zLabel = "f(T,E')", interpolation = 0 ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

class endlI952( endlNd.endlNd, endl3dmathClasses.endl3dmath ) : # Special I = 952 for NADS (i.e., total P(E,E') )
    """This class is for the I = 952 data. Special I value for NADS."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI952 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl3dmath class data structure."""

        self.name = ''
        endlNd.endlNd.__init__( self, f, 952, yo, C, I, S, h, points, bdflsFile = bdflsFile )

    def check( self, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 952 data and returns a list of 
        endlCheckerObject instances. Also, calls the endl3dmathmisc.check3dData function. 
        See endl3dmathmisc.check3dData for meaning of printWarning and printErrors."""

        ErrMsgs = []
        messages = endl3dmathmisc.check3dData( self.data, allowNegativeX = False, allowNegativeY = False, positiveZ = False, printWarning = False, \
            printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def set( self, data, checkDataType = 0, xLabel = "E (MeV)", yLabel = "E' (MeV)", zLabel = "P(E,E')", interpolation = 0 ) :

        endl3dmathClasses.endl3dmath.__init__( self, data, checkDataType = checkDataType, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

class endlI3( endlNd.endlNd, endl4dmathClasses.endl4dmath ) :
    """This class is for the I = 3 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI3 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl4dmath class data structure."""

        from fudge import gnd

        self.name = 'pointwise'
        endlNd.endlNd.__init__( self, f, 3, yo, C, I, S, h, points, i2 = 2, i3 = 3, bdflsFile = bdflsFile )
        self.xInterpolation = 'linear,linear,unitbase:x:linear,unitbase:y:linear'
        self.variablesUnits = 'energy_in[MeV];mu;energy_out[MeV];pdf(energy_out|mu,energy_in)'

    def check( self, normTolerance = 1e-5, printWarning = True, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 3 data and returns a list of
        endlCheckerObject instances. Also, calls the endl4dmathmisc.check4dData function.
        See endl4dmathmisc.check4dData for meaning of printWarning and printErrors."""

        normTolerance = max( normTolerance, normCheckTolerance )
        ErrMsgs = []
        messages = []
        Q = self.getQ( )
        if( len( self ) < 2 ) : messages.append( 'endlI3.check: (E, mu, Ep, P) data len = %d < 2' % len( self ) )
        for E, muEpP in self.data :
            if( len( muEpP ) < 2 ) : messages.append( 'endlI3.check: (mu, Ep, P) data len = %d < 2 for E = %e' % ( len( muEpP ), E ) )
            for mu, EpP in muEpP :
                if( len( EpP ) < 2 ) : messages.append( 'endlI3.check: Ep, P len = %d < 2 for E = %e and mu = %e' % ( len( EpP ), E, mu ) )
                if( ( self.C >= 50 ) and ( self.C < 58 ) ) : continue
                if( EpP[-1][0] > E + Q ) : messages.append( "endlI3.check: E' = %.8e > E + Q = %e for E = %e, mu = %e" % \
                    ( EpP[-1][0], E + Q, E, mu ) )
            d2 = endlmath.ZSum( muEpP )
            for mu, y in d2.data :
                if( abs( 1. - y ) > normTolerance ) : messages.append( 'endlI3.check: bad normalize = %.8e for E = %e and mu = %e' % ( y, E, mu ) )
        messages += endl4dmathmisc.check4dData( self.data, allowNegativeT = False, allowZeroT = allowZeroE, allowNegativeX = True, allowZeroX = True, \
            allowNegativeY = False, positiveZ = False, printWarning = printWarning, printErrors = printErrors, xCloseEps = xCloseEps, \
            maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def copyData( self ) :
        """Returns an endlI3 instance - that is a copy, and not a reference - of self."""

        return endlI3( None, self.yo, self.C, self.I, self.S, self.h, endl4dmathClasses.endl4dmath.copyData( self ) )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor4d to fix threshold."""

        fixThresholdFor4d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin, fixThresholdMode = fixThresholdMode,
            threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getMuEpPAtE( self, E, extrapolation = endl3dmathmisc.noExtrapolation ) :
        """Returns an endl3dmath object of P'( mu, E' ) = P( E, mu, E' ).
    E is compare to the incident energies in self and is considered to be a match if the
    difference is less than relEps. Also see muEpP."""

        MuEpP, i, n = None, 0, len( self.data )
        if( n == 0 ) : return( MuEpP )
        while( i < n ) :
            if( E == self.data[i][0] ) :
                MuEpP = self.data[i][1]
                break
            elif( E < self.data[i][0] ) :
                if( i == 0 ) :
                    if( extrapolation == endl3dmathmisc.flatExtrapolation ) : MuEpP = self.data[i][1]
                else :
                    E1, MuEpP1 = self.data[i-1]
                    E2, MuEpP2 = self.data[i]
                    MuEpP1 = endl3dmathClasses.endl3dmath( MuEpP1, checkDataType = 0 )
                    MuEpP2 = endl3dmathClasses.endl3dmath( MuEpP2, checkDataType = 0 )
                    f = ( E2 - E ) / ( E2 - E1 )
                    g = 1 - f
                    mu1Min, mu1Max = MuEpP1.data[0][0], MuEpP1.data[-1][0]
                    mu2Min, mu2Max = MuEpP2.data[0][0], MuEpP2.data[-1][0]
                    muMin = max( -1., mu1Min + g * ( mu2Min - mu1Min ) )        # We are going to use unit-base interpolation of mu.
                    muMax = min(  1., mu1Max + g * ( mu2Max - mu1Max ) )        # However, must be careful to handle case where mu ranges from
                    dMu1 = ( mu1Max - mu1Min )                                  # -1 to 1 inclusively.
                    mus = [ muMax * ( ( mu1Max - mu ) / dMu1 ) + muMin * ( ( mu - mu1Min ) / dMu1 ) for mu, EpP in MuEpP1.data ]
                    dMu2 = ( mu2Max - mu2Min )
                    for mu, EpP in MuEpP2.data :
                        doIt = True
                        mup = muMax * ( ( mu2Max - mu ) / dMu2 ) + muMin * ( ( mu - mu2Min ) / dMu2 )
                        for m in mus :
                            if( abs( m - mup ) < 1e-14 ) : doIt = False
                        if( doIt ) : mus.append( mup )
                    mus.sort( )
                    if( mus[0] < -1 ) :
                        while( mus[0] < -1 ) : del mus[0]
                        mus.insert( 0, -1 )
                    if( mus[-1] > 1 ) :
                        while( mus[-1] > 1 ) : del mus[-1]
                        mus.append( 1 )
                    MuEpP = []
                    dMu = muMax - muMin
                    for mu in mus :
                        mu1 = mu1Min * ( ( muMax - mu ) / dMu ) + mu1Max * ( ( mu - muMin ) / dMu )
                        mu2 = mu2Min * ( ( muMax - mu ) / dMu ) + mu2Max * ( ( mu - muMin ) / dMu )
                        d = f * MuEpP1.getAtX( mu1, unitBase = True, endl2dmathObject = True ) \
                            + g * MuEpP2.getAtX( mu2, unitBase = True, endl2dmathObject = True )
                        MuEpP.append( [ mu, d.data ] )
                break
            i += 1
        if( ( i == n ) and ( extrapolation == endl3dmathmisc.flatExtrapolation ) ) : MuEpP = self.data[-1][1]
        if( not( MuEpP is None ) ) : MuEpP = endl3dmathClasses.endl3dmath( MuEpP, checkDataType = 0, xLabel = "mu",
            yLabel = "E'", zLabel = "P(mu,E') per E' (1/MeV)" )
        return( MuEpP )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = self.tMin( )
        EMinNext = None
        if( len( self.data ) > 1 ) : EMinNext = self.data[1][0]
        return( EMin, EMinNext )

    def getThresholdsForChecker( self ) :
        "For internal use only."

        if( len( self ) > 1 ) : return [ self.data[0][0], self.data[1][0] ]
        if( len( self ) > 0 ) : return [ self.data[0][0] ]
        return []

    def E( self, i ) :
        "Returns the (i+1)^th incident energy value as a float or None if i is out-of-range."

        if ( i < len( self.data ) ) : return float( self.data[i][0] )
        return None

    def EIndexMu( self, iE, mu, unitBase = True ) :
        """Returns an endl2dmath object that is the interpolation of self's data at mu for energy index iE."""

        if( ( iE < 0 ) or ( iE >= len( self.data ) ) ) :
            raise Exception( '\nError in endlI3.EIndexMu: energy index iE = %d out of range.' % len( self.data ) )
        MuEpP = self.data[iE][1]
        if(   mu <= MuEpP[0][0] ) :
            data = MuEpP[0][1]
        elif( mu >= MuEpP[-1][0] ) :
            data = MuEpP[-1][1]
        else :
            data = endl3dmathmisc.interpolate3d( mu, MuEpP, unitBase = unitBase )
        return( endl2dmathClasses.endl2dmath( data ) )

    def EMax( self ) :
        "Returns the last incident energy value for self. Returns None if data is empty."

        return( self.tMax( ) )

    def EMin( self ) :
        "Returns the first incident energy value for self. Returns None if data is empty."

        return( self.tMin( ) )

    def muEpP( self, E, relEps = 1e-5 ) :
        """Returns an endl3dmath object of P'( mu, E' ) = P( E, mu, E' ).
    E is compare to the incident energies in self and is considered to be a match if the
    difference is less than relEps. Also see getMuEpPAtE."""

        n = len( self.data )
        for i in range( n ) :
            if ( abs( self.data[i][0] - E ) <= relEps ) :       # ???? needs interpolation.
                return endl3dmathClasses.endl3dmath( self.data[i][1], checkDataType = 0, xLabel = "mu",
                    yLabel = "E'", zLabel = "P(mu,E') per E' (1/MeV)" )
        return None

    def mapMuEpPToGrid( self, E, relEps = 1e-5 ) :
        """P( E = E, mu, E' ) with mu, and E' mapped to a grid determine from
    the mu and E' points.  Returns a endl3dmath object of list[ u, list[ E', P( u, E' ) ] ].
    See muEpP( ) for meaning of relEps."""

        def trim( d ) :

            d.sort( )
            a = []
            e = None
            for i in d :
                if ( i != e ) :
                    a.append( i )
                    e = i
            return( a )

        d3i = self.muEpP( E, relEps )                       # Get 3-d data (u, E', P).
        if ( d3i == None ) : return None
        EpGrid = []                                         # All E' points.
        EpL = []                                            # Lower E' end point of data for all u.
        EpU = []                                            # Upper E' end point of data for all u.
        for u_etal in d3i.data :
            if ( u_etal[1][0][1]  != 0. ) : EpL.append( u_etal[1][0][0] )
            if ( u_etal[1][-1][1] != 0. ) : EpU.append( u_etal[1][-1][0] )
            for etal in u_etal[1] : EpGrid.append( etal[0] )
        EpGrid = trim( EpGrid )                             # Sort and remove redundant Ep points.
        EpLU = []                                           # Additional E points to add to make end-points go to 0.
        EpL = trim( EpL )                                   # Sort and remove redundant EpL points.
        EpL = EpL[1:]                                       # Do not need the lowest global E point.
        i = 1
        for Ep in EpL :
            while ( EpGrid[i] < Ep ) : i += 1
            if ( ( Ep - EpGrid[i-1] ) > 10 * endlParameters.endlEpsx * Ep ) : EpLU.append( Ep * ( 1 - endlParameters.endlEpsx ) )
        EpU = trim( EpU )
        EpU = EpU[:-1]
        EpU.reverse( )
        i = len( EpGrid ) - 2
        for Ep in EpU :
            while ( EpGrid[i] > Ep ) : i -= 1
            if ( ( EpGrid[i+1] - Ep ) > 10 * endlParameters.endlEpsx * Ep ) : EpLU.append( Ep * ( 1 + endlParameters.endlEpsx ) )
        if ( EpLU != [] ) : EpGrid = EpGrid + EpLU
        EpGrid = trim( EpGrid )
        d3o = []
        for u_etal in d3i.data :                            # Loop over all u.
            d2 = []
            E1 = None
            etal = u_etal[1]
            i = 0
            EP = etal[i]
            for E in EpGrid :                               # Loop over E' grid.
                P = 0.
                if ( E == EP[0] ) :                         # E' is a Grid energy.
                    E1 = E
                    P1 = EP[1]
                    P = P1
                    i += 1
                    if ( i < len( etal ) ) :
                        EP = etal[i]
                    else :
                        E1 = None
                elif ( E1 != None ) :                       # Interpolate.
                    P = ( EP[1] * ( E - E1 ) + P1 * ( EP[0] - E ) ) / ( EP[0] - E1 )
                d2.append( [ E, P ] )
            d3o.append( [ u_etal[0], d2 ] )
        return endl3dmathClasses.endl3dmath( d3o, checkDataType = 0, xLabel = "mu", yLabel = "E'", zLabel = "Probability" )

    def normalize( self ) :
        "Normalizes the data so that the integral P(E, mu, E') dE' = 1."

        for E_Data in self.data :                   # Loop for each E.
            for mu_Data in E_Data[1] :              # Loop for each mu.
                s = 0
                Ep0 = None
                for Data in mu_Data[1] :            # Loop for each E'.
                    Ep1 = Data[0]
                    C1 = Data[1]
                    if ( Ep0 != None ) : s += ( Ep1 - Ep0 ) * ( C1 + C0 )
                    Ep0 = Ep1
                    C0 = C1
                s /= 2.
                for Data in mu_Data[1] :
                    if( s != 0. ) : Data[1] /= s    # Loop for each E'.

    def convertToI4( self, i1=None, lmax=0 ):
        "Converts the self and the I=1 data given as an argument to an endlI4 instance where the new object has Legendre orders 0, 1, ..., lmax"

        from xData import series1d as seriesModule

        # Set up the I=4 file, especially the header crap
        dummy = endlI3( None, self.yo, self.C, self.I, self.S, self.h, [] )
        dummy.setI( 4 )
        i4 = endlI4( None, self.yo, self.C, 4, self.S, dummy.h, [] )
        endlmisc.copyEndlHeader( i4, self )

        # Loop through l's expanding the product of the I=1 and I=3 files into Legendre polynomials
        for l in range( 0, lmax+1 ):
            thisLTerm = endl3dmathClasses.endl3dmath()
            legPoly = endl2dmathmisc.convertFunctionToLinLin( lambda x: ( 2 * l + 1 ) * seriesModule.Legendre( l, x ), -1.0, 1.0, 1e-4 )
            for iE in range( len ( self.data ) ):
                E = self.E( iE )
                i1_muP = endl2dmathClasses.endl2dmath( data = i1.data[ iE ][1] ) 
                i3_muEpP = self.muEpP( E, relEps = 1e-5 )
                # collect list of Ep values
                EpList = []
                for muEpPpair in i3_muEpP.data: EpList += filter( lambda x: x not in EpList, [ EpP[0] for EpP in muEpPpair[1] ] )
                EpList.sort()
                # Integrate in mu for each Ep value
                EpP = []
                for Ep in EpList :
                    i3_muP = endl2dmathClasses.endl2dmath( data = [ [ x, i3_muEpP.getValue( x, Ep ) ] for x in i3_muEpP.xArray() ] )
                    P = i1_muP.integrateThreeFunctions( i3_muP, legPoly, -1., 1. )
                    EpP.append( [ Ep, P ] )
                thisLTerm.data.append( [ E, EpP ] )
            i4.setlData( l, thisLTerm )    
        return i4

    def reduceToEEpP( self ) :
        "Integrates the mu dimension, returning an endl3dmath object of list[ E, list[ E', P( E, E' ) ] ]."

        d3 = []
        for e_etal in self.data :
            d = self.mapMuEpPToGrid( e_etal[0], relEps = 1e-14 )
            u = None
            for u_etal in d.data :
                if ( u == None ) :
                    d3p = []
                    for etal in u_etal[1] : d3p.append( [ etal[0], etal[1], 0. ] )
                else :
                    i = 0
                    for etal in u_etal[1] :
                        dd = d3p[i]
                        dd[2] = dd[2] + ( etal[1] + dd[1] ) * ( u_etal[0] - u )
                        dd[1] = etal[1]
                        i = i + 1
                u = u_etal[0]
            d2 = []
            for d in d3p : d2.append( [ d[0], d[2] / 2. ] )
            d3.append( [ e_etal[0], d2 ] )
        return endl3dmathClasses.endl3dmath( d3, checkDataType = 0, xLabel = "E", yLabel = "E'", zLabel = "Probability" )

    def reduceToEMuP( self, normalize = False ) :                          # class endlI3
        "Integrates the E' dimension, returning an unnormalized (normalize = False) or normalized (normalize = True) endlI1 object."

        d3 = []
        for e_etal in self.data :
            d2 = []
            for u_etal in e_etal[1] :
                Ep1 = None
                for etal in u_etal[1] :
                    if ( Ep1 == None ) :
                        s = 0
                    else :
                        s += ( etal[1] + P1 ) * ( etal[0] - Ep1 )
                    Ep1 = etal[0]
                    P1 = etal[1]
                d2.append( [ u_etal[0], s / 2. ] )
            d3.append( [ e_etal[0], d2 ] )
        dummy = endlI3( None, self.yo, self.C, self.I, self.S, self.h, [] )
        dummy.setI( 1 )
        I1 = endlI1( None, self.yo, self.C, 1, self.S, dummy.h, d3 )
        endlmisc.copyEndlHeader( I1, self )
        if( normalize ) : I1.normalize( )
        return I1

    def runningEpSum( self ) :
        """Does a running sum of the E' data, returning an endl4dmath object of the results."""

        d4 = []
        for E_uEpP in self.data :
            d3 = []
            for u_EpP in E_uEpP[1] :
                d3.append( [ u_EpP[0], endlmath.runningYSum( u_EpP[1] ).data ] )
            d4.append( [ E_uEpP[0], d3 ] )
        return endl4dmathClasses.endl4dmath( d4, tLabel = self.tLabel, xLabel = self.xLabel, yLabel = self.yLabel,
            zLabel = "running Sum of P( E, mu, E' ) vs E'" )

    def plot( self, E = None, xMin = None, xMax = None, yMin = None, yMax = None, zMin = None, zMax = None, \
        xyzlog = 0, tLabel = None, xLabel = None, yLabel = None, zLabel = None, title = None, \
        xrot = None, zrot = None ) :                                                                # class endlI3
        """Plots the data.  If E = None uses 4d plotting otherwise uses 3d plotting for the requested E.

xyzlog values and meaning::
    xyzlog   plot-type for x-y-z axis
   -----------------------------------
      0     linear-linear-linear
      1     log-linear-linear
      2     linear-log-linear
      3     log-log-linear
      4     linear-linear-log
      5     log-linear-log
      6     linear-log-log
      7     log-log-log"""

        if ( tLabel == None ) and ( self.tLabel != None ) : tLabel = self.tLabel
        if ( xLabel == None ) and ( self.xLabel != None ) : xLabel = self.xLabel
        if ( yLabel == None ) and ( self.yLabel != None ) : yLabel = self.yLabel
        if ( zLabel == None ) and ( self.zLabel != None ) : zLabel = self.zLabel

        if ( E != None ) :
            iE = None
            i = -1
            for E_xyz in self.data :
                i += 1
                if ( abs( E - E_xyz[0] ) < E * 1e-8 ) : iE = i
            if ( iE == None ) : raise Exception( "\nError in endlI3.plot: E value not found" )
            dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, \
                zMin = zMin, zMax = zMax, zLabel = zLabel, xrot = xrot, zrot = zrot )
            endl4dmathmisc.plot3dFrom4d( self.data, iE, dt, xyzlog = xyzlog )
        else :
            endl4dmathClasses.endl4dmath.plot( self, xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax, zMin = zMin, zMax = zMax, xyzlog = xyzlog, \
                tLabel = tLabel, xLabel = xLabel, yLabel = yLabel, title = title, tScaleLabel = "Incident Energy", xrot = xrot, zrot = zrot )

    def set( self, data, checkDataType = 0, tLabel = "E (MeV)", xLabel = "mu", yLabel = "E' (MeV)", zLabel = "pdf(E,mu,E') vs E' (1/MeV)", interpolation = 0 ) :

        endl4dmathClasses.endl4dmath.__init__( self, data, checkDataType = checkDataType, tLabel = tLabel, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

class endlI4( endlNd.endlNd, endl4dmathClasses.endl4dmath ) :
    """This class is for the I = 4 data."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI4 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl4dmath class data structure."""

        from fudge import gnd

        self.name = 'pointwise'
        endlNd.endlNd.__init__( self, f, 4, yo, C, I, S, h, points, i0 = 2, i1 = 0, i2 = 1, i3 = 3, bdflsFile = bdflsFile )
        self.xInterpolation = 'none,linear,unitbase:x:linear,unitbase:y:linear'
        self.variablesUnits = 'l;energy_in[MeV];energy_out[MeV];C_l(energy_out|energy_in)'

    def check( self, normTolerance = 1e-5, printWarning = False, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 4 data and returns a list of
        endlCheckerObject instances. Also, calls the endl4dmathmisc.check4dData function.
        See endl4dmathmisc.check4dData for meaning of printWarning and printErrors."""

        normTolerance = max( normTolerance, normCheckTolerance )
        ErrMsgs = []
        messages = []
        Q = self.getQ( )
        for l, EEpP in self.data :
            if( len( EEpP ) < 2 ) : messages.append( 'endlI4.check: (E, Ep, P) data len = %d < 2 for l = %d' % ( len( EEpP ), l ) )
            for E, EpP in EEpP :
                if( len( EpP ) < 2 ) : messages.append( 'endlI4.check: Ep, P len = %d < 2 for l = %d and E = %e' % ( len( EpP ), l, E ) )
                if( ( self.C >= 50 ) and ( self.C < 58 ) ) : continue
                if( EpP[-1][0] > E + Q ) : messages.append( "endlI4.check: E' = %.8e > E + Q = %e for l = %.0f, E = %e" % \
                    ( EpP[-1][0], E + Q, l, E ) )
        d2 = endlmath.ZSum( self.data[0][1] )
        for E, y in d2.data :
            if( abs( 1. - y ) > normTolerance ) : messages.append( 'endlI4.check: bad normalize = %.8e for l = 0 and E = %e' % ( y, E ) )
        messages += endl4dmathmisc.check4dData( self.data, allowNegativeT = False, allowZeroT = False, allowNegativeX = False, allowZeroX = allowZeroE, \
            allowNegativeY = False, positiveZ = False, printWarning = printWarning, printErrors = printErrors, xCloseEps = xCloseEps, \
            maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def convertToI1I3( self, nMu = 21, lMax = None ) :
        """This methods converts the I = 4 data into the endl I = 1 and I = 3 data using nMu equally spaced mu values.
        This method calls convertToUI3. The returned value is the tuple (I1, I3)."""

        I3 = self.convertToUI3( nMu = nMu, lMax = lMax )
        I1Data = []
        for E, muEpP in I3.data :
            I1MuP = []
            for mu, EpP in muEpP :
                I1MuP.append( [ mu, endlmath.YSum( EpP ) ] )
            I1Data.append( [ E, I1MuP ] )
        dummy = endlI4( None, self.yo, self.C, self.I, self.S, self.h, [] )
        dummy.setI( 1 )
        I1 = endlI1( None, self.yo, self.C, 1, self.S, dummy.h, I1Data )
        I1.normalize( )
        I3.normalize( )
        return( I1, I3 )

    def convertToUI3( self, nMu = 21, lMax = None ) :
        """This methods converts the I = 4 data into unnormalized I = 3 data using nMu equally spaced mu 
        values.  Note that the data is not normalized; hence, the normalized method should be called on the
        returned object. This method will only convert up to Legendre order l = min( lMax, series1d.maxLegendreOrder ) 
        of self's data."""

        from xData import series1d as seriesModule

        if( lMax == None ) : lMax = seriesModule.maxLegendreOrder
        lMax = max( 0, min( lMax, seriesModule.maxLegendreOrder ) )
            
        I3Data = []
        if ( len( self ) == 1 ) :
            mu_EpP = []
            for E_EpP in self.data[0][1] :
                EpPm1 = []
                EpPp1 = []
                for EpP in E_EpP[1] :
                    EpPm1.append( [ EpP[0], EpP[1] ] )
                    EpPp1.append( [ EpP[0], EpP[1] ] )
                I3Data.append( [ E_EpP[0], [ [ -1., EpPm1 ], [ 1., EpPp1 ] ] ] )
        else :
            EArray = []                                 # Array of incident energies
            for E, EpP in self.data[0][1] : EArray.append( E )
            for l, EEpP in self.data :                   # Check that all l's have the same incident energies
                if ( len( EEpP ) != len( EArray ) ) :
                    raise Exception( "\nError in endlI4.convertToUI3: incident energies do not align" )
                i = 0
                for E, EpP in EEpP :
                    if ( E != EArray[i] ) : raise Exception( "\nError in endlI4.convertToUI3: incident energies do not align (l = %d)" % l )
                    i += 1
            muArray = []                                # Array of mus
            rMu = range( nMu )
            s = 2. / ( nMu - 1 )
            for iu in rMu : muArray.append( -1. + s * iu )
            muArray[nMu - 1] = 1.
            l_LArray = []
            for l, EEpP in self.data :
                if( l > lMax ) : break
                LArray = []
                for mu in muArray : LArray.append( seriesModule.Legendre( l, mu ) )
                l_LArray.append( LArray )
            iE = -1
            for E in EArray :
                iE += 1
                Ep_Array = []                           # Construct E' grid for E.
                for l, EEpP in self.data :
                    if( l > lMax ) : break
                    for Ep, P in EEpP[iE][1] :
                        if ( Ep not in Ep_Array ) : Ep_Array.append( Ep )
                Ep_Array.sort( )
                mu_EpP = []                             # Construct mu, [ E', P ] grid
                for mu in muArray :
                    EpP = []
                    for Ep in Ep_Array : EpP.append( [ Ep, 0. ] )
                    mu_EpP.append( [ mu, EpP ] )
                iEp = -1
                for Ep in Ep_Array :
                    iEp += 1
                    il = -1
                    for l, EEpP in self.data :
                        if( l > lMax ) : break
                        il += 1
                        iEpP = 0
                        EpP_ = EEpP[iE][1]
                        for EpP in EpP_ :
                            if ( Ep <= EpP[0] ) : break
                            iEpP += 1
                        if ( Ep == EpP_[iEpP][0] ) :
                            P = EpP_[iEpP][1]
                        else :
                            P = ( EpP_[iEpP-1][1] * ( EpP_[iEpP][0]  - Ep ) + EpP_[iEpP][1] * ( Ep - EpP_[iEpP-1][0] ) ) / \
                                ( EpP_[iEpP][0] - EpP_[iEpP-1][0] )
                        for iu in rMu : mu_EpP[iu][1][iEp][1] += P * l_LArray[il][iu]
                I3Data.append( [ E, mu_EpP ] )
        dummy = endlI4( None, self.yo, self.C, self.I, self.S, self.h, [] )
        dummy.setI( 3 )
        I3 = endlI3( None, self.yo, self.C, 3, self.S, dummy.h, I3Data )
        endlmisc.copyEndlHeader( I3, self )
        return I3

    def copyData( self ) :
        """Returns an endlI4 instance that is a copy, and not a reference, of self."""

        return endlI4( None, self.yo, self.C, self.I, self.S, self.h, endl4dmathClasses.endl4dmath.copyData( self ) )

    def EMax( self ) :
        "Returns the last incident energy value for self. Returns None if data is empty."

        EMax = None
        for l, EEpP in self.data :
            if( EMax == None ) :
                if( len( EEpP ) > 0 ) : EMax = EEpP[0][0]
            else :
                if( len( EEpP ) > 0 ) : EMax = max( EMax, EEpP[0][0] )
        return( EMax )

    def EMin( self ) :
        "Returns the first incident energy value for self. Returns None if data is empty."

        EMin = None
        for l, EEpP in self.data :
            if( EMin == None ) :
                if( len( EEpP ) > 0 ) : EMin = EEpP[0][0]
            else :
                if( len( EEpP ) > 0 ) : EMin = min( EMin, EEpP[0][0] )
        return( EMin )

    def fixThreshold( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV = 10e-3, EMin = 0., 
            fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :
        """Calls endlIClasses.fixThresholdFor3d to fix threshold for each l-order data."""

        for l, EEpP in self.data :
            d = endl3dmathClasses.endl3dmath( data = EEpP, checkDataType = 0 )
            fixThresholdFor3d( d, thresholdCrossSectionIsZero, self.I, threshold, dThreshold_MeV = dThreshold_MeV, EMin = EMin, realSelf = self, fixThresholdMode = fixThresholdMode,
                threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getEMin_EMinNext( self ) :
        """Gets the first two Energy values from self's data. None is return for any absent value."""

        EMin = None
        EMinNext = None
        for l, EEpP in self.data :
            n = len( EEpP )
            if( n > 0 ) :
                EMin_ = EEpP[0][0]
                if( ( EMin == None ) or ( EMin_ < EMin ) ) : EMin = EMin_
                if( n > 1 ) :
                    EMinNext_ = EEpP[1][0]
                    if( ( EMinNext == None ) or ( EMinNext_ < EMinNext ) ) : EMinNext = EMinNext_
        return( EMin, EMinNext )

    def getlData( self, l ) :
        """Returns an endl3dmath object for the data at Legendre order l. If l is outside the domain of the data, then None is returned."""

        for l_, EEpP in self.data :
            if( l_ == l ) : 
                label = 'l = %d' % l
                if( hasattr, self, 'label' ) : label = self.label + ' for l = %d' % l
                return( endl3dmathClasses.endl3dmath( EEpP, xLabel = self.xLabel, yLabel = self.yLabel, zLabel = self.zLabel, label = label ) )
        return( None )

    def getl_EData( self, l, E, unitBase = True ) :
        """Returns an endl2dmath object for the data at Legendre order l and incident energy E. If l or E is outside the
        domain of the data, then None is returned.  If the requested E is not a point in the data then linear interpolation 
        is performed. If unitBase is true, the Ep, P(l,E,Ep) data are unit based linear interpolated."""

        EEpP = self.getlData( l ).data
        if( EEpP == None ) : return( None )
        nE = len( EEpP )
        if( nE == 0 ) : return( None )
        if( ( E < EEpP[0][0] ) or ( E > EEpP[-1][0] ) ) : return( None )
        return( endl2dmathClasses.endl2dmath( endl3dmathmisc.interpolate3d( E, EEpP, unitBase = unitBase ), xLabel = self.yLabel, yLabel = self.zLabel ) )

    def setlData( self, l, EEpP ) :
        """Set data for Legendre order l to EEpP. EEpP can be an endl3dmath object or a suitable python list. The l-value
        must be an existing l-value or one greater than the current maximum l-value."""

        data = endl3dmathmisc.get3dmathData( EEpP, 'endlI4.setlData', 'EEpP' )
        nl = len( self.data )
        if( ( l < 0 ) or ( l > nl ) ) : raise Exception( '\nError in endlI4.setlData: l = %d must be between 0 to lMax + 1 (= %d) inclusive.' % ( l, nl ) )
        if( l == nl ) :                         # Adding the next l order.
            self.data.append( [ l, data ] )
        else :
            for i in xrange( nl ) :
                if( l == self.data[i][0] ) :
                    self.data[i][1] = data
                    return

    def setl_EData( self, l, E, EpP ) :
        """Set data for Legendre order l and incident energy E to EpP. EpP can be an endl2dmath object or a suitable python list.
        The l-value must be an existing l-value or one greater than the current maximum l-value."""

        data = endl2dmathmisc.get2dmathData( EpP, 'endlI4.setl_EData', 'EpP' )
        nl = len( self.data )
        if( ( l < 0 ) or ( l > nl ) ) : raise Exception( '\nError in endlI4.setl_EData: l = %d must be between 0 to lMax + 1 (= %d) inclusive.' % ( l, nl ) )
        if( l == nl ) :
            self.data.append( [ l, [ [ E, data ] ] ] )
        else :
            for l_, EEpP in self.data :
                if( l == l_ ) :
                    nE = len( EEpP )
                    if( nE == 0 ) :
                        EEpP.append( [ E, data ] )
                    elif( E < EEpP[0][0] ) :
                        EEpP.insert( 0, [ E, data ] )
                    elif( E > EEpP[-1][0] ) :
                        EEpP.append( [ E, data ] )
                    else :
                        for i in xrange( nE ) :
                            if( E <= EEpP[i][0] ) : break
                        if( E == EEpP[i][0] ) :
                            EEpP[i][1] = data
                        else :
                            EEpP.insert( i, [ E, data ] )
                    return

    def getThresholdsForChecker( self ) :
        "For internal use only."

        thresholds = []
        for l, EAndOthers in self.data :
            if( len( EAndOthers ) > 0 ) : thresholds.append( EAndOthers[0][0] )
        return thresholds

    def l( self, i ) :
        """Returns the (i+1)^th Legendre l-value (e.g., self.l( 0 ) returns the first l-value) or None if i is out-of-range."""

        if ( i < len( self.data ) ) : return self.data[i][0]
        return None

    def normalize( self, doL_0_Only = True ) :
        """Normalizes the l = 0 data so that the integral P(l = 0, E, E') dE' = 1. If doL_0_Only is true, then only the
        l = 0 component is normalizes. Otherwise, the l = 0 sum is used to normalize all orders. The return value the tuple
        (n0, nm) where n0 (nm) is the number of sum that are zero (negative)."""

        n0 = 0
        nm = 0
        if ( self.data[0][0] == 0 ) :
            iE = -1
            for d0 in self.data[0][1] :         # Loop for each E.
                iE += 1
                s = 0
                Ep0 = None
                for d1 in d0[1] :               # Loop for each E'.
                    Ep1 = d1[0]
                    C1 = d1[1]
                    if ( Ep0 != None ) : s += ( Ep1 - Ep0 ) * ( C1 + C0 )
                    Ep0 = Ep1
                    C0 = C1
                s /= 2.
                if( s == 0. ) :
                    n0 += 1
                    s = 1.
                elif( s < 0. ) :
                    nm += 1
                for d1 in d0[1] : d1[1] /= s    # Loop for each E'.
                if( not doL_0_Only ) :
                    for l, EEpP in self.data :
                        if( l == 0 ) : continue
                        EpP = EEpP[iE][1]               # Note, index of iE and E not checked. 
                        for d1 in EpP : d1[1] /= s      # Loop for each E'.
        return( n0, nm )

    def plot( self, l = None, xMin = None, xMax = None, yMin = None, yMax = None, zMin = None, zMax = None, \
        xyzlog = 0, tLabel = None, xLabel = None, yLabel = None, zLabel = None, title = None, \
        xrot = None, zrot = None ) :
        """Plots the data for Legendre order l.

xyzlog values and meaning::
    xyzlog   plot-type for x-y-z axis
   -----------------------------------
      0     linear-linear-linear
      1     log-linear-linear
      2     linear-log-linear
      3     log-log-linear
      4     linear-linear-log
      5     log-linear-log
      6     linear-log-log
      7     log-log-log"""

        if ( tLabel == None ) and ( self.tLabel != None ) : tLabel = self.tLabel
        if ( xLabel == None ) and ( self.xLabel != None ) : xLabel = self.xLabel
        if ( yLabel == None ) and ( self.yLabel != None ) : yLabel = self.yLabel
        if ( zLabel == None ) and ( self.zLabel != None ) : zLabel = self.zLabel

        if ( l != None ) :
            dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, \
                zMin = zMin, zMax = zMax, zLabel = zLabel, xrot = xrot, zrot = zrot )
            endl4dmathmisc.plot3dFrom4d( self.data, l, dt, xyzlog = xyzlog )
        else :
            endl4dmathClasses.endl4dmath.plot( self, xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax, zMin = zMin, zMax = zMax, \
                xyzlog = xyzlog, tLabel = tLabel, xLabel = xLabel, yLabel = yLabel, title = title, tScaleLabel = "'l order'", \
                xrot = xrot, zrot = zrot )

    def reduceToEEpP( self ) :
        """reduceToEEpP( )\n    Returns the l = 0 term as an endl3dmath object of list[ E, list[ E', P( E, E' ) ] ].
    Note, the l = 0 term is P( E, E' )."""

        return endl3dmathClasses.endl3dmath( self.data[0][1], checkDataType = 0, xLabel = "E",
            yLabel = "E'", zLabel = "Probability( E, E' )" )

    def set( self, data, checkDataType = 0, tLabel = "l", xLabel = "E (MeV)", yLabel = "E' (MeV)", zLabel = "Legendre coef. C_l(E,E') vs E' (1/MeV)",
        interpolation = 0 ) :

        endl4dmathClasses.endl4dmath.__init__( self, data, checkDataType = checkDataType, tLabel = tLabel, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )

    def toString( self, format = None ) :
        """Returns a string with the data in the order (E, E', l, P) as required by ENDL."""

        s = '\n'.join( endlmisc.string4dData( self.data, i0 = 2, i1 = 0, i2 = 1, i3 = 3, fmt0 = format, fmt1 = format, fmt2 = format, fmt3 = format ) )
        return( s + '\n' )

    def toZAsFrame( self, newProjectileMass, newTargetMass, halflife, bdflsFile, ELevel = 0. ) :

        newI4, r = toZAsFrameMisc( endlI4, self, newProjectileMass, newTargetMass, ELevel, halflife, bdflsFile, False )
        for l, EEpP_ in newI4.data :
            for EEpP in EEpP_ : EEpP[0] *= r
        return( [ newI4 ] )

class endlI20( endlNd.endlNd, endl4dmathClasses.endl4dmath ) :
    """This class is for the I = 20 data. I 20 data is unresolved resonance probability table data which
    has four columns. These columns are 1) s incident energy, 2) Temperature, 3) probability and 
    4) cross section."""

    def __init__( self, f, yo, C, I, S, h, points, bdflsFile = None ) :
        """Constructor for the endlI20 class.  See the module endlNd.py for the meanings of
        f, yo, C, I, S and h. Points must be a valid endl4dmath class data structure."""

        self.name = 'URRProbabilityTable'
        endlNd.endlNd.__init__( self, f, 20, yo, C, I, S, h, points, i0 = 0, i1 = 1, i2 = 2, i3 = 3, bdflsFile = bdflsFile )
        self.variablesUnits = 'energy_in[MeV];temperature[keV];probability;cross_section[barn]'

    def check( self, normTolerance = 1e-4, printWarning = False, printErrors = True, xCloseEps = None, allowZeroE = False, maxAbsFloatValue = None, **arg ) :
        """Checks to see that the data is consistance with I = 20 data and returns a list of
        endlCheckerObject instances. Also, calls the endl4dmathmisc.check4dData function.
        See endl4dmathmisc.check4dData for meaning of printWarning and printErrors."""

        normTolerance = max( normTolerance, normCheckTolerance )
        ErrMsgs = []
        messages = []
        for E, TPXsec in self.data :
            if( len( TPXsec ) < 1 ) : messages.append( 'endlI20.check: (T, P, Xsec) data len = %d < 1 for E = %e' % ( len( EPXsec ), E ) )
            for T, PXsec in TPXsec :
                if( len( PXsec ) < 1 ) : 
                    messages.append( 'endlI20.check: (P, Xsec) len = %d < 1 for E = %e and T = %e' % ( len( PXsec ), E, T ) )
                else :
                    if( PXsec[0][0] < 0 ) : messages.append( 'endlI20.check: bad normalize = %.8e for E = %e and T = %e' % ( PXsec[0][0], E, T ) )
                    PPrior = 0.
                    sum = 0
                    for P, Xsec in PXsec :
                        if( Xsec < 0 ) : messages.append( 'endlI20.check: negative cross section = %e for E = %e and T = %e' % ( Xsec, E, T ) )
                        if( PPrior > P ) : messages.append( 'endlI20.check: bad cdf (p[i-1] = %e > p[i] = %e for E = %e and T = %e' % ( PPrior, P, E, T ) )
                        sum += ( P - PPrior ) * Xsec
                        PPrior = P
                    if( abs(P-1) > normTolerance ) : messages.append( 'endlI20.check: bad probability normalization = %e for E = %e and T = %e' % ( P, E, T ) )
                    if( ( sum != 0 ) and ( abs(sum-1) > normTolerance ) ) : 
                        messages.append( 'endlI20.check: bad normalization = %e for E = %e and T = %e' % ( sum, E, T ) )
        messages += endl4dmathmisc.check4dData( self.data, allowNegativeT = False, allowZeroT = allowZeroE, allowNegativeX = False, allowZeroX = True, \
            allowSameY = True, allowNegativeY = False, positiveZ = True, printWarning = printWarning, printErrors = printErrors, xCloseEps = xCloseEps, \
            maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages ) > 0 ) : ErrMsgs = [ endlmisc.endlCheckerObject( data = self, message = messages ) ]
        return( ErrMsgs )

    def copyData( self ) :
        """Returns an endlI20 instance that is a copy, and not a reference, of self."""

        return endlI20( None, self.yo, self.C, self.I, self.S, self.h, endl4dmathClasses.endl4dmath.copyData( self ) )

    def normalize( self ) :
        """Normalizes the data so that the cdf's are 1."""

        for E, TPXsec in self.data :
            for T, PXsecs in TPXsec :
                n = 1. / PXsec[-1][0]
                for PXsec in PXsecs : PXsec[0] /= n
                n = 0
                for P, Xsec in PXsecs : n += P * Xsec 
                if( n != 0 ) :
                    for PXsec in PXsecs : PXsec[1] /= n

    def set( self, data, checkDataType = 0, tLabel = "E (MeV)", xLabel = "Temperature (keV)", yLabel = "probability", zLabel = "cross section (barn)",
        interpolation = 0 ) :

        endl4dmathClasses.endl4dmath.__init__( self, data, checkDataType = checkDataType, tLabel = tLabel, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel,
            interpolation = interpolation )
        self.numberOfEnergies = 0
        self.numberOfTemperatures = 0
        self.numberOfProbabilities = 0
        if( len( data ) > 0 ) :
            self.numberOfEnergies = len( self )
            self.numberOfTemperatures = len( self.data[0][1] )
            self.numberOfProbabilities = len( self.data[0][1][0][1] )

    def toString( self, format = None ) :
        """Returns a string with the data in the order (E, T, P, Xsec) as required by ENDL."""

        s = '\n'.join( endlmisc.string4dData( self.data, i0 = 0, i1 = 1, i2 = 2, i3 = 3, fmt0 = format, fmt1 = format, fmt2 = format, fmt3 = format ) )
        return( s + '\n' )

def getThresholdsForChecker2d( self ) :

    if( len( self ) > 1 ) : return [ self.data[0][0], self.data[1][0] ]
    if( len( self ) > 0 ) : return [ self.data[0][0] ]
    return []

def fixThresholdTest( threshold, dThreshold_MeV, EMin, EMin_, fixThresholdMode ) :

    if( fixThresholdMode is None ) : fixThresholdMode = fixThresholdMode_RaiseOnly
    if( fixThresholdMode == fixThresholdMode_None ) :
        return( False )
    elif( fixThresholdMode == fixThresholdMode_RaiseOnly ) :
        fixThreshold = threshold - EMin_ > dThreshold_MeV
    else :
        fixThreshold = abs( threshold - EMin_ ) > dThreshold_MeV
    if( fixThreshold or ( EMin_ < EMin ) ) : return( True )
    return( False )

def fixThresholdFor2d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, thresholdValue, EMin = 0., fixThresholdMode = fixThresholdMode_RaiseOnly, 
        threshold_MeV_shiftWarning = 0.1 ) :

    if( len( self ) == 0 ) : return
    if( EMin > threshold ) : threshold = EMin
    EMin_ = self.xMin( )
    if( fixThresholdTest( threshold, dThreshold_MeV, EMin, EMin_, fixThresholdMode ) ) :
        dE = threshold - EMin_
        if( abs( dE ) > threshold_MeV_shiftWarning ) :
            endlmisc.printWarning( '    2d: Moving threshold %e by %e for ZA = %d, %s' % ( self.data[0][0], dE, self.ZA, `self` ) )
        if( dE < 0 ) :
            if( not( thresholdCrossSectionIsZero ) ) :
                endlmisc.printWarning( '    Warning 2d: adding point at threshold where prior cross section was not 0: %s' % `self` )
            if( thresholdValue is None ) : thresholdValue = self.data[0][1]
            self.data.insert( 0, [ threshold, thresholdValue ] )
        else :
            ELast = EMin_ + 1.                                 # Make sure the first point is always done.
            for Ey in self.data :
                if( Ey[0] > ELast ) : break
                Ey[0] += dE
                ELast = Ey[0]

def fixThresholdFor3d( self, thresholdCrossSectionIsZero, I, threshold, dThreshold_MeV, EMin = 0., realSelf = None, fixThresholdMode = fixThresholdMode_RaiseOnly,
        threshold_MeV_shiftWarning = 0.1 ) :

    if( len( self ) == 0 ) : return
    if( EMin > threshold ) : threshold = EMin
    EMin_ = self.xMin( )
    if( fixThresholdTest( threshold, dThreshold_MeV, EMin, EMin_, fixThresholdMode ) ) :
        self_ = self
        if( realSelf is not None ) : self_ = realSelf
        dE = threshold - EMin_
        if( abs( dE ) > threshold_MeV_shiftWarning ) :
            endlmisc.printWarning( '    3d: Moving threshold %e by %e for ZA = %d, %s' % ( self.data[0][0], dE, self_.ZA, `self_` ) )
        if( dE > 0. ) :                                     # Move lower data up.
            ELast = EMin_ + 1.                              # Make sure the first point is always done.
            EyzPrepend = None
            for Eyz in self.data :
                if( Eyz[0] > ELast ) : break
                Eyz[0] += dE
                ELast = Eyz[0]
            if( EyzPrepend is not None ) : self.data.insert( 0, EyzPrepend )
        else :                                              # Create point at threshold.
            if( not( thresholdCrossSectionIsZero ) ) :
                endlmisc.printWarning( '    Warning 3d: adding point at threshold where prior cross section was not 0: %s' % `self_` )
            delta = fixThreshold_deltaFunctionEpsilon
            if( I == 1 ) :
                if( ( self_.C == 10 ) or ( self_.S == 1 ) ) :
                    endlmisc.printWarning( '    3d: Adding isotropic at threshold for com data: %s' % `self_` )
                    self.data.insert( 0, [ threshold, [ [ -1, 0.5 ], [ 1.0, 0.5 ] ] ] )
                else :
                    endlmisc.printWarning( '    3d: Adding forward peaked delta function at threshold: %s' % `self_` )
                    self.data.insert( 0, [ threshold, [ [ 1 - delta, 0.0 ], [ 1.0, 2 / delta ] ] ] )
            else :
                if( I in [ 21, 22 ] ) :
                    endlmisc.printWarning( '    3d: I = %d data not support by fixThresholdFor3d: ZA = %d: %s' % ( I, self_.ZA, `self_` ) )
                else :
                    endlmisc.printWarning( '    3d: Creating data at threshold %e for I = %d data' % ( threshold, I ) )
                    self.data.insert( 0, [ threshold, [ [ delta, 2 / delta ], [ 2 * delta, 0.0 ] ] ] )

def fixThresholdFor4d( self, thresholdCrossSectionIsZero, threshold, dThreshold_MeV, EMin = 0., fixThresholdMode = fixThresholdMode_RaiseOnly, threshold_MeV_shiftWarning = 0.1 ) :

    if( len( self ) == 0 ) : return
    if( EMin > threshold ) : threshold = EMin
    EMin_ = self.EMin( )
    if( fixThresholdTest( threshold, dThreshold_MeV, EMin, EMin_, fixThresholdMode ) ) :
        dE = threshold - EMin_
        if( abs( dE ) > threshold_MeV_shiftWarning ) :
            endlmisc.printWarning( '    4d: Moving threshold %e by %e for ZA = %d, %s' % ( self.data[0][0], dE, self.ZA, `self` ) )
        if( dE > 0 ) :                                      # Move lower data up.
            ELast = EMin_ + 1.                              # Make sure the first point is always done.
            for Exyz in self.data :
                if( Exyz[0] > ELast ) : break
                Exyz[0] += dE
                ELast = Exyz[0]
        else :                                              # Create point at threshold.
            if( not( thresholdCrossSectionIsZero ) ) :
                endlmisc.printWarning( '    Warning 4d: adding point at threshold where prior cross section was not 0: %s' % `self` )
            delta = fixThreshold_deltaFunctionEpsilon
            if( self.I == 3 ) :
                yo = self.yo
                if( yo > 9 ) : yo -= 10
                mass1, mass2, mass3 = self.bdflsFile.mass( self.yi ), self.bdflsFile.mass( self.ZA ), self.bdflsFile.mass( yo )
                EpLab = mass1 * mass3 * threshold / ( mass1 + mass2 )**2
                if( EpLab - delta < 0 ) :
                    EpP = [ [ 0.0, 2 / delta ], [ delta, 0 ] ]
                else :
                    EpP = [ [ EpLab - delta, 0.0 ],
                            [         EpLab, 1 / delta ],
                            [ EpLab + delta, 0.0 ] ]
                self.data.insert( 0, [ threshold, [ [ 1 - delta, EpP ],
                                                    [       1.0, EpP ] ] ] )
            else :
                endlmisc.printWarning( '    4d: Creating data at threshold not supported for %s' % `self` )

def I1sToCommonGrids( I1List, muGridPerE = True ) :
    """This routine takes a list of I = 1 data and puts their E and mu data onto a common grid. Note, this function may add
    points for which integral( d_mu P(E,mu) ) = 0 and not 1 (that is, the normalization method may fail)."""

    Es = []
    muGrids = {}
    muMasterGrid = []
    for I1 in I1List :
        for E, muP in I1.data :
            if( E not in Es ) :
                Es.append( E )
                muGrids[E] = []
            muGrid = muGrids[E]
            for mu, P in muP :
                if( mu not in muGrid ) : muGrid.append( mu )
                if( mu not in muMasterGrid ) : muMasterGrid.append( mu )
    Es.sort( )
    for E in Es : muGrids[E].sort( )
    muMasterGrid.sort( )
    thisMuGrid = muMasterGrid
    for I1 in I1List :
        EMin, EMax = I1.xMin( ), I1.xMax( )
        if( ( EMin == None ) or ( EMax == None ) ) : raise Exception( 'I = 1 data does not contain at least two energy points: %s' % `I1` )
        for E in Es :
            if( muGridPerE ) : thisMuGrid = muGrids[E]
            if( EMin <= E <= EMax ) : 
                muP = endl2dmathClasses.endl2dmath( I1.getEData( E ) )
            else :
                muP = endl2dmathClasses.endl2dmath( [ [ thisMuGrid[0], 0 ], [ thisMuGrid[-1], 0 ] ] )
            for mu in thisMuGrid : muP.setValue( mu, muP.getValue( mu ) )
            I1.setEData( E, muP )

def toZAsFrameMisc( endlI, old, newProjectileMass, newTargetMass, ELevel, halflife, bdflsFile, firstParameterEIn ) :

    import endl2
    h = copy.deepcopy( old.h )
    new = endlI( None, old.yo, old.C, old.I, old.S, h, points = old.copyData( ), bdflsFile = bdflsFile ) 
    new.setZA( endl2.yoToZA( old.getYi( ) ) )
    new.setYi( endl2.ZAToYo( old.getZA( ) ) )
    new.setMass( newTargetMass )
    new.setELevel( ELevel )
    new.setHalflife( halflife )
    new.setX1( 0. )
    r = newProjectileMass / newTargetMass
    if( firstParameterEIn ) :
        for xy in new.data : xy[0] *= r
    return( new, r )
