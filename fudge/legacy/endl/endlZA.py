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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
This module contains the class endlZA.
"""

ZNameErrorCounter = 0

import os
import tempfile
import shutil
import copy
import re
from fudge import fudgeDefaults, fudgeParameters
from fudge.core.utilities import fudgeExceptions, brb
import endl2dmathClasses
import bdfls
import fudgeDocumentationFile
import endlFile
import endl2
import endlmisc
import endlIClasses
import endl_C
import endl_y
import endl_I
import endl_Z
import endlReactionParameters
try :
    from fudge.structure import xensl
except :
    endlmisc.printWarning( 'Could not import xensl.py' )

class endlZA :
    """This class contains the information for an ENDL za.
    """

    def __init__( self, ZA, yi, database = None, release = 'current', workDir = None, suffix = "", readOnly = 0, bdflsFile = None ) :
        """Creates a new endlZA object."""

        if( bdflsFile == None ) : bdflsFile = bdfls.getDefaultBdfls( )
        self.bdflsFile = bdflsFile
        self.sZA = endlmisc.strZASuffix( ZA, suffix )
        ZA_, suffix_ = endlmisc.intZASuffix( self.sZA )
        self.ZA = ZA_
        self.Z, self.A = endl2.ZandAFromZA( self.ZA )
        self.suffix = suffix_
        self.targetName   = endl2.nameForYoOrZA(   self.sZA, NameASeperator = "", ZAOnly = 1, AddNatural = 1, m_to_m1 = True, suffixSeperator = '_' )
        self.targetSymbol = endl2.endlToGNDName( self.sZA )
        self.name = self.targetSymbol
        self.yiTags = endlmisc.incidentParticleTags( yi )
        self.yi = self.yiTags[0]
        self.projectileName = endl2.endlToGNDName( self.yi )
        self.readOnly = readOnly or fudgeParameters.ReadOnly
        gappsDir = os.path.realpath( '/usr/gapps/data/nuclear' )
        self.documentation = None
        if( database == None ) :                    # This is a designer isotope.
            self.DesignerIsotope = 1
            self.source = "Designer isotope"
            if( workDir == '' ) : workDir = None
            if ( self.readOnly ) :
                self.workDir = None
            else :
                if( workDir == None ) :
                    self.workDir = tempfile.mkdtemp( )
                else :
                    self.workDir = workDir
                if( gappsDir == os.path.realpath( self.workDir )[:len( gappsDir )] ) :
                    if( 'development' not in os.path.realpath( self.workDir ) ) :
                        raise Exception( 'cannot make workDir point to /usr/gapps/data/nuclear, except in a development sub-directory' )
                self.workDir = endlmisc.appendyiDirIfNeeded( self.workDir, self.yi )
                self.workDir = os.path.join( self.workDir, self.sZA )
                if( not os.path.isdir( self.workDir ) ) : os.makedirs( self.workDir )
            if ( fudgeParameters.VerboseMode > 0 ) : print "workDir =", self.workDir, "\nSource directory =", self.source
            self.filelist = []
        else :
            self.DesignerIsotope = 0
            ( self.source, Status ) = endlmisc.getFullPath( database, fudgeDefaults.ENDL_DATABASE_DIR )
            if( Status == 0 ) : ( self.database, Status ) = endlmisc.getFullPath( os.path.join( database, 'ascii' ), \
                os.path.join( fudgeDefaults.TRANSLATED_DATABASE_DIR, release ) )
            if ( Status != 2 ) : raise Exception( "\nError in endlZA.__init__: ZA sub-directory is not a directory or does not exist" )
            if( workDir == '' ) : workDir = self.source
            self.source = endlmisc.appendyiDirIfNeeded( self.source, self.yi )
            self.source = os.path.join( self.source, self.sZA )
            if( os.path.isdir( self.source ) ) :                # Does the ZA directory exist
                if ( self.readOnly ) :
                    self.workDir = None
                else :
                    if( workDir == None ) :
                        self.workDir = tempfile.mkdtemp( )
                    else :
                        self.workDir = workDir
                    if( gappsDir == os.path.realpath( self.workDir )[:len( gappsDir )] ) :
                        if( 'development' not in os.path.realpath( self.workDir ) ) :
                            raise Exception( 'cannot make workDir point to /usr/gapps/data/nuclear, except in a development sub-directory' )
                    self.workDir = endlmisc.appendyiDirIfNeeded( self.workDir, self.yi )
                    self.workDir = os.path.join( self.workDir, self.sZA )
                    if( not os.path.isdir( self.workDir ) ) : os.makedirs( self.workDir )
                if ( fudgeParameters.VerboseMode > 0 ) : print "workDir =", self.workDir, "\nSource directory =", self.source
                F = os.listdir( self.source )
                if( ( self.workDir != None ) and os.path.samefile( self.source, self.workDir ) ) :
                    if( fudgeParameters.VerboseMode > 0 ) : 
                        endlmisc.printWarning( 'Are you sure you know what you are doing as database and work directories are the same?' )
                else :
                    if( not self.readOnly ) :
                        for f in F :
                            fp = os.path.join( self.source, f )
                            try :
                                if( os.path.isfile( fp ) ) : shutil.copyfile( fp, os.path.join( self.workDir, f ) )
                            except :
                                endlmisc.printWarning( 'Could not copy source file %s' % fp )
                F.sort( )
                self.filelist = []
                for f in F :        # Process only yo##c##i###s### files and create list self.filelist[]
                    if( re.match( r"yo\d\dc\d\di\d\d\ds\d\d\d$", f ) ) :
                        if( f[8:11] in [ '030', '032' ] ) : continue        # Skip fluorescence data that is not support in Fudge.
                        self.filelist.append( endlFile.endlFile( 0, f, self.source + os.sep + f, self.ZA, self.yi, bdflsFile = self.bdflsFile ) )
                documentationFileName = os.path.join( self.source, 'documentation.txt' )
                self.documentation = None
                if( os.path.exists( documentationFileName ) ) : self.documentation = fudgeDocumentationFile.fudgeDocumentationFile( documentationFileName )
            else :
                raise Exception( "\nError in endlZA.__init__: database %s does not contain isotope = %d" % ( self.source, ZA_ ) )

    def __getattr__( self, name ) :

        global ZNameErrorCounter
        if( name in [ 'ZName', 'ZSymbol' ] ) :
            if( ZNameErrorCounter == 0 ) :
                print "\n ------------------------- Warning -------------------------"
                print "  You have just used '%s' which is being deprecated in version 3." % name
                print "  ZName and ZSymbol are deprecated members of endlZA. Please"
                print "  use targetName and targetSymbol instead. This message is only"
                print "  printed once."
                print "\n  Traceback information to follow:"
                import traceback
                traceback.print_stack( )
                print "\n -----------------------------------------------------------"
            ZNameErrorCounter += 1
            if( name == 'ZName' ) : return( self.targetName )
            return( self.targetSymbol )
        raise AttributeError( "endlZA instance has no attribute '%s'" % name )

    def __getitem__( self, index ) :
        """Calls findData( yo, C, I, S, X1, X2, X3, X4, Q ) and returns the matching item.
For example, if za is an endlZA object then either of the two methods below 
can be used to get a reference to the neutron in-elastic cross section,
>>> i1 = za[ ( 1, 11, 1, 0 ) ]
>>> i1 = za.findData( 1, 11, 1, 0 )
"""

        if ( len( index ) < 3 ) :
            raise Exception( "\nError in endlZA.__getitem__: index must contain at least yo, C and I, ( yo, C, I, S, X1, X2, X3, X4, Q )" )
        elif ( len( index ) > 9 ) :
            raise Exception( "\nError in endlZA.__getitem__: index must contain only ( yo, C, I, S, X1, X2, X3, X4, Q )" )
        yo = index[0]
        C = index[1]
        I = index[2]
        S_XnQ = [ None, None, None, None, None, None ]
        i = 0
        while ( len( index ) > i + 3 ) :
            S_XnQ[i] = index[i + 3]
            i += 1
        return self.findData( yo, C, I, S_XnQ[0], S_XnQ[1], S_XnQ[2], S_XnQ[3], S_XnQ[4], S_XnQ[5] )

    def addEndlFile( self, endlfile, replace = 0, checkHeaders = 1, halflife = None, templateData = None ) :
        """Adds endlFile instance endlfile to self.  If replace is 0 and file exist in self, then
        a raise is executed.  If replace is 1 and file exist, a warning message is printed
        and endlfile replaces the existing file.  If replace is niether 0 or 1, any existing
        file in self is replaced by endlfile without any notification.  If checkHeaders is true,
        the headers of endlfile are compared to other headers in self, and modified if needed.
        A reference to endlfile is returned."""

        if( templateData == None ) :
            for f in self.filelist :
                if( f.C != endlfile.C ) : continue
                for l in f.levels :
                    if( len( l.h[0] ) >= 68 ) :
                        templateData = l
                        break
            if( templateData == None ) :
                for f in self.filelist :
                    for l in f.levels :
                        if( len( l.h[0] ) >= 68 ) :
                            templateData = l
                            break
        i = 0
        if( replace ) :
            for f in self.filelist :
                if ( f.name == endlfile.name ) :
                    if( replace ) : endlmisc.printWarning( "Warning in endlZA.addEndlFile: file exist and is being replaced" )
                    del self.filelist[i]
                    break
                i += 1
        f = self.addFile( endlfile.yo, endlfile.C, endlfile.I, endlfile.S, halflife = halflife )
        if( endlfile.levels == [] ) : endlfile.read( )
        for l in endlfile.levels : f.addEndlData( l )
        if( checkHeaders ) :
            if( templateData != None ) :
                if( ( l.C != 10 ) and ( l.C != templateData.C ) ) : endlmisc.printWarning( "Warning in endlZA.addEndlFile: template reaction not the same as new data, so Q and temperature may not be correct." )
                for l in f.levels :
                    l.setZA( templateData.getZA( ) )
                    l.setMass( self.bdflsFile.mass( l.getZA( ) ) )
                    l.setYi( templateData.getYi( ) )
                    l.setHalflife( templateData.getHalflife( ) )
                    l.setTemperature( templateData.getTemperature( ) )
                    if( l.C == templateData.C ) : l.setQ( templateData.getQ( ) )
#                    if( len( l.h[0] ) > 34 ) :
#                        l.h[0] = templateData.h[0][:25] + l.h[0][25:34] + templateData.h[0][34:]
#                    else :
#                        l.h[0] = copy.deepcopy.templateData.h[0]
#                    l.setYo( endlfile.yo )
            else :
                if( l.C != 10 ) : endlmisc.printWarning( "Warning in endlZA.addEndlFile: no template data found, so Q and temperature may not be correct." )
                for l in f.levels :
                    l.setZA( f.ZA )
                    l.setMass( self.bdflsFile.mass( f.ZA ) )
                    l.setYi( f.yi )
                    l.setHalflife( self.bdflsFile.halflife( f.ZA ) )
        return f

    def addFile( self, yo, C, I, S, halflife = None, printWarnings = True ) :
        "Adds the (yo, C, I, S) file to endlZA.  If the file already exist a raise is executed."

        if ( ( ( I == 10 ) or ( I == 11 ) ) and printWarnings ) :
            endlmisc.printWarning( "Warning in endlZA.addFile: I-values of 10 or 11 may be overwritten in processing" )
        if ( ( I == 0 ) and ( yo != 0 ) ) : raise Exception( "\nError in endlZA.addFile: For I = 0, yo must also be 0" )
        i = 0
        c = 1
        name = "yo%2.2dc%2.2di%3.3ds%3.3d" % ( yo, C, I, S )
        for f in self.filelist :
            if ( f.name == name ) : raise fudgeExceptions.ENDL_addFileException_FileExist( "Error in endlZA.addFile: file already exist" )
            if ( f.name  > name ) : c = 0
            if ( c ) : i += 1
        f = endlFile.endlFile( 1, name, "Designer File", self.ZA, self.yi, halflife = halflife, bdflsFile = self.bdflsFile )
        self.filelist.insert( i, f )
        return( f )

    def changeFileYoCIS( self, file, yo, C, I, S, deleteOldFile = True ) :
        """Changes endlFile object file's yo-, C-, I- and S-values to yo, C, I and S.
        If deleteOldFile is true, the old file is deleted (a new file is not created unless method save is called)."""

        if( ( file.yo == yo ) and ( file.C == C ) and ( file.I == I ) and ( file.S == S ) ) : return
        oldName = file.name
        file.setYiZAYoCISs( yo = yo, C = C, I = I, S = S )
        if( deleteOldFile ) :
            cmd = "rm -rf " + self.workDir + "/" + oldName
            if ( self.workDir != None ) : 
                if ( fudgeParameters.VerboseMode > 0 ) : print "executing system command >>>", cmd
                os.system( cmd )

    def check( self, dQ_MeV = 10e-3, dThreshold_MeV = 10e-3, xCloseEps = None, checkForEnDepData = False, normTolerance = 1e-5, allowZeroE = False,
        checkList = {}, crossSectionsOnly = False, maxAbsFloatValue = None ) :
        "Checks all read in data in self."

        checkListKeys = { 'checkEMin' : True, 'checkMissing' : True, 'checkSubChecks' : True, 'energyBalance' : 10e-3, 'maxGammaMultiplicity' : 100.,
            'crossSectionEnergyMax' : 20 }
        for checkKey in checkList : checkListKeys[checkKey] = checkList[checkKey]
        maxGammaMultiplicity = checkListKeys['maxGammaMultiplicity']
        ErrMsgs = []
        targetParameters = { 'T' : None, 'ELevel' : None, 'halflife' : None }
        for f in self.filelist :
            if( f.levels == [] ) : f.read( )
        if( checkListKeys['checkSubChecks'] ) :
            for f in self.filelist : 
                allowZeroE_ = allowZeroE
                if( f.I in [ 941, 942 ] ) : allowZeroE_ = True
                ErrMsgs += f.check( targetParameters = targetParameters, dQ_MeV = dQ_MeV, dThreshold_MeV = dThreshold_MeV, xCloseEps = xCloseEps,
                    allowZeroE = allowZeroE_, normTolerance = normTolerance, maxAbsFloatValue = maxAbsFloatValue )
        reactionsDatas = self.findReactionsDatas( )
        crossSectionEnergyMax = checkListKeys['crossSectionEnergyMax']
        if( type( crossSectionEnergyMax ) != type( None ) ) :
            crossSectionEnergyMax = float( crossSectionEnergyMax )
            I0s = self.findDatas( I = 0 ) 
            for I0 in I0s :
                if( I0.xMax( ) < crossSectionEnergyMax ) :
                    if( I0.getValue( I0.xMax( ) ) != 0 ) :
                        ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = I0.yo, C = I0.C, S = I0.S, I = 0, X1 = I0.X1, \
                            X2 = I0.X2, X3 = I0.X3, X4 = I0.X4, Q = I0.Q, message = 'non zero termination of cross section data below EMax' ) )
        for reactionDatas in reactionsDatas :
            I0 = reactionDatas[0]
            C, S, X1, X2, X3, X4, Q = I0.C, I0.S, I0.X1, I0.X2, I0.X3, I0.X4, I0.Q
            if( I0.I != 0 ) :
                ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 0, C = C, S = S, I = 0, X1 = X1, \
                    X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'missing I = 0 data' ) )
                continue
            energyDepedentQ = None
            if( C == 15 ) :
                promptAndDelayed = []
                for data in reactionDatas :
                    if( data.I == 7 ) : promptAndDelayed.append( data.S )
                    if( data.I == 12 ) : energyDepedentQ = data
                if( 0 not in promptAndDelayed ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 1, C = C, S = 0, \
                    I = 7,  X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'missing I = 7, S = 0 data' ) )
                if( 7 not in promptAndDelayed ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 1, C = C, S = 7, \
                    I = 7,  X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'missing I = 7, S = 7 data' ) )
                if( 17 in promptAndDelayed ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 1, C = C, S = 17, \
                    I = 7,  X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'unsupported I = 7, S = 17 data found' ) )
            EMinI0 = I0.xMin( )
            if( checkListKeys['checkEMin'] ) :
                for data in reactionDatas :
                    if( data.I < 20 ) :             # 20 is URR probability tables.
                        if( data.columns == 4 ) :
                            E = data.EMin( )
                        else :
                            E = data.xMin( )
                        if( E < ( EMinI0 - dThreshold_MeV ) ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, data = data, \
                            message = "first E value = %e is less than I = 0's first E value = %e" % ( E, EMinI0 ) ) )
                        if( E > ( EMinI0 + dThreshold_MeV ) ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, data = data, \
                            message = "first E value = %e is greater than I = 0's first E value = %e" % ( E, EMinI0 ) ) )
            yoCounts = 10 * [0]
            yos = []
            if( I0.C == 15 ) :
                yoCounts[1] = -1
                yos = [ 1 ]
                yosRes = [ 1 ]
            else :
                residualZA, yosTuple, QDummy = endl2.residualZA_yos_Q( self.yi, self.ZA, C, printQWarning = False, bdflsFile = self.bdflsFile )
                yos += yosTuple                     # This is done because a list, [], is wanted and not a tuple, ().
                if( len( yos ) > 0 ) :
                    if( yos[0] < 0 ) :
                        if( yos[0] == -3 ) :        # Fission
                            yos[0] = 1
                        elif( yos[0] == -4 ) :      # many yos[1]
                            yos = [ yos[1] ]
                    for i in xrange( len( yos ) ) : yos[i] = endl2.ZAToYo( yos[i] )
                yosRes = copy.deepcopy( yos )
                try :
                    iyo = endl2.ZAToYo( residualZA )
                    yosRes.append( iyo )
                except :
                    pass
                for yo in yos :
                    if( yo >= 0 ) :
                        iyo = endl2.ZAToYo( yo )
                        yoCounts[iyo] += 1
                for data in reactionDatas :
                    yo = data.yo
                    if( 9 < yo < 20 ) : yo -= 10 
                    if( yo == 0 ) : continue
                    if( yo not in yosRes ) :
                        if( ( yo == 7 ) and ( 10 < C < 50 ) ) : continue
                        ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, data = data, \
                            message = "yo = %d invalid for this reaction" % data.yo ) )
            yosTotal = 0
            for yo in xrange( len( yoCounts ) ) :
                yosTotal += yoCounts[yo]
                if( yoCounts[yo] != 0 ) : yoCounts[yo] = yo            # Now if yo present, set yoCounts[yo] to yo's ID.
            try :
                iyo = endl2.ZAToYo( residualZA )
                yoCounts[0] = iyo + 10
            except :
                yoCounts[0] = 0
            for data in reactionDatas :
                if( data.yo not in yosRes ) :
                    if( ( data.yo == 7 ) and ( 10 < C < 50 ) and ( data.I != 10 ) ) : yoCounts[7] = 7
            for yo in xrange( 20 ) :
                datas = self.findDatas( yo = yo, C = C, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )
                Is = {}
                for data in datas :
                    if( data.I not in Is ) : Is[data.I] = 0
                    Is[data.I] += 1
                for I in Is :
                    if( Is[I] > 1 ) : 
                        datas = self.findDatas( yo = yo, C = C, S = S, I = I, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )
                        ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, data = datas[0], message = "multiple data found" ) )
            totalEnergyDeposition = endl2dmathClasses.endl2dmath( [] )
            energyDepositionToGammasPresent = False
            for yo in yoCounts :                                        # yoCounts[0] is the residual yo if one exists.
                if( yo > 0 ) and not crossSectionsOnly :
                    yoMsg = None
                    datas = self.findDatas( yo = yo, C = C, I = ( 1, 3, 4, 9 ), S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )
                    lenI4 = 0
                    Is = []
                    for data in datas :
                        if( data.I == 4 ) : lenI4 = len( data )
                        Is.append( data.I )
                    I1present = 1 in Is
                    I3present = 3 in Is
                    I4present = 4 in Is
                    I9present = 9 in Is
                    if( C == 46 and yo == 7 and I1present and ( not I3present and not ( I4present and lenI4 == 1 ) ) ) :
                        yoMsg = 'missing I = 3 or 4 (l=0) data to go with I = 1 data'
                    elif( ( C != 10 ) and ( S == 0 ) and I1present and not( I3present or I4present ) ) :
                        yoMsg = 'missing I = 3 or 4 data to go with I = 1 data'
                    elif( I3present ) :
                        if( not I1present ) :
                             yoMsg = 'missing I = 1 data to go with I = 3 data'
                        else :
                            I1 = self.findData( yo = yo, C = C, I = 1, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, ifMultipleReturnFirst = True )
                            I3 = self.findData( yo = yo, C = C, I = 3, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, ifMultipleReturnFirst = True )
                            I1XArray = I1.xArray( )
                            I3TArray = I3.tArray( )
                            if( I1XArray != I3TArray ) :
                                yoMsg = 'mismatch in I = 1 and I = 3 incident energies'
                            else :
                                for i, EmuP1 in enumerate( I1.data ) :
                                    E1, muP1 = EmuP1
                                    E3, muEpP3 = I3.data[i]
                                    if( len( muP1 ) != len( muEpP3 ) ) :
                                        yoMsg = 'mismatch in I = 1 and I = 3 number of mus for E = %e' % E1
                                    else :
                                        for j, muP in enumerate( muP1 ) :
                                            if( muP[0] != muEpP3[j][0] ) :
                                                yoMsg = 'mismatch in I = 1 and I = 3 mus for E = %e' % E1
                                                break
                    elif( not I4present ) :
                        okay =          not ( ( self.ZA  == 3006 ) and ( yo == 1 ) and ( C == 22 ) ) # Special case for break up of residual in ENDL
                        okay = okay and not ( ( self.ZA  == 6012 ) and ( yo == 1 ) and ( C == 27 ) )
                        if( I0.S == 1 ) :
                            if( yo in yos ) : okay = False
                        if( ( yosTotal > 1 ) and okay ) :               # yosTotal does not count the residual, so 2-body must have yosTotal = 1.
                            yoMsg = 'missing I = 3 or 4 data'
                        elif( not I1present ) :
                            yoMsg = 'missing I = 1 data'
                        else :
                            if( ( self.A == 0 ) and ( ( I0.C < 50 ) or ( I0.C > 59 ) ) and ( self.yi != yo ) ) : 
                                yoMsg = 'natural target cannot have only I = 1 data for this reaction'
                    if( yoMsg != None ) :
                        ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = yo, C = C, S = S, X1 = X1, X2 = X2, \
                            X3 = X3, X4 = X4, Q = Q, message = '%s for yo = %2d' % ( yoMsg, yo ) ) )
                    if( I3present and I4present ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = yo, C = C, S = S, X1 = X1, \
                            X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'I = 3 and 4 present for yo = %2d' % yo ) )
                    if( ( 10 < C < 50 ) and ( yo == 7 ) and ( not I9present ) ) :
                        ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = yo, C = C, I = 9, S = S, X1 = X1, \
                            X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'missing I = 9 data for yo =  7' ) )
                    if( checkForEnDepData ) :
                        try :
                            data = self.findData( yo = yo, C = C, I = 10, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, ifMultipleReturnFirst = True )
                            totalEnergyDeposition = totalEnergyDeposition + data
                            if( yo == 7 ) : energyDepositionToGammasPresent = True
                        except :
                            ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = yo, C = C, I = 10, S = S, X1 = X1, 
                                X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'missing I = 10 for yo = %2d' % yo ) )
            gammaMultiplicities = self.findDatas( yo = 7, C = C, I = 9, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )
            for gamma in gammaMultiplicities :
                if( data.yMax( ) > maxGammaMultiplicity ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 7, C = C, I = 9, 
                    S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'large gamma multiplicity = %d' % data.yMax( ) ) )
            if( checkForEnDepData and ( yoCounts[0] == 0 ) and ( ( C < 50 ) or ( C > 57 ) ) ) :
                try :
                    data = self.findData( yo = 0, C = C, I = 11, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, ifMultipleReturnFirst = True )
                    totalEnergyDeposition = totalEnergyDeposition + data
                except :
                    ErrMsgs.append( endlmisc.endlCheckerObject( yi = self.yi, ZA = self.ZA, suffix = self.suffix, yo = 0, C = C, I = 11, S = S, X1 = X1, X2 = X2, \
                            X3 = X3, X4 = X4, Q = Q, message = 'missing I = 11 for yo = %2d' % yo ) )
            if( checkForEnDepData and ( ( C < 50 ) or ( C > 57 ) ) ) :
                if( type( energyDepedentQ ) == type( None ) ) : energyDepedentQ = endl2dmathClasses.endl2dmath( [ [ I0.xMin( ), Q ], [ I0.xMax( ), Q ] ] )
                incidentEnergy = endl2dmathClasses.endl2dmath( [ [ I0.xMin( ), I0.xMin( ) ], [ I0.xMax( ), I0.xMax( ) ] ] )
                dELevels = 0.
                if( S == 1 ) :
                    if( ( C in [ 11, 40, 41, 42, 44, 45 ] ) and not energyDepositionToGammasPresent ) : dELevels -= I0.getX1( )
                dE = totalEnergyDeposition - ( incidentEnergy + energyDepedentQ + dELevels )
                dEMax = float( checkListKeys['energyBalance'] )
                i = -1
                for E, deltaE in dE.data :
                    i += 1
                    if( abs( deltaE ) > dEMax ) : ErrMsgs.append( endlmisc.endlCheckerObject( yi = 0, ZA = self.ZA, suffix = self.suffix, yo = 0, C = C, I = None, S = S, \
                        X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q, message = 'energy balance issue at E(%6d) = %e (dE = %e)' % ( i, E, deltaE ) ) )
        ErrMsgs.sort( )
        i = 0
        j = 1
        sortOrder = endlmisc.endlCheckerObject.sortOrder
        endlmisc.endlCheckerObject.sortOrder = [ 'database', 'yi', 'ZA', 'suffix', 'C', 'S', 'X1', 'yo', 'Q', 'I', 'X2', 'X3', 'X4' ]
        while( j < len( ErrMsgs ) ) :
            if( ErrMsgs[i] == ErrMsgs[j] ) :
                if( type( ErrMsgs[i].message ) != type( [] ) ) : ErrMsgs[i].message = [ ErrMsgs[i].message ]
                if( type( ErrMsgs[j].message ) == type( [] ) ) :
                    ErrMsgs[i].message += ErrMsgs[j].message
                else :
                    ErrMsgs[i].message.append( ErrMsgs[j].message )
                del ErrMsgs[j]
            else :
                i += 1
                j = i + 1
        endlmisc.endlCheckerObject.sortOrder = sortOrder
        return( ErrMsgs )

    def clean( self ) :
        "Removes all I=10 and I=11 data files from self."

        r = range( len( self.filelist ) )
        r.reverse( )
        for i in r :
            if ( self.filelist[i].I == 10 ) or ( self.filelist[i].I == 11 ) : del self.filelist[i]

    def CList( self ) :
        "Returns the list of sorted C-values for self."

        cl = []
        for l in self.filelist :
            if( l.C not in cl ) : cl.append( l.C )
        cl.sort( )
        return( cl )

    def convertS7ForOldProcessing( self ) :
        """This method removes all S = 7, fission delayed nuclear data, data except I = 7 data 
    and collapses the I = 7 data to a single data set, since ndfgen cannot handle non I = 7, 
    S = 7 data and can only handle a single I = 7, S = 7 data set."""

        delayedFileList = self.files( C = 15, S = 7 )
        for delayedFile in delayedFileList : 
            if( delayedFile.I != 7 ) : self.removeFile( delayedFile.yo, delayedFile.C, delayedFile.I, delayedFile.S )
        delayedFileList = self.files( C = 15, S = 7, I = 7 )
        if( len( delayedFileList ) == 1 ) :
            delayedFile = delayedFileList[0]
            if( delayedFile.levels == [] ) : self.read( C = 15, S = 7, I = 7 )
            datas = self.findDatas( C = 15, S = 7, I = 7 )
            if( len( datas ) > 1 ) :
                firstData = True
                for data in datas :
                    if( firstData ) :
                        delayedData = data
                        firstData = False
                    else :
                        delayedData = delayedData + data
                        delayedFile.unreadData( X1 = data.X1 )
                datas[0].set( delayedData )
                datas[0].setX1( 0. )
        elif( len( delayedFileList ) > 1 ) :
            raise Exception( "\nError in endlZA.convertS7ForOldProcessing: multiple I = 7, S = 7 data files." )

    def convertI4ForOldProcessing( self, nMu = 21, lMax = None, skipIf_deuteron_S8 = True ) :
        """This method finds all I = 4 data with Legendre order greater than 0 in self and converts
    it into I = 1 and 3 data so the mcfgen can process it.  See method convertToUI3 in class
    endlI4 for meaning of nMu and lMax.  If skipIf_deuteron_S8 is true, than for ZA == 1002 the
    S = 8 data is not converted as mcfgen was designed to ignore this data in ENDL94 and its
    descendent data sets."""

        I4s = self.findDatas( I = 4 )
        I4sToRemove = []
        for I4 in I4s :
            if( ( skipIf_deuteron_S8 ) and ( self.ZA == 1002 ) and ( I4.S == 8 ) ) : continue
            if( len( I4 ) > 1 ) :
                if ( fudgeParameters.VerboseMode > 0 ) : print 'converting yo = %2d, C = %2d, S = %3d to I = 1 and 3' % ( I4.yo, I4.C, I4.S )
                I3 = I4.convertToUI3( nMu = nMu, lMax = lMax )
                I1 = I3.reduceToEMuP( normalize = True )
                I3.normalize( )

                I3File = self.findFile( yo = I4.yo, C = I4.C, I = 3, S = I4.S )
                if( I3File == None ) : I3File = self.addFile( yo = I4.yo, C = I4.C, I = 3, S = I4.S )
                I3File.addEndlData( I3 )

                I1File = self.findFile( yo = I4.yo, C = I4.C, I = 1, S = I4.S )
                if( I1File == None ) : I1File = self.addFile( yo = I4.yo, C = I4.C, I = 1, S = I4.S )
                I1File.addEndlData( I1 )

                yoCIS = [ I4.yo, I4.C, I4.I, I4.S ]
                if( yoCIS not in I4sToRemove ) : I4sToRemove.append( yoCIS )

        for yo, C, I, S in I4sToRemove : self.removeFile( yo, C, I, S )

    def SListForC( self, C ) :
        "Returns the list of S-values for C-value of self."

        sl = []
        for l in self.filelist :
            if( ( l.C == C ) and ( l.S not in sl ) ) : sl.append( l.S )
        return( sl )

    def cs( self, C, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None ) :
        """Calls findData( 0, C, 0, S, X1, X2, X3, X4, Q ) and returns the item. Deprecated method, will not be in version 4."""

        return self.findData( 0, C, 0, S, X1, X2, X3, X4, Q )

    def file( self, yo = None, C = None, I = None, S = None ) :
        "Deprecated method, use findFile instead. This method will not be in version 4."

        return( self.findFile( yo = yo, C = C, I = I, S = S ) )

    def findFile( self, yo = None, C = None, I = None, S = None ) :
        """Returns reference to file matching yo, C, I and S. If no file matches, None is returned. If multiple
        files match, an FUDGE_Exception is raised."""

        a = self.findFiles( yo, C, I, S )
        if ( a == [] ) : return( None )
        if ( len( a ) > 1 ) :
            endlmisc.printWarning( "  yo   C   I   S" )
            for f in a : endlmisc.printWarning( " %3d %3d %3d %3d" % ( f.yo, f.C, f.I, f.S ) )
            raise fudgeExceptions.FUDGE_Exception( "Error in endlZA.file: multiple files found (see 'yo   C   I   S' list above)" )
        return a[0]

    def files( self, yo = None, C = None, I = None, S = None ) :
        "Deprecated method, use findFiles instead. This method will not be in version 4."

        return( self.findFiles( yo = yo, C = C, I = I, S = S ) )

    def findFiles( self, yo = None, C = None, I = None, S = None ) :
        """Returns a list of files matching yo, C, I and S"""

        a = []
        for f in self.filelist :
            if ( f.IsyoCIS( yo, C, I, S ) ) : a.append( f )
        return a

    def findItems( self, yo = None, C = None, I = None, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None ) :
        "Deprecated method, use findDatas instead. This method will not be in version 4."

        return( self.findDatas( yo = yo, C = C, I = I, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q ) )

    def findDatas( self, yo = None, C = None, I = None, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None ) :
        """Returns a list of items matching yo, C, I, S, X1, X2, X3, X4 and Q.  yo, C, I and S
        can be a list of integers (e.g., findDatas( C = [ 10, 11 ] ). The value None means to match all."""

        items = []
        for f in self.filelist : 
            if ( f.IsyoCIS( yo, C, I, S ) ) :
                if ( f.levels == [] ) : raise Exception( "\nError in endlZA.findDatas: one or more data files not read" )
                for l in f.levels :
                    if ( ( ( X1 == None ) or ( X1 == l.X1 ) ) and
                         ( ( X2 == None ) or ( X2 == l.X2 ) ) and
                         ( ( X3 == None ) or ( X3 == l.X3 ) ) and
                         ( ( X4 == None ) or ( X4 == l.X4 ) ) and
                         ( ( Q  == None ) or ( Q  == l.Q  ) ) ) : items.append( l )
        return( items )

    def findItem( self, yo = None, C = None, I = None, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None ) :
        "Deprecated method, use findData instead. This method will not be in version 4."

        return( self.findData( yo = yo, C = C, I = I, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q ) )

    def findData( self, yo = None, C = None, I = None, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None, ifMultipleReturnFirst = False ) :
        """Returns a single item matching yo, C, I, S, X1, X2, X3, X4 and Q."""

        f = self.findDatas( yo, C, I, S, X1, X2, X3, X4, Q )
        if ( f == [] ) :
            s = "\nError in endlZA.findData: no match for (yo, C, I, S, X1, X2, X3, X4, Q) = (%s, %s, %s, %s, %s, %s, %s, %s, %s)" % \
                ( `yo`, `C`, `I`, `S`, `X1`, `X2`, `X3`, `X4`, `Q` )
            raise Exception( s )
        elif ( ( len( f ) > 1 ) and ( not ifMultipleReturnFirst ) ) :
            print "The following values missing "
            print "\nList of (yo, C, I, S, X1, X2, X3, X4, Q) values"
            print "     yo   C    I   S         X1             X2             X3             X4             Q"
            missing = [ 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
            q_0 = f[0]
            for q in f : 
                missing = [ missing[0] or ( q_0.yo != q.yo ), missing[1] or ( q_0.C != q.C ), 
                    missing[2] or ( q_0.I != q.I ), missing[3] or ( q_0.S != q.S ), missing[4] or ( q_0.X1 != q.X1 ),
                    missing[5] or ( q_0.X2 != q.X2 ), missing[6] or ( q_0.X3 != q.X3 ), missing[7] or ( q_0.X4 != q.X4 ),
                    missing[8] or ( q_0.Q != q.Q ) ]
                print "    (%2d, %2d, %3d, %2d, %14.6e %14.6e %14.6e %14.6e %14.6e)" % ( q.yo, q.C, q.I, q.S, q.X1, q.X2, q.X3, q.X4, q.Q )
            print "\nValues not specified are -",
            if( yo == None ) : print "yo -",
            if(  C == None ) : print "C -",
            if(  I == None ) : print "I -",
            if(  S == None ) : print "S -",
            if( X1 == None ) : print "X1 -",
            if( X2 == None ) : print "X2 -",
            if( X3 == None ) : print "X3 -",
            if( X4 == None ) : print "X4 -",
            if( Q  == None ) : print "Q -",
            print "\nAdditional values required are -", 
            if( missing[0] ) : print "yo -",
            if( missing[1] ) : print "C -",
            if( missing[2] ) : print "I -",
            if( missing[3] ) : print "S -",
            if( missing[4] ) : print "X1 -",
            if( missing[5] ) : print "X2 -",
            if( missing[6] ) : print "X3 -",
            if( missing[7] ) : print "X4 -",
            if( missing[8] ) : print "Q -",
            print "\n"
            raise Exception( "\nError in endlZA.findData: multiple objects found matching (yo, C, I, S, X1, X2, X3, X4, Q) (see above Traceback)" )
        return( f[0] )

    def findReactionsDatas( self ) :
        """Returns a list of reactionDatas for each reaction in self where reactionDatas is a list of
        all data for a reaction."""

        reactionsDatas = []
        reactions = self.reactionList( )
        for reaction in reactions :
            datas = self.findDatas( C = reaction.C, S = reaction.S, X1 = reaction.X1, X2 = reaction.X2, X3 = reaction.X3, X4 = reaction.X4, Q = reaction.Q )
            if( datas[0].C == 15 ) :
                datas += self.findDatas( C = 15, S = 7, X2 = reaction.X2, X3 = reaction.X3, X4 = reaction.X4, Q = reaction.Q )
            Is = []
            for data in datas : Is.append( [ data.yo, data.I, data.S, data ] )
            Is.sort( )
            datas = []
            for yo, I, S, data in Is : datas.append( data )
            reactionsDatas.append( datas )

        return( reactionsDatas )

    def fixHeaders( self, yi = None, ZA = None, ELevel = None, temperature = None, dTemperatureEps = 1e-4 ) :
        """This method makes the headers for all the data in self be consistent (e.g., makes all the
        temperatures the same). Items the are made consistent are ZA, yi, mass, ELevel, halflife, temperature
        amoung all the data and yo, C, I and S within a file. yi, ZA, ELevel and temperature can be inputted.
        For input parameters of None, the C = 10, I = 0 file values are used. If this file does not exist,
        a raise is executed, Mass and halflife are retrieved from the bdfls file."""

        def getAnyData( anyData ) :
            """Attempts to find a dataset that is complete enough to determine. If one is not found a raise is executed."""

            if( anyData != None ) : return( anyData )
            anyDatas = self.findDatas( )
            if( len( anyDatas ) > 0 ) : return( anyDatas[0] )
            raise fudgeExceptions.ENDL_fixHeadersException_NoC10Data( 'No dataset is complete enough for determining required information.' )

        anyData = None
        self.readIfNotRead( )
        if( ( yi == None ) or (  ZA == None ) or ( ELevel == None ) or ( temperature == None ) ) :
            anyData = getAnyData( anyData )
            if( yi == None ) : yi = anyData.getYi( )
            if( ZA == None ) : ZA = anyData.getZA( )
            if( ELevel == None ) : ELevel = anyData.getELevel( )
            if( temperature == None ) : temperature = anyData.getTemperature( )
        mass = self.bdflsFile.mass( ZA )
        if( mass == None ) :
            anyData = getAnyData( anyData )
            mass = anyData.getMass( )
        halflife = self.bdflsFile.halflife( ZA )
        if( halflife == None ) :
            try :
                anyData = getAnyData( anyData )
                halflife = anyData.getHalflife( )
            except :
                halflife = 0.
                if( AFromZA( self.ZA ) == 0 ) : halflife = 1e50
        datas = self.findDatas( )
        for data in datas :
            data.setMass( mass )
            data.setELevel( ELevel )
            if( halflife == None ) :
                halflife = data.getHalflife( )
            else :
                data.setHalflife( halflife )
            data.setTemperature( temperature )
        files = self.findFiles( )
        for file in files : file.setYiZAYoCISs( yi = yi, ZA = ZA )

    def fixupReactions( self, fixReactions = None, fixReactionsEps = 1e-3, fixReactionsQuery = True ) :
        """If the parameter fixReactions is one (or a list containing one or more) of the following 
        'Q', 'X1, 'X2, 'X3' or 'X4', this method will try to find data sets which are not considered 
        to be for the same reaction due to slightly different value(s) in the members as given by 
        fixReactions, but are within fixReactionsEps (e.g., fixReactions = ['Q', 'X1']). If 
        fixReactionsQuery is true, you will be asked for your option on reaction parameters that 
        this method deems the same."""

        def compareValue( v1, v2, doIt, eps ) :
            """For internal use only."""

            if( doIt ) :
                theSame = abs( v2 - v1 ) <= ( abs( v1 ) + abs( v2 ) ) * eps
            else :
                theSame = v1 == v2
            return( theSame )

        if( fixReactions != None ) :
            doIt = 7 * [ False ]
            doIts = { 'X1' : 2, 'X2' : 3, 'X3' : 4, 'X4' : 5, 'Q' : 6 }
            for r in fixReactions :
                if( r not in [ 'Q', 'X1', 'X2', 'X3', 'X4' ] ) : 
                    raise Exception( '\nError endlZA.fixupReactions: %s not a valid fixReactions option' % r )
                doIt[doIts[r]] = True
            reactions = self.reactionList( )
            COld = None
            oldIndex = 0
            n = len( reactions )
            for reaction in reactions :
                oldIndex += 1
                index = oldIndex
                while( index < n ) :
                    nextReaction = reactions[index]
                    if( ( nextReaction.C == reaction.C ) and ( nextReaction.S == reaction.S ) ) : # Same C and S.
                        theSame = True
                        for i in [ 'X1', 'X2', 'X3', 'X4', 'Q' ]:
                            theSame = theSame and compareValue( reaction[i], nextReaction[i], doIt[doIts[i]], fixReactionsEps )
                        if( theSame ) :
                            if( fixReactionsQuery ) :
                                while( True ) :
                                    print 'Reaction set one:'
                                    datas1 = self.findDatas( C = reaction.C, S = reaction.S, Q = reaction.Q, \
                                        X1 = reaction.X1, X2 = reaction.X2, X3 = reaction.X3, X4 = reaction.X4 )
                                    for data in datas1 : print data
                                    print
                                    print 'Reaction set two:'
                                    datas2 = self.findDatas(C = nextReaction[0], S=nextReaction[1], Q = nextReaction[6], \
                                        X1 = nextReaction[2], X2=nextReaction[3], X3=nextReaction[4], X4=nextReaction[5])
                                    for data in datas2 : print data
                                    response = brb.Pause( 'Are they the same (y/n)? : [y] ' )
                                    if( response in [ '', 'n', 'y' ] ) :
                                        if( response == 'n' ) : theSame = False
                                        break
                                    else :
                                        print 'Please respond with y or n'
                            else :
                                datas2 = self.findDatas(C = nextReaction[0], S=nextReaction[1], Q = nextReaction[6], \
                                    X1 = nextReaction[2], X2=nextReaction[3], X3=nextReaction[4], X4=nextReaction[5])
                        if( theSame ) :
                            for data in datas2 :
                                data.setX1( reaction.X1 )
                                data.setX2( reaction.X2 )
                                data.setX3( reaction.X3 )
                                data.setX4( reaction.X4 )
                                data.setQ( reaction.Q )
                                reactions[index] = endlReactionParameters.endlReactionParameters( 0, 0, 0, 0, 0, 0, 0 )
                    else :
                        break
                    index += 1
        
    def fixQandThresholds( self, specialCases = 1, fixThresholdMode = endlIClasses.fixThresholdMode_RaiseOnly, dThreshold_MeV = 10e-3, EFloor = 1e-11,
        threshold_MeV_shiftWarning = 0.1 ) :
        """This methods fixes all the Qs and thresholds for all the data in self. Also, this method makes sures 
        that the first energy point is greater than or equal to EFloor for all data. If the threshold is moved 
        by more than threshold_MeV_shiftWarning, a warning is printed."""

        try :
            C1 = self.findData( C = 1 )
            C1.setQ( 0. )
        except :
            pass
        I0s = self.findDatas( I = 0 )

        S1Cs = []
        for I0 in I0s :
            if( I0.S == 1 ) : S1Cs.append( I0.C )

        for I0 in I0s :
            if( ( I0.C < 8 ) or ( 49 < I0.C < 58 ) ) : continue
            Ss = [ I0.S ]
            if( I0.C == 15 ) : Ss.append( 7 )
            datas = self.findDatas( C = I0.C, S = Ss, X1 = I0.X1, X2 = I0.X2, X3 = I0.X3, X4 = I0.X4, Q = I0.Q )

            yos = endl_C.endl_C_yoInfo( I0.C )                                  # Attempt to handle charged particle "Coulomb barrier" threshold.
            CoulombBarrier = False
            if( len( yos ) > 0 ) :
                if(   yos[0] > 1000 ) :                                         # Charged particles if first out.
                    CoulombBarrier = True
                elif( ( yos[0] == endl_C.yoX ) and ( yos[1] > 1000 ) ) :        # nX for a charged particle.
                    CoulombBarrier = True
                elif( ( yos[0] == endl_C.yo_eq_yi ) and ( self.yi > 1 ) ) :     # yo = yi
                    CoulombBarrier = True

            Q, threshold = I0.determineQandThreshold( specialCases = specialCases )
            if( ( ( Q == 0 ) and ( I0.S == 0 ) and ( I0.C in S1Cs ) ) or ( CoulombBarrier ) ) :
                EMin = I0.xMin( )
                for data in datas :
                    if( data.I > 9 ) : continue
                    try :
                        EMin_ = data.xMin( )
                    except :
                        EMin_ = data.EMin( )
                    EMin = min( EMin, EMin_ )
            else :
                EMin = None
                for data in datas :
                    if( data.I > 9 ) : continue
                    if( data.I in [ 2, 5, 6, 8 ] ) :
                        print 'DATA I = ', data.I
                        continue
                    EMin_, EMinNext_ = data.getEMin_EMinNext( )
                    if( EMin is None ) : EMin, EMinNext = EMin_, EMinNext_
                    EMin, EMinNext = min( EMin, EMin_ ), min( EMinNext, EMinNext_ )

                EMin = EFloor
                if( ( threshold > 0. ) and ( EMinNext <= threshold ) ) : EMin = threshold

            thresholdCrossSectionIsZero = I0.data[0][1] == 0
            for data in datas :
                data.setQandThreshold( thresholdCrossSectionIsZero, Q, threshold, fixThresholdMode = fixThresholdMode, 
                        dThreshold_MeV = dThreshold_MeV, EMin = EMin, threshold_MeV_shiftWarning = threshold_MeV_shiftWarning )

    def getDocumentation( self ) :
        """Return None if there is not documentation file, otherwise, the text of the documentation file is returned."""

        if( self.documentation == None ) : return( None )
        return( self.documentation.getText( ) )

    def setDocumentation( self, text ) :
        """Sets self's documentation text to text."""

        if( self.documentation is None ) : self.documentation = fudgeDocumentationFile.fudgeDocumentationFile( os.path.join( self.source, 'documentation.txt' ) )
        self.documentation.setText( text )

    def heat( self, T, interpolationAccuracy = 0.002, heatAllPoints = False, doNotThin = True, EMin = 1e-11,
        heatBelowThreshold = True, removeClosePoints = True, heatAllEDomain = True ) :
        """Calls heat method for all cross section data (i.e., I = 0). and sets the temperature in the header for all data to T."""

        self.read( printWarning = False )
        datas = self.findDatas( )
        for data in datas :
            if( data.I == 0 ) : data.set( data.heat( T, interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                EMin = EMin, heatBelowThreshold = heatBelowThreshold, removeClosePoints = removeClosePoints, heatAllEDomain = heatAllEDomain ) )
            data.setTemperature( T )

    def info( self ) :
        "Prints Source, Z, mass and halflife for self."

        m = self.bdflsFile.mass( self.ZA )
        l = self.bdflsFile.halflife( self.ZA )
        print "Source =", self.source
        print "suffix =", self.suffix
        print "Working directory =", self.workDir
        print "yi = %d %s %s, Z name = %s, mass = %e, halflife = %e" % ( self.yi, self.yiTags[2], self.yiTags[3], endl_Z.endl_ZLabel( self.ZA / 1000 ), m, l )

    def ll( self, yo = None, C = None, I = None, S = None ) :
        "Prints a long listing of the data in self."

        print "Source =", self.source
        print "          Size   File name     ( yo,  C,   I,   S )    ------  yo ------  C-desc Comment"
        print "-----------------------------------------------------------------------------------"
        sC = "%%%ds" % endl_C.endl_CLabelMaxWidth( )
        sY = "%%%ds" % endl_y.endl_yLongLabelMaxWidth( )
#       sI = "%%%ds" % endl_I.endl_ILabelMaxWidth( )
        Counter = 0
        for l in self.filelist :
            if ( l.IsyoCIS( yo, C, I, S ) ) :
                s = "------"
                r = "r"
                if( l.levels == [] ) :
                    r = "-"
                else :
                    if ( len( l.levels ) == 1 ) : s = "%6d" % len( l.levels[0] )
                print "%3d %s%s %s" % ( Counter, r, l.mode, s ), 
                print l.name, "[( %2d, %2d, %3d, %3d )]" % ( l.yo, l.C, l.I, l.S ),
                s = endl_y.endl_yLongLabel( l.yo )
                if ( s == None ) : s = " "
                print sY % s,
                print sC % endl_C.endl_CLabel( l.C ),
#               print sI % endl_I.endl_ILabel( l.I ),
                print endl_I.endl_ILabel( l.I ),
                if   ( l.S == 1 ) : print "for levels",
                elif ( l.S == 2 ) : print "for pre-equilibrium ",
                elif ( l.S == 3 ) : print "for gamma production",
                elif ( l.S == 7 ) : print "for delayed",
                print
                if ( l.levels != [] ) and ( len( l.levels ) > 1 ) :
                    for q in l.levels :
                        print "    %10d                [( %2d, %2d, %3d, %3d, %e, %e, %e, %e, %e )]" % \
                            ( len( q ), l.yo, l.C, l.I, l.S, q.X1, q.X2, q.X3, q.X4, q.Q )
            Counter = Counter + 1

    def llcs( self, C = None, S = None ) :
        """Prints a long listing of cross sections, same as ll( yo = None, C, I = 0, S ).
         Deprecated method, will not be in version 4."""

        self.ll( None, C, 0, S )

    def ls( self, yo = None, C = None, I = None, S = None, w = 3 ) :
        "Prints a short listing of the data in self."

        print "Source =", self.source
        d = []
        Counter = 0
        for l in self.filelist :
            if ( l.IsyoCIS( yo, C, I, S ) ) :
                d.append( "%3d %s" % ( Counter, l.__repr__( ) ) )
            Counter = Counter + 1
        brb.tlist( d, w )

    def lscs( self, C = None, S = None, w = 3 ) :
        """Prints a short listing of cross sections, same as ls( yo = None, C = C, I = 0, S = S, w = w ).
         Deprecated method, will not be in version 4."""

        self.ls( None, C, 0, S )

    def plotcs( self, C = 1, S = None ) :
        "Plots cross section data for C_Value = C and S_Value = S.  Deprecated method, will not be in version 4."

        for l in self.filelist :
            if ( l.IsyoCIS( None, C, 0, S ) and ( l.levels != [] ) ) :
                sS = ""
                if ( l.S == 1 ) : sS = "for levels"
                l.levels.plot( ylabel = endl_C.endl_CLabel( l.C ) + " Cross section " + sS + " (barn)" )

    def reactionList( self ) :
        """Returns a list of all possible reactions (excluding C = 1). Each element of the list is a
        list of [ C, S, X1, X2, X3, X4, Q ] for a reaction."""

        csxqs = []
        CList = self.CList( )
        for C in CList :
            if( C == 1 ) : continue
            sxqs = self.SXList( C )
            for sxq in sxqs : csxqs.append( endlReactionParameters.endlReactionParameters( C, sxq[0], sxq[1], sxq[2], sxq[3], sxq[4], sxq[5] ) )
        return( csxqs )

    def read( self, yo = None, C = None, I = None, S = None, printWarning = True ) :
        "Reads in all data matching yo, C, I and S. If printWarning is True and data has already be read in, a warning message is printed."

        if ( self.DesignerIsotope ) : raise fudgeExceptions.ENDL_DesignerIsotopeReadException( "Error endlZA.read: cannot read designer isotope data from file" )
        for l in self.filelist :
            if ( l.IsyoCIS( yo, C, I, S ) ) :
                if ( l.levels != [] ) :
                    if( printWarning ) : endlmisc.printWarning( "Warning in endlZA.read: yo = %d, C = %d, I = %d and S = %d already read." % \
                        ( l.yo, l.C, l.I, l.S ) )
                else :
                    l.read( )

    def readIfNotRead( self, yo = None, C = None, I = None, S = None ) :
        "Reads in all data matching yo, C, I and S which is not currently read in."

        if ( self.DesignerIsotope ) : return
        for l in self.filelist :
            if ( l.IsyoCIS( yo, C, I, S ) ) :
                if ( l.levels == [] ) : l.read( )

    def readcs( self, C = None, S = None ) :
        "Reads in all unread cross section data with C_Value = C and S_Value = S. Deprecated method, will not be in version 4."

        self.read( yo = None, C = C, I = 0, S = S )

    def removeFile( self, yo = None, C = None, I = None, S = None ) :
        """Removes self's reference to all files matching yo, C, I and S and deletes them from the working directory.  Also see unreadFile."""

        n = len( self.filelist ) - 1
        while( n >= 0 ) :
            if ( self.filelist[n].IsyoCIS( yo, C, I, S ) ) :
                if( self.workDir != None ) :
                    cmd = "rm -rf " + self.workDir + "/" + self.filelist[n].name
                    if ( fudgeParameters.VerboseMode > 0 ) : print "executing system command >>>", cmd
                    os.system( cmd )
                del self.filelist[n]
            n -= 1

    def save( self, yo = None, C = None, I = None, S = None, path = None, file = None ) :
        """Saves all data matching yo, C, I and S to path = path and file = file.  In general file should not be specified."""

        if ( path == None ) :
            if ( self.workDir == None ) : raise Exception( "\nError in endlZA.save: project does not have a work directory." )
            path = self.workDir
        if( not os.path.exists( path ) ) : os.makedirs( path )
        for l in self.filelist :
            if ( l.IsyoCIS( yo, C, I, S ) ) :
                if ( l.levels != [] ) : 
                    file_ = file
                    if ( file_ == None ) : file_ = l.name
                    Name = os.path.join( path, file_ )
                    if ( fudgeParameters.VerboseMode > 0 ) : print "Saving data to file", Name
                    f = open( Name, "w" )
                    l.save( f )
                    f.close( )
        if( not ( self.documentation is None ) ) : self.documentation.write( os.sep.join( ( path, 'documentation.txt' ) ) )

    def savecs( self, C = None, S = None, path = None ) :
        "Saves all cross section data matching C and S. Deprecated method, will not be in version 4."

        self.save( yo = None, C = C, I = 0, S = S, path = path )

    def sumcs( self ) :
        "Returns a endl2dmath object for the sum of all cross secitions for reactions 10 <= C < 50."

        def sortfunc4sumcs_( x, y ) :
            "For internal use only."

            if ( len( x )  < len( y ) ) : return -1
            if ( len( x ) == len( y ) ) : return 0
            return 1

        cs = self.findDatas( 0, None, 0 )
        cs.sort( sortfunc4sumcs_ )
        s = 0
        for i in cs :
            if ( 10 <= i.C < 50 ) : s = s + i                           # Only want 10 <= C < 50.
        return s

    def SXList( self, C ) :
        "Returns a sorted list of S- and X-values (and Q) for all reaction of type C in self, excluding S = 7."

        sxl = []
        for f in self.filelist :
            if( f.C == C )  :
                if ( f.levels == [] ) : f.read( )
                for l in f.levels :
                    if( ( C == 15 ) and ( l.S == 7 ) ) : continue
                    sx = [ l.S, l.X1, l.X2, l.X3, l.X4, l.Q ]
                    if ( sx not in sxl ) : sxl.append( sx )
        sxl.sort( )
        return( sxl )

    def unreadFile( self, yo = None, C = None, I = None, S = None ) :
        """Removes self's reference to all files matching yo, C, I and S.  Also see removeFile."""

        for f in self.filelist :
            if ( f.IsyoCIS( yo, C, I, S ) ) : f.levels = []

    def toZAsFrame( self, workDir = None, I3ToI4lMax = 0 ) :
        """This methods boosts the data in self to invert the projectile/target particles."""

        class dummy :

            def __init__( self ) :

                self.yo = 100

        def addFile( newData, halflife ) :

            if( ( newData.yo != 0 ) and ( newData.C in [ 8, 9, 10 ] ) ) :
                yo = newData.yo + 10
                if( newData.yo > 10 ) : yo = newData.yo - 10
                newData.setYo( yo )
            file = newTarget.findFile( yo = newData.yo, C = newData.C, I = newData.I, S = newData.S )
            if( file is None ) : file = newTarget.addFile( newData.yo, newData.C, newData.I, newData.S, halflife, printWarnings = True )
            file.addEndlData( newData )

        if( self.ZA not in [ 1, 1001, 1002, 1003, 2003, 2004 ] ) : raise Exception( 'Cannot convert to %d as projectile' % self.ZA )
        self.readIfNotRead( )
        newTarget = endlZA( endl2.yoToZA( self.yi ), self.name, bdflsFile = self.bdflsFile, workDir = workDir, readOnly = False )
        newProjectileMass = self.bdflsFile.mass( self.name )
        newTargetMass = self.bdflsFile.mass( self.yi )
        halflife = self.bdflsFile.halflife( self.yi )
        datas = self.findReactionsDatas( )
        for data in datas :
            for datum in data :
                if( datum.I == 0 ) :
                    addFile( datum.toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )[0], halflife )
                    break
            data.append( dummy( ) )
            yo = -1
            for datum in data :
                if( datum.yo != yo ) :
                    if( yo >= 0 ) :
                        yoDatas.sort( )
                        Is = [ I for I, yoData in yoDatas ]
                        newDatas = []
                        if( Is == [ 1 ] ) :
                            newDatas = yoDatas[0][1].toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )
                        elif( Is == [ 4 ] ) :
                            newDatas = yoDatas[0][1].toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )
                        elif( Is == [ 1, 4 ] ) :                # The logic in this elif assumes toZAsFrame returns one data set each.
                            print 'WARNING from endlZA.toZAsFrame: Discarded I = 1 data from I = [ 1, 4 ] set for oldYi = %d, oldYo = %d and C = %d' % \
                                ( self.yi, yoDatas[0][1].yo, yoDatas[0][1].C )
                            newDatas_ = yoDatas[0][1].toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )
                            newDatas = newDatas_ + yoDatas[1][1].toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )
                        elif( Is == [ 1, 3 ] ) :
                            print 'WARNING from endlZA.toZAsFrame: Converting I = 1,3 to I = 4, lMax = %d for oldYi = %d, oldYo = %d and C = %d' % \
                                ( I3ToI4lMax, self.yi, yoDatas[0][1].yo, yoDatas[0][1].C )
                            I4 = yoDatas[1][1].convertToI4( yoDatas[0][1], I3ToI4lMax )
                            newDatas = I4.toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )
                        else :
                            raise Exception( 'Distribution combination I = %s not supported' % ( Is ) )
                        for newData in newDatas : addFile( newData, halflife )
                    yoDatas = []
                    if( datum.yo != 0 ) : yo = datum.yo
                if( yo == 100 ) : break
                if( datum.I in [ 1, 3, 4 ] ) : yoDatas.append( [ datum.I, datum ] )
                if( datum.I == 9 ) : addFile( datum.toZAsFrame( newProjectileMass, newTargetMass, halflife, self.bdflsFile )[0], halflife )
        return( newTarget )
        
    def toGND( self, evaluationLibrary, evaluationVersion, xenslIsotopes = None, testing = False, verbose = 0 ) :

        import fudge.legacy.converting.endlToGND as endlToGND
        return( endlToGND.toGND( self, evaluationLibrary, evaluationVersion, xenslIsotopes = xenslIsotopes, verbose = verbose ) )

    def adjustcs( self, C = 10, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None, f = 5.e-5, csRelChangeNotice = .1 ) :
        """Adjust the cross section given by C, S, X1, X2, X3, X4 and Q  so that the total cross section remains fixed.  
        If C = 1 (i.e., total cross section), then the total is adjusted to keep the others at there current values."""

        s = self.sumcs( )                                           # Calulate the new sum.
        t = self.cs( 1 )                                            # Get the current total cross section.
        if ( C == 1 ) :                                             # Adjust the total cross section.
            old_a = t.copyData( )
            t.set( s )
            a = t
        else :                                                      # Adjust some other cross section.
            d = t - s                                               # The difference.
            s = ( t + s ) / 2.
            for i in range( len( d ) ) :                            # Set small relative differences to zero.
                if ( abs( d.data[i][1] ) < f * s.data[i][1] ) : d.data[i][1] = 0.
            d.trim( )                                               # Remove trailing and ending zeros.
            if ( len( d ) == 0 ) : return
            a = self.cs( C, S, X1, X2, X3, X4, Q )                      # Get cross section to adjust.
            old_a = a.copyData( )
            a.set( a + d )                                      # Adjust cross section.
            c = 0
            emin = 0.
            for e in a.data :                                   # Check for negative cross sections.
                if ( e[1] < 0. ) :
                    c = c + 1
                    emin = min( emin, e[1] )
                    if ( c < 5 ) : endlmisc.printWarning( "Warning in endlZA.adjustcs: negative cross section %e %e set to 0." % (e[0], e[1]) )
                    e[1] = 0.
            if ( c > 0 ) : endlmisc.printWarning( "Warning in endlZA.adjustcs: %d negative cross sections set to 0 (%e)." % ( c, emin ) )
        s = ( old_a + a ) / 2.
        c = 0
        i = 0
        j = 0
        rcmax = 0.
        for e in s.data :                                           # Check for large, compared to csRelChangeNotice,
            if ( e[1] != 0. ) :                                     # relative difference between a and old_a.
                v, i = a.GetValueNIndex( e[0], i )
                vold, j = old_a.GetValueNIndex( e[0], j )
                rc = ( v - vold ) / e[1]
                if ( abs( rc ) > csRelChangeNotice ) :
                    if ( abs( rc ) > abs( rcmax ) ) : rcmax = rc
                    c = c + 1
                    if ( c < 5 ) :
                        endlmisc.printWarning( "Warning in endlZA.adjustcs: relative change of %e." % rc )
                        endlmisc.printWarning( "       Value at E = %e changed from %e to %e" % ( e[0], vold, v ) )
        if ( c > 0 ) : endlmisc.printWarning( "Warning in endlZA.adjustcs: %d relative changes greater than %e found (%e)." \
            % ( c, csRelChangeNotice, rcmax ) )

    def adjustcsWithWeights( self, CsNWeights, f = 5.e-5, total = None ) :
        """Adjusts the cross sections given by CsNWeights, so that the total cross section 
        remains fixed. CsNWeights is a list of C-value and weight. For example, a CsNWeights 
        of [ [ 10, 6 ], [ 11, 4 ], [ 12, 2 ] ] would adjust the C = 10 cross section by
        50% of the different, the C = 11 cross section by 33.3% of the different and the C = 12
        cross section by 16.7% of the different. Relative total cross section less than f are set to 0."""

        s = self.sumcs( )                                       # Calulate the new sum.
        if( total == None ) : total = self.cs( 1 )              # Get the current total cross section.
        d = total - s                                           # The difference.
        s = ( total + s ) / 2.
        i = 0
        for xy in d.data :                                      # Set small relative differences to zero.
            if ( abs( xy[1] ) < f * s.data[i][1] ) : xy[1] = 0.
            i += 1
        d.trim( )                                               # Remove trailing and ending zeros.
        if ( len( d ) == 0 ) : return
        csNWs = []
        totalWeight = 0
        x = 1e-8
        for CNW in CsNWeights :
            csList = self.findDatas( C = CNW[0], I = 0 )
            for cs in csList :
                w = CNW[1]
                if( CNW[1] == None ) : w = cs
                totalWeight = totalWeight + w
                csNWs.append( cs )
        d = d / totalWeight
        for cs in csNWs :
            cs.set( cs * ( 1. + d ) )

    def syncLevelEnergies( self, checkOnly = False, isomerMinimumHalflife = None ) :
        """Synchronizes the level energies of the target, and discrete states of the residual
        with a database containing the accepted level energies for all nuclei. The residual
        level energies are synchronized for both the inelastic S=1 data (X1), and isomers (X4).
        Checks that the isomer level labels m? are consistent with the levels given in the database.
        If checkOnly = True then only check against database and print differences,
        else if checkOnly = False then make the changes to the data.
        If isomerMinimumHalflife is given then sets this variable in the Levels class before performing sync."""
        
        if isomerMinimumHalflife : xensl.isomerMinimumHalflife = isomerMinimumHalflife
        targetLevels = xensl.Levels( self.ZA )
        if self.suffix :
            targetEnergy = targetLevels( self.suffix ).energy
            # check isomer level is actually labeled as a metastable state
            if 'm' in self.suffix :
                for m in targetLevels.metastables.keys( ) :
                    if targetEnergy == targetLevels( m ).energy : break
                else:
                    endlmisc.printWarning( 'WARNING: Isomer level energy %s, not a known isomer for target ZA=%s' % ( targetEnergy, self.sZA ) )
        else :
            targetEnergy = 0.

        for d in self.findDatas( ) :
            # check target energy levels
            if targetLevels( d.ELevel ).energy != targetEnergy :
                endlmisc.printWarning( 'WARNING: Isomer label %s with level energy %s for ZA=%s, does not match level energy %s given in database' % ( self.suffix, d.ELevel, self.ZA, targetEnergy ) )
            if d.ELevel != targetEnergy :
                endlmisc.printWarning( 'Changing level energy from %s to %s for target ZA=%s\n%s' % ( d.ELevel, targetEnergy, self.sZA, d ) )
                if not checkOnly : d.ELevel = targetEnergy
            residualZA = endl2.residualZA_yos_Q( d.yi, d.ZA, d.C )[0]
            if not residualZA : continue
            residualLevels = xensl.Levels( residualZA )
            # check S=1 levels
            if d.S == 1 and d.X1 != residualLevels( d.X1 ).energy :
                endlmisc.printWarning( 'Changing discrete state energy from %s to %s for residual ZA=%s\n%s' % ( d.X1, residualLevels( d.X1 ).energy, residualZA, d ) )
                if not checkOnly : d.X1 = residualLevels( d.X1 ).energy
            # check residual isomer levels
            if d.X4 != residualLevels( d.X4 ).energy :
                endlmisc.printWarning( 'Changing isomer level energy from %s to %s for residual ZA=%s\n%s' % ( d.X4, residualLevels( d.X4 ).energy, residualZA, d ) )
                if not checkOnly : d.X4 = residualLevels( d.X4 ).energy
            # check isomer level is actually labeled as a metastable state
            if d.X4 :
                for m in residualLevels.metastables.keys( ) :
                    if residualLevels( d.X4 ).energy == residualLevels( m ).energy : break
                else:
                    endlmisc.printWarning( 'WARNING: Isomer level energy %s, not a known isomer for residual ZA=%s\n%s' % ( d.X4, residualZA, d ) )
