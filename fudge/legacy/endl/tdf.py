#! /usr/bin/env python
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



""" This module contains the methods to run tdfgen, producing tdf files for light ion thermonuclear reactions. Methods include
prepping the endlZA input and running the code (processTDF_Reaction), and various wrappers to provide running options.


"""

import os,sys,getopt,math
from fudge.legacy.endl import endlProject
from fudge import fudgeDefaults
#import endlProject,endl2,bdfls,endl_Z,endl2dmathClasses,fudgeFileMisc,endlmisc                          ### for live version of fudge at /usr/apps/fudge/current/Src/
from fudge.legacy.endl import  endl2, bdfls,endl_Z,endl2dmathClasses,fudgeFileMisc,endlmisc  ### for private copies of fudge
from endl2dmathClasses import endl2dmath



### find path to tdfgen
#fullpath = os.path.realpath(__file__)
#print fullpath
#thisdir,filename = os.path.split(fullpath)
#print thisdir,filename
#containerdir,subdir = os.path.split(thisdir)
#print containerdir,subdir
#if os.path.isfile(os.path.join(containerdir,subdir,'tdfgen')):
#    pathToTdfgen = os.path.join(containerdir,subdir,'tdfgen')
#elif os.path.isfile(os.path.join(containerdir,'bin','tdfgen')):
#    pathToTdfgen = os.path.join(containerdir,'bin','tdfgen')
#else:
#    pass ## throw exception that we can't find tdfgen
#    raise OSError( "cant find tdfgen in %s" % ( containerdir ))
#print  pathToTdfgen   


#yi, ZA, C = 4, 1003, 12 # H3__H3_n1__n1__He4
#yi, ZA, C = 3, 1002, 11 # H2__H2_n1__He3
#yi, ZA, C = 3, 1002, 40 # H2__H2_H1__H3
#yi, ZA, C = 5, 1002, 40 # He3__H2_H1__He4
#yi, ZA, C = 3, 2003, 40 # H2__He3_H1__He4
#yi, ZA, C = 4, 1002, 11 # H3__H2_n1__He4
#yi, ZA, C = 3, 1003, 11 # H2__H3_n1__He4
#yi, ZA, C = 1, 3006, 42 # n1__Li6_H3__He4
#yi, ZA, C = 2, 3006, 44 # H1__Li6_He3__He4
#yi, ZA, C = 2, 3007, 45 # H1__Li7_He4__He4
#yi, ZA, C = 3, 3006, 11 # H2__Li6_n1__Be7
#yi, ZA, C = 3, 3006, 40 # H2__Li6_H1__Li7
#yi, ZA, C = 3, 3006, 45 # H2__Li6_He4__He4

#yi, ZA, C = 2, 5011, 37 # H1__B11_He4__He4__He4
#yi, ZA, C = 3, 1002, 46 # H2__H2_g__He4
#yi, ZA, C = 3, 3006, 39 # Li6__H2_H1__H3__He4
#yi, ZA, C = 1, 4009, 46 # n1__Be9_g__Be10
#yi, ZA, C = 5, 4009, 46 # He3__Be9_g__C12
#yi, ZA, C = 1, 2003, 46 # n1__He3_g__He4
#yi, ZA, C = 1, 1001, 10 # n1__H1_n1__H1 (elastic)


default_C_list = ( 11, 40, 41, 42, 12, 44, 45 )
default_C_list = ( 11, 40, 41, 42, 12, 44, 45 , 39) # try to get Li6__H2_H1__H3__He4 too!

default_lib_path = "/usr/gapps/data/nuclear/endl_official/"

default_ZAs = [ 1, 1001, 1002, 1003, 2003, 2004, 3006, 3007 ]

tdfgenpath = os.path.join(fudgeDefaults.DefaultFudgePath,'bin','tdfgen')

def run_command( command, verbose = False ):
    if verbose: print command
    try:
        os.system( command )
    except Exception as err:
        os.sys.exit( 'command failed: %s ' % command )
                
def processTDF_Reaction( target, C, S = None, X1 = None, X2 = None, X3 = None, X4 = None, Q = None, outputDir = None, 
        workDir = None, movePoints = None, movePointsEps = 1e-4, addPoints = None, libraryName = "ENDL2008", libraryVersion = "2", verbose = False, overWrite=True, dryrun=False ) :

    endl_bdfls = bdfls.getDefaultBdfls( )
    AMUToMeV = endl_bdfls.constant( 4 )
    dataList = target.findDatas( C = C, I = 0, S = S, X1 = X1, X2 = X2, X3 = X3, X4 = X4, Q = Q )
    for xsec in dataList:
        xsec.set( xsec.thin( interpolationAccuracy = 1e-4 ) )
        if( movePoints != None ) :
            for x, y in xsec.data :
                for x_, y_ in movePoints :
                    if( abs( x - x_ ) < movePointsEps * x ) :
                        xsec.setValue( x, y_ )
                        break
        if( addPoints != None ) :
            for x, y in addPoints : xsec.setValue( x, y )

## thicken the data if needed ##
        x1 = None
        data = []
        for x2, y2 in xsec.data :
            if( not( x1 is None ) ) :
                if( ( x2 / x1 ) > 1.2 and y2 > 0.0 and y1 > 0.0) :
                    x = x1
                    a = math.log( y2 / y1 ) / math.log( x2 / x1 )
                    while( True ) :
                        x = x * 1.1
                        if( x > 0.999 * x2 ) : break
                        y = y1 * math.pow( x / x1, a )
                        data.append( [ x, y ] )
            data.append( [ x2, y2 ] )
            x1, y1 = x2, y2
        xsec.set( endl2dmathClasses.endl2dmath( data ) )

        yi = endl2.ZAToYo( target.yi )
        yiZA = endl2.yoToZA( yi )
        yiZ, yiA = endl2.ZandAFromZA( yiZA )

        electronMass = endl_bdfls.mass( 9 )
        residualZA, yos, Q = endl2.residualZA_yos_Q( target.yi, target.ZA, C )
        projectileMass = endl_bdfls.mass( target.yi ) - yiZ * electronMass
        targetMass = endl_bdfls.mass( target.ZA ) - target.Z * electronMass
        projectileZ = yiZ
        targetZ = target.Z
        projectileA = yiA
        targetA = target.A

        ZA = target.ZA
        ZAZ, ZAA = endl2.ZandAFromZA( ZA )

        if( projectileMass > targetMass ) :
            reaction = '%s%d+%s%d' % ( endl_Z.endl_ZSymbol( yiZ ), yiA, endl_Z.endl_ZSymbol( ZAZ ), ZAA )
            #print "Skipping reaction "+reaction+", it is in inverse kinematics"
            return
        else :
            reaction = '%s%d+%s%d=' % ( endl_Z.endl_ZSymbol( ZAZ ), ZAA, endl_Z.endl_ZSymbol( yiZ ), yiA )

        outGoing = []
        for yo in yos : 
            if( yo < 1000 ) :
                outGoing.append( [ endl2.yoToZA( yo ) ] )
            else :
                outGoing.append( [ yo ] )
        outGoing.append( [ residualZA ] )

        for i in outGoing :
            iZA = i[0]
            Z, A = endl2.ZandAFromZA( iZA )
            thisMass = endl_bdfls.mass( iZA ) - Z * electronMass
            if residualZA == iZA: thisMass += xsec.X1
            i.insert( 0, thisMass )
            i.append( Z )
            i.append( A )
        outGoing.sort( )
        s = ''
        for mass, iZA, Z, A in outGoing :
            reaction += '%s%s%d' % ( s, endl_Z.endl_ZSymbol( Z ), A )
            s = '+'
        
        if xsec.S != 0:
            level_index = dataList.index(xsec)
            if level_index > 0:
                reaction += '_e%i' % ( level_index )
        
        if dryrun:
            #print target.yi, target.ZA, C, reaction
            return target.yi, target.ZA, C, reaction

        outputStr  = [ '## Fudge generated data for tdfgen version:1.0.0' ]
        outputStr.append( '## Data generated from:fudge' )
        outputStr.append( '' )
        outputStr.append( '## Reaction:%s' % reaction )
        outputStr.append( '' )
        outputStr.append( '# Masses of particles in MeV.' )
        outputStr.append( '## Mass of projectile:%.12e' % ( projectileMass * AMUToMeV ) )
        outputStr.append( '## Mass of target:%.12e' % ( targetMass * AMUToMeV ) )
        outputStr.append( '## Z of projectile:%.12e' % ( projectileZ ) )
        outputStr.append( '## A of projectile:%.12e' % ( projectileA ) )
        outputStr.append( '## Z of target:%.12e' % ( targetZ ) )
        outputStr.append( '## A of target:%.12e' % ( targetA ) )

        outputStr.append( '' )
        outputStr.append( '## Number of final particles:%d' % len( outGoing ) )
        for mass, iZA, Z, A in outGoing : outputStr.append( '## %.12e' % ( mass * AMUToMeV ) )

        outputStr.append( '' )
        outputStr.append( '## Lab or CM frame:Lab' )


        xsec_ = []
        x1 = None
        f = 1.03
        for x2, y2 in xsec.data :
            if( x1 != None ) :
                x = f * x1
                while( f * x <= x2 ) :
                    xsec_.append( [ x, xsec.getValue( x ) ] )
                    x = f * x
            xsec_.append( [ x2,  y2 ] )
            x1, y1 = x2, y2
        xsec_ = endl2dmathClasses.endl2dmath( xsec_ )
        outputStr.append( '## Number of data points:%d' % len( xsec_ ) )
        outputStr.append( '# E(MeV)    Sigma( barn )' )
        outputStr.append( '#------------------------' )
        outputStr.append( xsec_.toString( ) )

        outputStr = '\n'.join( outputStr )

        inputFile = fudgeFileMisc.fudgeTempFile( dir = workDir )
        inputFile.write( outputStr )

        inputName = inputFile.getName( )
        if verbose: print inputName

        #
        if not overWrite:
            if os.path.exists( reaction + ".tdf" ): raise IOError( reaction + ".tdf exists!" )

        # Assemble list of commands to run during processing
        commands = [ "cp " + inputName + " " + reaction + ".tdfgen" ]
        print "starting tdfgen"
        if libraryVersion in [ None, '' ]:
            commands.append( '%s -i %s -o %s -d %s' % ( tdfgenpath, reaction + ".tdfgen", reaction + ".tdf", libraryName ) )
        else:
            commands.append( '%s -i %s -o %s -d %s -n %s' % ( tdfgenpath, reaction + ".tdfgen", reaction + ".tdf", libraryName, libraryVersion ) )

        # Now process
        if verbose: print
        for command in commands:  run_command( command, verbose )
        if verbose: print

        # Cleanup
        if not os.path.exists( reaction + ".tdf" ): raise RuntimeError( "tdfgen did not produce the output file %s !" % ( reaction + ".tdf" ) )
        if outputDir != None:
            if not os.path.exists( outputDir ): raise OSError( "output dir %s does not exist" % ( outputDir ) )
            run_command( 'mv %s %s' % ( reaction + ".tdf", outputDir ), verbose )
        import glob
        if glob.glob( '*core*' ) != []: raise RuntimeError( "tdfgen coredumped! better check it out!" )
        run_command( 'rm %s ' % ( reaction + ".tdfgen" ), verbose )
        
        #print target.yi, target.ZA, C, reaction
        return target.yi, target.ZA, C, reaction

# ------------- process_one_reaction -----------------
def process_one_reaction( yi, ZA, C, libraryName = 'endl2009.0', libraryVersion = '', libPath = default_lib_path, outputDir = None, workDir = 'work', verbose = False, dryrun = False  ):
    '''
       process one reaction, use this for testing
    
       projectile designator: yi = ( 1,n; 2,p; 3,d; 4,t; 5,3He; 6,He; 7,g )
       target ZA = 1000 * Z + A                                            
       I = 0 is cross section data.                                        
    '''
    
    if verbose:
        print '---- process_one_reaction ----'
        print 'yi:',yi
        print '    za:',ZA
        print "        C:",C
    e = endlProject( os.path.join(libPath,libraryName,libraryVersion), projectile = yi ) 
    target = e.readZA( ZA )
    if verbose: target.ll()         
    target.read( I = 0 )
    reactionname = processTDF_Reaction( target, C, workDir = workDir, outputDir = outputDir, libraryName = libraryName, libraryVersion = libraryVersion, verbose = verbose, dryrun = dryrun )
    if dryrun:
        return reactionname
            

# ------------- process_one_evaluation -----------------
def process_one_evaluation( yi, ZA, libraryName = 'endl2008.2', libraryVersion = '', libPath = default_lib_path, outputDir = None, workDir = 'work', verbose = False, dryrun = False ):
    '''process all reactions for an evaluation that we possibley can'''
    e = endlProject( database = os.path.join(libPath,libraryName,libraryVersion), projectile = yi )

    try:	target = e.readZA( ZA )
    except: return 

    if verbose:
        print '---- process_one_evaluation ----'
        print 'yi:',yi
        print '    za:',ZA
    
    if verbose: target.ll()         
    target.read( I = 0 )
    CList = filter( lambda x: x in default_C_list, target.CList() )
    
    cutOffEnergy = 0.010 # 10 keV cutoff, if threshold > this value, then don't make a TDF file

    #if CList == []: print "No C's to process for yi =", yi, "za =", ZA,"!"
    reactionnames = []
    for C in CList:
        if target.findDatas( I=0, C=C )[0].data[0][0] > cutOffEnergy:
            #print "        C:", C, "first energy point too high"
            continue
        if verbose: print "        C:",C
        try:
            reactionname = processTDF_Reaction( target, C, workDir = workDir, outputDir = outputDir, libraryName = libraryName, libraryVersion = libraryVersion, verbose = verbose, dryrun = dryrun )
            reactionnames.append(reactionname)
        except RuntimeError, err:
            print err
            print "\n.... continuing\n"
    if dryrun:
        #print reactionnames
        return reactionnames        

# ------------- process_sublibrary -----------------
def process_sublibrary( yi, libraryName = 'endl2008.2', libraryVersion = '', libPath = default_lib_path, outputDir = None, workDir = 'work', verbose = False, dryrun = False ):
    for ZA in [ 1, 1001, 1002, 1003, 2003, 2004, 3006, 3007 ]:
        if verbose: 
            print '---- process_sublibrary ----'
            print '    za:',ZA
        process_one_evaluation( yi, ZA, libraryName = libraryName, libraryVersion = libraryVersion, libPath = libPath, outputDir = outputDir, workDir = workDir, verbose = verbose, dryrun = dryrun )


# ------------- process_everything -----------------
def process_everything( libraryName = 'endl2008.2', libraryVersion = '', libPath = default_lib_path, outputDir = None, workDir = 'work', verbose = False, dryrun = False ):
    for yi in range( 1, 7 ):
        if verbose: 
            print '---- process_everything ----'
            print 'yi:',yi
        process_sublibrary( yi, libraryName = libraryName, libraryVersion = libraryVersion, libPath = libPath, outputDir = outputDir, workDir = workDir, verbose = verbose, dryrun = dryrun )

def process_names(libraryName = 'endl2008.2', libraryVersion = '', libPath = default_lib_path, outputDir = None, workDir = 'work', verbose = False, dryrun = True, ZAlist=default_ZAs  ):
    reactionnames = []
    for yi in range( 1, 7 ):
        for ZA in ZAlist:
            A = process_one_evaluation( yi, ZA, libraryName = libraryName, libraryVersion = libraryVersion, libPath = libPath, outputDir = outputDir, workDir = workDir, verbose = verbose, dryrun = dryrun )
            reactionnames.extend(A)
    #print reactionnames  
    return reactionnames        

# ------------- usage -----------------
def usage():
    '''print usage message'''
    print "Usage: tdf.py (options) C1 C2 C3 ..."
    print 
    print '    here the arguments C1, C2, C3 ... is the list of C numbers to process when specifying a yi and za'
    print
    print '    valid options:'
    print '         -h, --help      prints this message'
    print '         -v              prints tdf version number'
    print '         -V, --verbose   set vebose flag'
    print "         -o, --output    sets the output directory (requires name, otherwise defaults to '.')"
    print "         --yi            sets the projectile ( 1,n; 2,p; 3,d; 4,t; 5,3He; 6,He; 7,g ) (defaults to _all of them_)"
    print "         --za            sets the ZA of the target (enter as Z*100+A, defaults to [ 1, 1001, 1002, 1003, 2003, 2004, 3006, 3007 ] )"
    print "         --library       sets the library (default is endl2008.2)"
    print "         --libversion    sets the library version (default is '', is appended onto the library string)"
    print "         --libpath       sets the libpath"
    print "         --workdir       sets fudge's working directory"
    print "         --all           process all reactions in all evaluations for all projectiles in library"
    print "         --dryrun        runs through data lookup and name munging, but does write out any files"

# ------------- version -----------------
def version():
    '''prints version information, uses tdfgen to do it'''
    os.system( 'tdfgen -v' )


# -------------- main ----------------
def main( arglist ):

    try:
        opts, args = getopt.getopt( arglist, "ho:vV", [ "help", "all", "verbose", "output=", "yi=", "za=", "work=", "library=", "libversion=", "libpath=", "dryrun", "tdfgenpath=" ] )

    except getopt.GetoptError, err:
        # print help information and exit:
        print str( err ) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    output = None
    verbose = False
    yi = None
    za = None
    library = 'endl2009.0'
    libversion = ''
    libpath = default_lib_path
    workdir = None
    all = False
    dryrun = False
    tdfgenpath = '/usr/apps/fudge/current/bin/'
        
    for o, a in opts:
        if o in ( "-V", "verbose" ):
            verbose = True
        elif o in ( "-v", "version" ):
            version()
            sys.exit()
        elif o in ( "-h", "--help" ):
            usage()
            sys.exit()
        elif o in ( "-o", "--output" ):
            output = a
        elif o == "--yi":
            yi = int( a )
        elif o == "--za":
            za = int( a )
        elif o == "--library":
            library = a
        elif o == "--libversion":
            libversion = a 
        elif o == "--libpath":
            libpath = a
        elif o == "--work":
            workdir = a
        elif o == "--tdfgenpath":
            tdfgenpath = a
        elif o == "--all":
            all = True
        elif o == "--dryrun":
            dryrun = True
        else:
            assert False, "unhandled option"

    if verbose:
        print '---- options ----'
        print "output:",output
        print "yi:",yi
        print "za:",za
        print "library:",library
        print "libversion:",libversion
        print "libpath:",libpath
        print "workdir:",workdir
        print "tdfgenpath:",tdfgenpath

    if all:
        process_everything( libraryName = library, libraryVersion = libversion, libPath = libpath, outputDir = output, workDir = workdir, verbose = verbose, dryrun = dryrun )
    elif za == None and yi != None:
        process_sublibrary( yi, libraryName = library, libraryVersion = libversion, libPath = libpath, outputDir = output, workDir = workdir, verbose = verbose, dryrun = dryrun )
    elif za != None and yi != None:
        if len( args ) == 0:    
            process_one_evaluation( yi, za, libraryName = library, libraryVersion = libversion, libPath = libpath, outputDir = output, workDir = workdir, verbose = verbose, dryrun = dryrun )
        else:
            for C in args:
                process_one_reaction( yi, za, int( C ), libraryName = library, libraryVersion = libversion, libPath = libpath, outputDir = output, workDir = workdir, verbose = verbose, dryrun = dryrun )
                
    else:
        print "must specify a yi with the --yi option or use the --all option"
        sys.exit()
        

# ------------------------------
if __name__ == "__main__":
    main( sys.argv[1:] )
    
