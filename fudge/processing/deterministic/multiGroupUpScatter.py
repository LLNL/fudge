# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, sys, subprocess

import fudge as fudgeModule
from LUPY import subprocessing as subprocessingModule
from LUPY import times as timesModule
from fudge.processing.deterministic import transferMatrices as transferMatricesModule

def SnElasticUpScatter(style, tempInfo, comment=None):
    """
    Generate input and call processing code to generate a transfer matrix for two-body angular distribution.
    If the distribution is actually made up of two different forms in different energy regions, this function
    calls itself in the two regions and sums the result.
    """

    logFile = tempInfo['logFile']
    workDir = 'upscatter.work'
    
    if not os.path.exists(workDir):
        os.mkdir(workDir)
    projectileGroupBoundaries = tempInfo['groupBoundaries']
    s = '\n'.join(transferMatricesModule.GBToString('Projectile', projectileGroupBoundaries, None).split()[8:])  ### exclude some header stuff we dont want to parse in C
    dataFile = open('%s/groupStructure.dat'%workDir, 'w')
    dataFile.write(s)
    dataFile.close()
    
    ### run the calcUpscatter code -- ./calcUpscatterKernel 42092 3.1e-8 ==> upscatterEMuEp.out
    cmd = os.path.join(os.path.split(os.path.dirname(fudgeModule.__file__))[0], 'bin', 'calcUpscatterKernel')
    cmd += '  %s %g %g %g %g %d ' % ('groupStructure.dat', tempInfo['targetMassRatio'], tempInfo['temperature'],
            tempInfo['minEval'], tempInfo['maxEval'], tempInfo['legendreOrder'])
    logFile = executeCommand(cmd, workDir, 'legendre')
    
    ### call the transferMatrix parseOutput on the outfile
    TM1, TME, paras = transferMatricesModule.parseOutputFile('%s/upscatterLegendre.out'%workDir, 0, firstLine = "upscatter: version 1 \n")
    maxG = paras['maxIncidentEnergyGroup']
    
    if TM1 is not None:
        negative_l0Counter = checkNegative_l0(TM1, "0", logFile)
    if TME is not None:
        negative_l0Counter += checkNegative_l0(TME, "E", logFile)
    if negative_l0Counter > 0:
        print('WARNING: %d negative l=0 elements found in transfer matrix' % negative_l0Counter)
    APE = parseAveEnergyOutputFile('%s/AveEnergy.out'%workDir)

    return TM1, TME, APE, paras['maxIncidentEnergyGroup']
        
def executeCommand(cmd, workDir, workFile):

    transferMatrixExecute = cmd.split()[0]
    if( os.path.exists( transferMatrixExecute ) ) :
        srcPath = os.path.abspath( './' )
    else :
        transferMatrixExecute = os.path.split(transferMatrixExecute)[1]
        srcPath = os.path.join(sys.prefix, 'bin', transferMatrixExecute)
        if( os.path.exists( srcPath ) ) :
            transferMatrixExecute = srcPath
        else :
            srcPath = os.path.join( os.path.dirname( os.path.dirname( fudgeModule.__file__ ) ), 'bin' )
            transferMatrixExecute = os.path.join( srcPath, transferMatrixExecute )

    cmd = [transferMatrixExecute] + cmd.split()[1:]

    if not os.path.exists(workDir):
        os.mkdir(workDir)
    current_directory = os.getcwd()
    os.chdir(workDir)
    fullFileName = 'upscatter_%s' % workFile 
    infoFile = '%s.log' % fullFileName
    
    t0 = timesModule.Times( )
    try :
        fCmd = open('%s.sh' % workFile, 'w')
        fCmd.write(' '.join(cmd) + '\n')
        fCmd.close()
        status, stdout, stderr = subprocessingModule.executeCommand( cmd, stdout = infoFile , stderr = subprocess.STDOUT )
    except :
        fErr = open( fullFileName + ".err", "w" )
        fErr.close( )
        raise
        
    fOut = open( infoFile, 'a' )
    fOut.write( str( t0 ) + "\n" )
    if( status != 0 ) :
        print("status = ", status)
        raise Exception( 'Upscatter failed for %s %s' % ( cmd, fullFileName ) )
    fOut.close()
    
    os.chdir(current_directory)
    return os.path.join(workDir, infoFile)

### checker for resulting legendre matrices
def checkNegative_l0( TM_EEpL, weight, infoFile ) :

    fOut = open( infoFile, 'a' )
    negative_l0Counter, largestNegative = 0, 0.
    for ig in sorted( TM_EEpL.keys( ) ) :
        TM_EpL = TM_EEpL[ig]
        for ih in sorted( TM_EpL.keys( ) ) :
            l0Cell = TM_EpL[ih][0]
            if( l0Cell < 0 ) :
                negative_l0Counter += 1
                largestNegative = min( largestNegative, l0Cell )
                fOut.write( 'negative l=0 for weight %s at row %3d column %3d: %s\n' % ( weight, ig, ih, l0Cell ) )
    fOut.close()
    return( negative_l0Counter )


def parseAveEnergyOutputFile(fileName) :
    try :
        f = open( fileName )
    except :
        raise Exception( 'Could not open file = %s' % fileName )
    ls = f.readlines( )
    f.close( )
    
    APE = []
    for s in ls:
        APE.append([float(x) for x in s.split()[-2:]])
    
    return APE    
        


def rescaleCrossSection( groupedCrossSec, UpScatterMatrix ) :

    for g,sigma in enumerate(groupedCrossSec):
        summedarray = [ UpScatterMatrix[g][h][0] for h in range(len(UpScatterMatrix[g])) ]
        sumXS = sum(summedarray) 
        if sumXS != 0.0:
            ratio = sigma/sumXS
            for l in range(len(UpScatterMatrix[g][0])):
                for h in range(len(UpScatterMatrix[g])) :
                    if (UpScatterMatrix[g][h][l] != 0.0):
                        UpScatterMatrix[g][h][l] = UpScatterMatrix[g][h][l]*ratio

    return UpScatterMatrix

    
def analyzeMatrix(TM, tempInfo):
    #import numpy as np
        
    def plotMatrix(TM, tempInfo, lastElement):
        
        
        #s = '\n'.join(transferMatricesModule.GBToString( 'Projectile', gb, None ).split()[8:])  ### exclude some header stuff we dont want to parse in C
    
        import numpy as np
        from fudge.vis.matplotlib import plot_matrix as plotModule
        #xR,yR,lR = len(TM),len(TM[0]),len(TM[0][0])
        xR,yR,lR = lastElement+5,lastElement+5,len(TM[0][0])
        gb = tempInfo['groupBoundaries'].boundaries.values.values[0:xR]
        TM_A = np.zeros((xR,yR))
        for g in range(xR):
            for h in range(yR):
                TM_A[g,h] = TM[g][h][0]
        plotModule.plot_matrix(TM_A,zlog=True,xyTitle=('Outgoing Energy Group','Incident Energy Group')) #,energyBoundariesX=gb)  #,xylog=True,energyBoundariesX=,energyBoundariesY=
        plotModule.pyplot.show()
        tempStr = str(tempInfo['temperature']).replace(' ','').replace('/','over')
        plotModule.pyplot.savefig('upscatterTM_L0_%s.png'%tempStr)

    xR,yR,lR = len(TM),len(TM[0]),len(TM[0][0])
    #print(xR,yR,lR)
    
    elements = []        
    lastElement = 0        
    for g in range(xR-1,0,-1):
        if TM[g][g][0] > 0.0 : 
            lastElement = g
            break
                
    #print('last nonzero matrix element = ', lastElement)
    plotMatrix(TM, tempInfo, lastElement)
            
            
    return lastElement
