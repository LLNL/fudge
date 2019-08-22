import os

def readXYData(filename, commentCharacter='#'):
    '''
    Read in a file containing 2 columns of x, y, dy
    Lines beginning with commentCharacter are ignored
    '''
    if not os.path.exists( filename ): raise IOError( "XY data file %s not found"%filename)
    results = []
    for line in open(filename).readlines():
        try:
            if line.startswith(commentCharacter): continue
            results.append( map( float, line.split()[0:2] ) )
        except ValueError: pass
    return results 


def readXYdYData(filename, commentCharacter='#'):
    '''
    Read in a file containing 3 columns of x, y, dy
    Lines beginning with commentCharacter are ignored
    '''
    if not os.path.exists( filename ): raise IOError( "XYdY data file %s not found"%filename)
    results = []
    for line in open(filename).readlines():
        try:
            if line.startswith(commentCharacter): continue
            results.append( map( float, line.split()[0:3] ) )
        except ValueError: pass
    return results 


def readEvaluation( filename, targ=None, proj=None, verbose=True, skipBadData=True, reconstructResonances=True, continuumSpectraFix=False,
                    skipCovariances=False, verboseWarnings=True, printBadNK14=False, ignoreBadDate=True, acceptBadMF10FissionZAP=True):
    '''
    Read in an evaluation in either Fudge/GNDS, ENDF or AMPX/BOF format and return result as Fudge classes
    '''
    import fnmatch
    
    # OK, try AMPX first
    if fnmatch.fnmatch(filename,'*.ampx*'): 
        try:
            import ampx2fudge, ampx
            ampx_za, bounds = ampx.readEvaluation( filename, str(targ), str(proj) )
            return [ 
                ampx2fudge.convertAmpxNuclideToFudgeReactionSuite(ampx_za, bounds), 
                fudge.gnds.covariances.covarianceSuite.covarianceSuite() ]
        except ImportError:
            print "WARNING: Could not load AMPX module.  Is it in your path?"
    
    else:
        firstline = open(filename).readline()
        
        # Is the file a GNDS file?
        if firstline.startswith( "<?xml" ) or firstline.startswith("<reactionSuite "):
            import fudge.gnds
            RS = fudge.gnds.reactionSuite.readXML( filename)
            try: 
                CS = fudge.gnds.covariances.covarianceSuite.readXML( filename.replace('.gnds.', '.gndsCov.'))
            except: 
                CS = fudge.gnds.covariances.covarianceSuite.covarianceSuite()
            return [ RS, CS ]

        # Maybe its an ENDF file?
        elif    firstline.endswith(' 0  0    0\n') or \
                firstline.endswith(' 0  0    0\r') or \
                firstline.endswith(' 0  0    0\r\n') or \
                firstline.endswith(' 0  0\n') or \
                firstline.endswith(' 0  0\r') or \
                firstline.endswith(' 0  0\r\n') or \
                filename.endswith('.endf'):
            from fudge.legacy.converting import endfFileToGNDS
            results = endfFileToGNDS.endfFileToGNDS( filename,
                                                   toStdOut=verbose,
                                                   skipBadData=skipBadData,
                                                   reconstructResonances=reconstructResonances,
                                                   continuumSpectraFix=continuumSpectraFix,
                                                   doCovariances=not skipCovariances,
                                                   verboseWarnings=verboseWarnings,
                                                   printBadNK14=printBadNK14,
                                                   ignoreBadDate=ignoreBadDate,
                                                   acceptBadMF10FissionZAP=acceptBadMF10FissionZAP)

            if type( results ) == dict:
                return [ results['reactionSuite'], results['covarianceSuite'] ]
            elif type( results ) == tuple:
                return [ results[0], results[1] ]
            else:
                raise TypeError( "endfFileToGNDS.endfFileToGNDS() returned a "+str(type(results))+", I don't know what to do with it" )

        # Failed!
        else: print "WARNING: Unknown file type, not reading %s"% filename
