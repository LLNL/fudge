import os


def read_columnar_data(filename, number_columns=0, comment_character='#'):
    """
    Read in a file containing variable number of columns
    Lines beginning with comment_character are ignored
    """
    if not os.path.exists(filename):
        raise IOError("Data file %s not found" % filename)
    results = []
    with open(filename) as datafile:
        for line in datafile.readlines():
            try:
                if line.startswith(comment_character) or line.strip() == '':
                    continue
                results.append(list(map(float, line.split()[0:number_columns])))
            except ValueError:
                pass
    return results


def readXYData(filename, comment_character='#'):
    """
    Read in a file containing 2 columns of x, y, dy
    Lines beginning with commentCharacter are ignored
    """
    return read_columnar_data(filename, number_columns=2, comment_character=comment_character)


def readXYdYData(filename, comment_character='#'):
    """
    Read in a file containing 3 columns of x, y, dy
    Lines beginning with commentCharacter are ignored
    """
    return read_columnar_data(filename, number_columns=3, comment_character=comment_character)


def readXdXYdYData(filename, comment_character='#'):
    """
    Read in a file containing 4 columns of x, dx, y, dy
    Lines beginning with commentCharacter are ignored
    """
    return read_columnar_data(filename, number_columns=4, comment_character=comment_character)




def readEvaluation(filename, targ=None, proj=None, verbose=True, skipBadData=True,
                   reconstructResonances=True, continuumSpectraFix=False, skipCovariances=False,
                   verboseWarnings=True, printBadNK14=False, ignoreBadDate=True, acceptBadMF10FissionZAP=True):
    """
    Read in an evaluation in either Fudge/GNDS, ENDF or AMPX/BOF format and return result as Fudge classes
    """
    import fnmatch

    # OK, try AMPX first
    if fnmatch.fnmatch(filename, '*.ampx*'):
        try:
            import ampx2fudge
            import ampx
            import fudge
            ampx_za, bounds = ampx.readEvaluation(filename, str(targ), str(proj))
            return [
                ampx2fudge.convertAmpxNuclideToFudgeReactionSuite(ampx_za, bounds),
                fudge.covariances.covarianceSuite.covarianceSuite(None, None, None)]
        except ImportError:
            print("WARNING: Could not load AMPX module.  Is it in your path?")

    else:
        firstline = open(filename).readline()

        # Is the file a GNDS file?
        if firstline.startswith("<?xml") or firstline.startswith("<reactionSuite "):
            import fudge.reactionSuite
            import fudge.covariances
            RS = fudge.reactionSuite.readXML(filename)
            try:
                CS = fudge.covariances.covarianceSuite.readXML(filename.replace('.gnds.', '.gndsCov.'))
            except:
                CS = fudge.covariances.covarianceSuite.covarianceSuite(None, None, None)
            return [RS, CS]

        # Maybe its an ENDF file?
        elif firstline.endswith(' 0  0    0\n') or \
                firstline.endswith(' 0  0    0\r') or \
                firstline.endswith(' 0  0    0\r\n') or \
                firstline.endswith(' 0  0\n') or \
                firstline.endswith(' 0  0\r') or \
                firstline.endswith(' 0  0\r\n') or \
                filename.endswith('.endf'):
            from brownies.legacy.converting import endfFileToGNDS
            results = endfFileToGNDS.endfFileToGNDS(filename,
                                                    toStdOut=verbose,
                                                    skipBadData=skipBadData,
                                                    reconstructResonances=reconstructResonances,
                                                    continuumSpectraFix=continuumSpectraFix,
                                                    doCovariances=not skipCovariances,
                                                    verboseWarnings=verboseWarnings,
                                                    printBadNK14=printBadNK14,
                                                    ignoreBadDate=ignoreBadDate,
                                                    acceptBadMF10FissionZAP=acceptBadMF10FissionZAP)

            if type(results) == dict:
                return [results['reactionSuite'], results['covarianceSuite']]
            elif type(results) == tuple:
                return [results[0], results[1]]
            else:
                raise TypeError("endfFileToGNDS.endfFileToGNDS() returned a " + str(
                    type(results)) + ", I don't know what to do with it")

        # Failed!
        else:
            print("WARNING: Unknown file type, not reading %s" % filename)
