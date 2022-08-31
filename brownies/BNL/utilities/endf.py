observableNameMap = {
    3: "Cross Section ( Sigma(E) )",
    23: "Cross Section ( Sigma(E) )",
    4: "Angular Distribution ( dSigma(E)/d(cos Theta) )",
    5: "Energy Distribution ( dSigma(E)/d(E') )",
    6: "Energy-Angle Distribution ( dSigma(E)/d(E')d(cos Theta) )"}

ENDFStandardsReactions = (('H1', 2), ('C0', 2), ('AU197', 102), ('U235', 18), ('U238', 18), ('U238', 102),
                          ('PU239', 18))


def isStandards(projectile, target, mt):
    if projectile != 'n':
        return False
    if mt not in [2, 102, 18]:
        return False
    if target not in ['H1', 'C0', 'AU197', 'U235', 'U238', 'PU239']:
        return False
    for x in ENDFStandardsReactions:
        if (target, mt) == x:
            return True
    return False


def hasMT(reactionSuite, mt):
    return mt in getEvaluationMTs(reactionSuite)


def getPlottableReactions(reactionSuite):
    return reactionSuite.reactions + reactionSuite.summedReactions + reactionSuite.fissionComponents


def getEvaluationMTs(reactionSuite, mtFilter=None):
    reactionList = getPlottableReactions(reactionSuite)
    if mtFilter is None:
        return [r.ENDF_MT for r in reactionList]
    else:
        return [r.ENDF_MT for r in reactionList if r.ENDF_MT in mtFilter]


def getReactions(reactionSuite, MT):
    return [r for r in getPlottableReactions(reactionSuite) if r.ENDF_MT == MT]


def getPointwiseCrossSection(reac):  # FIXME remove? logic is broken + appears to be unused
    if 'linear' in reac.crossSection.forms:
        return reac.crossSection['linear']
    elif 'pointwise' in reac.crossSection.forms:
        return reac.crossSection['pointwise']
    elif 'piecewise' in reac.crossSection.forms:
        return reac.crossSection['piecewise'].toPointwise_withLinearXYs(accuracy=1e-3, upperEps=1e-8)
    elif 'resonancesWithBackground' in reac.crossSection.forms:  # in case --noReconstruct option is used:
        return reac.crossSection['resonancesWithBackground'].background.toPointwise_withLinearXYs(accuracy=1e-3,
                                                                                                  upperEps=1e-8)
    return None


def getUncertainty(theCovariance, theData):
    """ extract absolute uncertainty vector (ie, in units of barn) from covariance """
    from fudge.covariances.covarianceMatrix import CovarianceMatrix
    from fudge.covariances.mixed import MixedForm
    theUncert = None
    for form in (CovarianceMatrix.moniker, MixedForm.moniker):
        if form in theCovariance.forms:
            theUncert = theCovariance.forms[form].getUncertaintyVector(theData, relative=False)
    if theUncert is None:
        print('Cannot plot uncertainty in any of these forms: ' + str([f for f in theCovariance.forms]))
        return
    return theUncert


def getAngularDistribution(reac, product):
    results = []
    for p in reac.outputChannel.particles:
        if p.label != product:
            print("    Skipping distributions for", p)
        else:
            print("    Processing distributions for", p, "    (" + str(p.distributions) + ")")
            for component in p.distributions.components.values():
                compName = str(component).split('/')[-1]
                if compName != "angular":
                    print('        Skipping ' + compName + ' distribution component     (' + str(component) + ')')
                    if compName == 'uncorrelated':
                        print('*** Check uncorrelated data -- we might be able to plot it too')
                    continue
                else:
                    print('        Processing ' + compName + ' distribution component     (' + str(component) + ')')
                    print('        Native form is', component.nativeData)
                    results.append(component.forms[component.nativeData])
                    break
    return results
