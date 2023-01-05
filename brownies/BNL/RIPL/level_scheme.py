import copy, os, pprint, math, traceback, sys, collections
from pqu.PQU import PhysicalQuantityWithUncertainty as PhysicalQuantity
from brownies.legacy.endl.structure import masses
from fudge.resonances.common import Spin, Parity
from PoPs.chemicalElements.misc import symbolFromZ as elementSymbolFromZ
from PoPs.chemicalElements.misc import idFromZA as nucleusNameFromZA
from PoPs.chemicalElements.misc import ZAInfo_fromString
from brownies.BNL.plot_evaluation.plotio import readEvaluation
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule


ENDFFILEPATH = os.environ['HOME'] + "/Sites/endfb71.dev/neutrons"

neutronMass = PhysicalQuantity(masses.massTable[1], 'amu')
protonMass = PhysicalQuantity(masses.massTable[1001], 'amu')
quantityDeterminationMap = {'c': 'chosen', None: 'unknown', '': 'unknown', 'u': 'uncertain', 'n': 'unique',
                            'g': 'wild guess'}
cSquared = PhysicalQuantity(1.00000000000000000000,'c**2')
hbar = PhysicalQuantity("6.58211928e-16","eV*s")
hbarc = PhysicalQuantity("197.32697","MeV*fm")
e2 = PhysicalQuantity("1.43998e-10","keV*cm")
muN2 = PhysicalQuantity("1.59234e-38","keV*cm**3")
b = PhysicalQuantity('1e-24','cm**2')


def double_factorial(n):
    """
    Calculate n!!

    :param n: an integer
    :return: n!!
    """
    if not isinstance(n,int): raise TypeError('n must be an int')
    product = 1
    if n%2 == 1: stop=1
    else: stop=0
    for i in range(n,stop,-2): product *= i
    return product


def is_transition_allowed( multipolarity, startJ, startPi, stopJ, stopPi ):
    """
    Test whether a transition of given multipolarity is allowed between the given states' JPi

    :param multipolarity: multipolarity is a tuple of form (X,L) where X='E' or 'M' and L is an int
    :param startJ: initial state spin
    :param startPi: inital state parity
    :param stopJ: final state spin
    :param stopPi: final state parity
    :return: if transition is allow
    """
    if multipolarity[0] not in ["E","M"] or type(multipolarity[1])!=int:
        raise ValueError( "multipolarity be of form (X,L) where X='E' or 'M' and L is an int" )
    if str(startPi)==str(stopPi):
        PiTest=1
    else:
        PiTest=-1
    if isinstance(startJ, Spin):
        startJ=startJ.value
    if isinstance(stopJ, Spin):
        stopJ=stopJ.value
    if multipolarity[0] == "E" and PiTest != pow(-1,multipolarity[1]):
        return False
    elif multipolarity[0] == "M" and PiTest != pow(-1,multipolarity[1]+1):
        return False
    return abs(startJ-stopJ)<=multipolarity[1] and startJ+stopJ>=multipolarity[1]


def WeisskopfSingleParticleEstimateBXL(multipolarity,A):
    """
    Computes Weisskopf's single particle for B(XL)

    :param multipolarity: multipolarity is a tuple of form (X,L) where X='E' or 'M' and L is an int
    :param A: number of nucleons in nucleus in question
    :return: the estimate of B(XL), it is dimensionless
    """
    if multipolarity[0] not in ["E","M"] or type(multipolarity[1])!=int:
        raise ValueError( "multipolarity be of form (X,L) where X='E' or 'M' and L is an int" )
    L=multipolarity[1]
    R=1.2e-13*pow(A,1./3.) # This is nuclear radius in cm
    x = R*R/1e-24  # This is R^2/b where b = 10^-24 cm
    BXL = 3.0/(3.0+L)
    BXL *= BXL/math.pi
    if multipolarity[0]=='E':
        return BXL*pow(x,L)/4.0
    else:
        return BXL*10.0*pow(x,L-1)


def gammaHalflife(multipolarity,BXL,Eg):
    if multipolarity[0] not in ["E","M"] or type(multipolarity[1])!=int:
        raise ValueError( "multipolarity be of form (X,L) where X='E' or 'M' and L is an int" )
    if not isinstance(Eg,PhysicalQuantity): raise TypeError("Eg must be a PhysicalQuantity")
    L=multipolarity[1]
    prefactor = math.log(2)*L*pow(double_factorial(2*L+1),2.0)*hbar*pow(hbarc/Eg,2*L+1)/8.0/math.pi/(L+1)/BXL
    if multipolarity[0]=='E':
        return prefactor/e2/pow(b,L)
    else:
        return prefactor/muN2/pow(b,L-1)


# ---------------------------------------------------------------------------------
#
# Basic data types
#
# ---------------------------------------------------------------------------------

class nucleus:
    '''
    Simple nucleus implementation, capturing (I hope) all the common elements from the RIPL, TALYS and GNDS formats.
    Member data from xParticle::

        - name : GNDS name
        - genre : GNDS genre
        - attributes : GNDS XML attributes dict.  One of these has key 'JPi' with value of Jpi corresponding to the nucleus' ground state or None
        - levels : GNDS level dict.  Key is the level index.  Level "0" is the ground state.  getJPi() will evaluate to this if attributes[ 'JPi' ] == None

    Other member data added here::

        - Z
        - A
        - ZA
        - suffix : GNDS suffix
        - symbol
        - mass : mass as PhysicalQuantity
        - Sn : neutron separation energy as PhysicalQuantity
        - Sp : proton separation energy as PhysicalQuantity
        - levelDensity : a level density instance
        - resonanceSpacingSWave : D0  (see below)
        - uncResonanceSpacingSWave : dD0  (see below)
        - resonanceSpacingPWave : D1  (see below)
        - uncResonanceSpacingPWave : dD1  (see below)
    '''

    def __init__(self, **kw):
        '''
        Easiest is just to call with GNDS name as sole argument like this:

            >>> nucleus( name="C12" )
        '''
        self.name=kw['name']
        self.attributes=kw.get('attributes', {})
        self.Z, self.A, self.ZA, lev = ZAInfo_fromString(self.name)
        self.mass=kw.get('mass', None)
        if self.mass is None:
            try:
                # i don't understand fudge mass setting, it is very confusing
                self.mass = PhysicalQuantity(masses.massTable[self.ZA], 'amu')
            except KeyError:
                self.mass = None
        try:
            self.Sn = (neutronMass + PhysicalQuantity(masses.massTable[self.ZA - 1], 'amu') - self.mass).inUnitsOf(
                'MeV/c**2') * cSquared
        except (KeyError, TypeError):
            self.Sn = None
        try:
            self.Sp = (protonMass + PhysicalQuantity(masses.massTable[self.ZA - 1001], 'amu') - self.mass).inUnitsOf(
                'MeV/c**2') * cSquared
        except (KeyError, TypeError):
            self.Sp = None

        self.symbol = elementSymbolFromZ[self.Z]
        self.lastLevelInCompleteScheme = None
        self.lastLevelWithAssignedJPi = None
        self.levelDensity = kw.get('levelDensity', None)
        self.resSpacingSWave = kw.get('D0', None)
        self.uncResSpacingSWave = kw.get('dD0', None)
        self.resSpacingPWave = kw.get('D1', None)
        self.uncResSpacingPWave = kw.get('dD1', None)
        self.endfMap = {}
        self.endfFile=None
        self.resonances = None
        self.BInWeisskopfUnits={}
        self.BInWeisskopfUnits[('E',1)]=kw.get( 'BE1', None )
        self.BInWeisskopfUnits[('E',2)]=kw.get( 'BE2', None )
        self.BInWeisskopfUnits[('M',1)]=kw.get( 'BM1', None )
        self.BInWeisskopfUnits[('E',3)]=kw.get( 'BE3', None )
        self.BInWeisskopfUnits[('M',2)]=kw.get( 'BM2', None )
        self.BInWeisskopfUnits[('E',4)]=kw.get( 'BE4', None )
        self.BInWeisskopfUnits[('M',3)]=kw.get( 'BM3', None )
        self.BInWeisskopfUnits[('E',5)]=kw.get( 'BE5', None )
        self.BInWeisskopfUnits[('M',4)]=kw.get( 'BM4', None )
        self.__levels=[]

    def addLevel(self, lev):
        self.__levels.append(lev)

    @property
    def levels(self):
        return self.__levels

    def readResonanceData(self, endfFile=None, verbose=False, keepGoing=True):
        """
        Reads the resonances from an ENDF File in which the compound nucleus formed
        by the target+projectile matches the nucleus in self.

        :param endfFile: the filename of the evaluation to load.  If None, get from default area
        :param verbose:
        :param keepGoing: if False, throw exception if endfFile not found
        :returns: None
        """
        if endfFile is None: endfFile = ENDFFILEPATH + os.sep + 'n-%s_%s_%s.endf' % (
            str(self.Z).zfill(3), elementSymbolFromZ(self.Z), str(self.A - 1).zfill(3))
        if not os.path.exists(endfFile):
            if keepGoing:
                print ('evaluation %s not found' % endfFile)
            else: raise RuntimeWarning('evaluation %s not found' % endfFile)
        else:
            endfMap = readEvaluation(endfFile, verbose=verbose, skipBadData=True)
            proj = endfMap[0].projectile.getZ_A_SuffixAndZA()
            targ = endfMap[0].target.getZ_A_SuffixAndZA()
            if proj[0] + targ[0] != self.Z or proj[1] + targ[1] != self.A:
                raise ValueError(
                    "target (%s) + projectile (%s) system doesn't form the requested compound system (%s)" % (
                        str(endfMap[0].projectile), str(endfMap[0].target),
                        str(nucleusNameFromZA(self.Z * 1000 + self.A)) ))
            self.endfFile = endfFile
            self.endfMap = endfMap
            if hasattr(self.endfMap[0], 'resonances'): self.resonances = endfMap[0].resonances

    def addResonancesToLevels(self, verbose=False):
        if self.resonances is None: return
        if self.resonances.resolved is None: return

        # Set up region lists
        if self.resonances.resolved.multipleRegions:
            regions = self.resonances.resolved.regions
        else:
            regions = [self.resonances.resolved]

        iLevel=self.getNumberLevels()-1
        iRes=-1

        # Process each region
        for RRidx, RR in enumerate(regions):
            if RR.nativeData.moniker in [ 'SingleLevel_BreitWigner', 'MultiLevel_BreitWigner', 'Reich_Moore' ]:
                if verbose:
                    pprint.pprint(RR.nativeData.resonanceParameters.data)
                for level in RR.nativeData.resonanceParameters.data:
                    iLevel += 1
                    iRes+=1
                    J = Spin(str(int(abs(level[2])*2.0))+'/2')
                    if level[2] == 0.0: Pi = Parity(-1)
                    else: Pi = Parity(+1)
                    self.addLevel(
                        nuclearLevel(
                            name=self.symbol+str(self.A)+'_r%i'%(iRes),
                            index=iLevel,
                            energy=(PhysicalQuantity(level[0],'eV')+self.Sn).inUnitsOf('MeV'),
                            spin=J,
                            parity=Pi))
            elif RR.nativeData.moniker == 'R_Matrix_Limited':
                for iSG, SG in enumerate(RR.nativeData.spinGroups):
                    if verbose:
                        print ("spin group %i (J=%s, Pi=%s)"%(iSG,str(SG.spin),str(SG.parity)))
                    if verbose:
                        pprint.pprint(SG.resonanceParameters.data)
                    for level in SG.resonanceParameters.data:
                        iLevel += 1
                        iRes+=1
                        self.addLevel(
                            nuclearLevel(
                                name=self.symbol+str(self.A)+'_r%i'%(iRes),
                                index=iLevel,
                                energy=(PhysicalQuantity(level[0],'eV')+self.Sn).inUnitsOf('MeV'),
                                spin=SG.spin,
                                parity=SG.parity))

        if self.resonances.unresolved is not None:
            pass

    def getNumberLevels(self):
        return len(self.levels)

    def getNumberGammas(self):
        return sum([self.levels[x].getNumberGammas() for x in self.levels])

    def getNumberDecays(self):
        return sum([self.levels[x].getNumberDecays() for x in self.levels])

    def getAllowedEMTransitions(self, iStart, multipolarity=None, BInWeisskopfUnits={}, BRmin=None):
        """
        return all possible EM transitions out of level iStart for a given multipolarity.
        If multipolarity is not given, this routine will do E1, E2 and M1 only

        :param iStart: the level index
        :param multipolarity: request a given multipolarity, as a string e.g. "E2"
        :param BInWeisskopfUnits:
        :param BRMin:  minimum BR to consider.  For example, setting BRmin=0.01 means no
                       BR's smaller than 1% will be generated.  Bigger BR's will be rescaled to sum to 1.0
        :return: a list of all possible EM transitions for a given multipolarity
        """
        #BInWeisskopfUnits.update(self.BInWeisskopfUnits)

        # Assemble all the multipolarities to consider
        if multipolarity is not None: multipolarityList = [multipolarity]
        else: multipolarityList=[('E',1), ('E',2), ('M',1)]#, ('E',3), ('M',2)] # defaults

        # A copy of current level to add gammas to
        dummyLevel = copy.copy(self.levels[iStart])
        dummyLevel.gammas=[]

        # Compute all allowed gammas
        for multipolarity in multipolarityList:
            for i in range(0,iStart):
                if is_transition_allowed(
                        multipolarity,
                        self.levels[iStart].spin,
                        self.levels[iStart].parity,
                        self.levels[i].spin,
                        self.levels[i].parity):

                    # Compute BXL
                    BXL=BInWeisskopfUnits[multipolarity]*WeisskopfSingleParticleEstimateBXL(multipolarity,self.A)

                    # Gamma energy
                    Eg = self.levels[iStart].energy-self.levels[i].energy

                    # Compute T1/2 here
                    Thalf=gammaHalflife(multipolarity,BXL,Eg)
                    try:
                        Thalf.convertToUnit('ns')
                    except TypeError as err:
                        raise TypeError( str(multipolarity) )

                    # Create a new gamma instance
                    g=nuclearLevelGamma(energy=Eg, finalLevel=self.levels[i], halflife=Thalf, multipolarity=multipolarity)

                    # Add it to list
                    dummyLevel.addGamma(g)

        # Sort the gammas
        dummyLevel.gammas.sort(cmp=lambda x,y:cmp(y.finalLevel.name,x.finalLevel.name))

        # Fix BRs
        if len(dummyLevel.gammas)!=0:
            Gtot = 1.0/dummyLevel.gammas[0].halflife
            for g in dummyLevel.gammas[1:]: Gtot+=1.0/g.halflife
            for g in dummyLevel.gammas:
                g.setEmissionProbabilities(Pg=1.0/(Gtot*g.halflife),ICC=0.0,Pe=0.0,nonRadiativeProbability=0.0)

        # Filter out gammas with small BR's, then renormalize
        if BRmin is not None:
            newGammas = []
            BRTotal = 0.0
            for g in dummyLevel.gammas:
                if g.gammaEmissionProbability >= BRmin:
                    newGammas.append(g)
                    BRTotal += g.gammaEmissionProbability
            for g in newGammas:
                g.gammaEmissionProbability/=BRTotal
            dummyLevel.gammas = newGammas


        return dummyLevel

    def getCumulativeLevelDistribution(self, J=None, Pi=None):
        myAxes = axesModule.Axes(labelsUnits={0: ('Excitation energy', 'MeV'), 1: ('Number of levels', '')})
        theData = [[0.0, 1.0]]
        def JPiMatches(level):
            Jmatch = J is None or str(level.spin)=='?'
            Pimatch = Pi is None or str(level.parity)=='?'
            if Jmatch or Pimatch: return True
            return int(float(level.spin)-J)<0.001 and str(level.parity)==Pi
        levelEnergyList = [float(self.levels[level].energy.inUnitsOf('MeV').value) for level in self.levels if JPiMatches(level)]
        levelEnergyList.sort()
        for thisEnergy in levelEnergyList[1:]:
            if abs(thisEnergy - theData[-1][0]) < 1e-6 * thisEnergy:
                theData[-1][1] = theData[-1][1] + 1  # Deal with corner case where two levels are degenerate
            else:
                theData.append([thisEnergy, theData[-1][1] + 1])
        return XYs1dModule.XYs1d(interpolation='flat',axes=myAxes, data=theData, accuracy=1e-5)

    def getLevelSpacingDistribution(self):
        myAxes = axesModule.defaultAxes(dimension=2, independentInterpolation='linear', dependentInterpolation='flat',
                                  labelsUnits={0: ( 'Excitation energy', 'MeV' ), 1: ( 'Level spacing, D', 'MeV' )})
        theData = []
        levelList = self.levels.keys()
        levelList.sort()
        thisEnergy = float(self.levels[levelList[0]].energy.inUnitsOf('MeV').value)
        for level in levelList[1:]:
            lastEnergy = copy.copy(thisEnergy)
            thisEnergy = float(self.levels[level].energy.inUnitsOf('MeV').value)
            theData.append([thisEnergy, thisEnergy - lastEnergy])
        return XYs1dModule.XYs(myAxes, data=theData, accuracy=1e-5)

    def getLastLevelInCompleteSchemeIndex(self):
        if self.lastLevelInCompleteScheme == None: raise NotImplementedError()
        return self.lastLevelInCompleteScheme

    def getLastLevelWithAssignedJPiIndex(self):
        if self.lastLevelWithAssignedJPi == None: raise NotImplementedError()
        return self.lastLevelWithAssignedJPi

    def getJPiList(self):
        result = []
        for x in self.levels.values():
            theSpin=x.spin
            theParity=x.parity
            if type(theSpin)==float:
                theSpin=Spin(str(int(2.0*x.spin))+'/2')
            if type(theParity)==int:
                theParity=parity(x.parity)
            tmp = (theSpin, theParity)
            if not tmp in result: result.append(tmp)
        result.sort()
        return result

    def statusString(self):
        return "Nucleus: " + str(self.name) + \
               ", ZA = " + str(self.ZA) + \
               ", mass = " + str(self.mass) + \
               ", Sn = " + str(self.Sn) + \
               ", Sp = " + str(self.Sp) + \
               ', numLevels = ' + str(self.getNumberLevels()) + \
               ', D0 = ' + str(self.resSpacingSWave) + ' +/- ' + str(self.uncResSpacingSWave) + \
               ', D1 = ' + str(self.resSpacingPWave) + ' +/- ' + str(self.uncResSpacingPWave)

    def levelReport(self, showGammas=False):
        belowSn=self.Sn is not None
        belowSp=self.Sp is not None
        levs = self.levels
        levs.sort(key=lambda x: x.energy)
        result =[ "           index         name              E_level     J Pi   Resonance?" ]
        for ilev,lev in enumerate(levs):
            if belowSn and lev.energy>=self.Sn:
                result.append( (8*' ')+' '+(56*"-")+" Sn = "+str(self.Sn) )
                belowSn=False
            if belowSp and lev.energy>=self.Sp:
                result.append( (8*' ')+' '+(56*"-")+" Sp = "+str(self.Sp) )
                belowSp=False
            result.append( '            '+str(ilev).rjust(4)+' '+lev.levelReport(showGammas=showGammas))
            if '_r' in lev.name: result[-1]+='   *'
        return '\n'.join(result)

    def levelReportMathematica(self,showGammas=False,maxLevelIndex=41,style='other',modifiedLevels=[17,26,27,37,39]):
        '''
        Before cut-n-pasting, do:
                <<"LevelScheme`"
        or
                Get["LevelScheme`"]
        :param showGammas:
        :param maxLevelIndex:
        :param style:
        :param modifiedLevels:
        :return:
        '''
        def parityTranslator(parity):
            if str(parity)=='+': return '+1'
            elif str(parity)=='-': return '-1'
            return 'None'
        def spinTranslator(spin):
            # FIXME: code logic for use of DiagonalFractionBox[a,b]
            return str(spin)
        def getColor(i):
            if i in modifiedLevels: return 'Red'
            return 'Black'
        def levelKey(name): return str(name).replace('_','')
        def levelFormatter(_lev,_leftHorPos,_rightHorPos,_yScale=1.0):
            spin=spinTranslator(_lev.spin)
            parity=parityTranslator(_lev.parity)
            _ilev=int(str(_lev.name).split('_e')[-1])
            return 'Lev[%s,%f,%f,%f,LabL -> LabelJP["%s", %s], Color->%s] ' % ( levelKey(_lev.name), _leftHorPos, _rightHorPos, _lev.energy.value*_yScale, spin, parity,getColor(_ilev) )
        def mathematicaCommand(_objects,A,symbol,xmin,xmax,ymin,ymax,fontsize=10,pageSizeX=None,pageSizeY=None):
            if pageSizeX is None: pageSizeX=xmax-xmin
            if pageSizeY is None: pageSizeY=3*(ymax-ymin)
            return '''Get["LevelScheme`"]
Figure[
    {
        ManualLabel[{7, -0.1}, "\!\(\*SuperscriptBox[\(\[InvisiblePrefixScriptBase]\), \(%A\)]\)%SYM", FontSize -> 35],
        %OBJECTLIST
    },
    PlotRange -> {{%XMIN, %XMAX}, {%YMIN, %YMAX}},
    ImageSize -> 72*{%SIZEX, %SIZEY}
]'''\
                .replace("%OBJECTLIST",',\n'.join(_objects))\
                .replace("%FONTSIZE", str(fontsize))\
                .replace("%XMIN",str(xmin))\
                .replace("%XMAX",str(xmax))\
                .replace("%YMIN",str(ymin))\
                .replace("%YMAX",str(ymax))\
                .replace("%SIZEX",str(pageSizeX))\
                .replace("%SIZEY",str(pageSizeY))\
                .replace("%A",str(A))\
                .replace("%SYM",symbol)
        def jPiAsNumber(_lev):
            try:
                jPi = float(_lev.spin.value)
            except TypeError:
                jPi = 0.0
            if str(_lev.parity) not in ['+','?']:
                jPi*=-1
            return jPi
        def deltaJPi(_levs):
            jPiList = []
            for _ilev,_lev in _levs: jPiList.append(jPiAsNumber(_lev))
            return max(jPiList)-min(jPiList)
        def maxEnergy(_levs):
            _maxE=0.0
            for _ilev,_lev in _levs[0:maxLevelIndex]:
                _maxE=max(_maxE,_lev.energy.value)
            return _maxE
        def sign(x):
            if x<0.0: return -1
            return +1

        # Get a sorted list of levels
        levs = self.levels.items()
        levs.sort(cmp=lambda x,y: cmp(x[1].energy, y[1].energy))

        if style=='betaDecay':
            # Style setup
            objects = []
            objects.append( "SetOptions[Trans, BackgroundT -> Automatic, NudgeT -> 2, Thickness -> 1, FontSize -> %FONTSIZE, OrientationT -> Automatic, FillColor -> LightGray]" )
            objects.append( "    SetOptions[Lev, Thickness -> 2, LabR -> Automatic, WingRiseWidth -> 10, WingTipWidth -> 10, FontSize -> %FONTSIZE, Color -> Black]" )
            maxE = maxEnergy(levs)
            lineWidth=self.getNumberGammas()/30.
            # Do levels first
            for ilev,lev in levs[0:maxLevelIndex]: objects.append( '        '+levelFormatter(lev,0.,lineWidth) )
            # Now set up the transitions
            if showGammas:
                objects.append("        AutoLevelInit[%f, -0.1, -0.1]" % (lineWidth*0.925) )
                for ilev,lev in levs[0:maxLevelIndex]:
                    if len(lev.gammas)>0:  objects.append( '            '+"AutoLevel[%s]"% levelKey(lev.name) )
                    for gam in lev.gammas: objects.append( '            '+'AutoTrans[%s,LabT->"%s",Color->%s]' % ( levelKey(gam.finalLevel), str(gam.energy), getColor(ilev)) )
            return mathematicaCommand(objects,self.A,self.symbol,0.,lineWidth*1.3,-0.5,maxE+0.5)
        else:
            objects=[]
            objects.append( "SetOptions[Trans, BackgroundT -> Automatic, NudgeT -> 2, Thickness -> 1, FontSize -> %FONTSIZE, OrientationT -> Automatic, ArrowType->ShapeArrow, FillColor -> LightGray]" )
            objects.append( "    SetOptions[Lev, Thickness -> 1, LabR -> Automatic, FontSize -> %FONTSIZE, Color -> Black]" )
            maxE = maxEnergy(levs)
            yScale=3.0
            xMin=0.0
            xMax=0.0
            xGap=0.1
            xWidth=5.0
            # Do levels first
            for ilev,lev in levs[0:maxLevelIndex]:
                jPi=jPiAsNumber(lev)
                x1=jPi*(xGap+xWidth)
                x2=jPi*(xGap+xWidth)+sign(jPi)*xWidth
                xStart=min(x1,x2)
                xStop=max(x1,x2)
                xMin=min(xMin,xStart)
                xMax=max(xMax,xStop)
                objects.append( '        '+levelFormatter(lev,xStart,xStop) )
                if showGammas:
                    for gam in lev.gammas:
                        gamXStart=xStart #0.5*(xStart+xStop)
                        jPiFinal=jPiAsNumber(gam.finalLevel)
                        gamXStop=jPiFinal*(xGap+xWidth) #+0.5*sign(jPiFinal)*xWidth
                        objects.append("            Trans[%s,%s,Width->%f,Color->%s]"%(levelKey(lev.name),levelKey(gam.finalLevel),2.0*xWidth*gam.emTransitionProbability,getColor(ilev)))
            return mathematicaCommand( objects, self.A, self.symbol, xmin=(xMin-xGap), xmax=(xMax+xGap), ymin=-0.7, ymax=(maxE+0.5), pageSizeX=deltaJPi(levs), pageSizeY=yScale*maxE )

    def toRIPLFormat(self):
        '''
        The format consists of three types of records: (i) identification,
        (ii) level, and (iii) gamma. Each isotope begins with an
        identification record such as the one shown below for the case of
        Mg-22 (heading has been introduced to facilitate reading but does not
        show up in the actual file):

          SYMB    A     Z     Nol    Nog    Nmax    Nc     Sn[MeV]     Sp[MeV]
          22Mg    22    12    17     18      9      4     19.382000    5.497000

        SYMB   : mass number with symbol of the element
        A      : mass number
        Z      : atomic number
        Nol    : number of levels in the decay scheme
        Nog    : number of gamma rays in the decay scheme
        Nmax   : maximum number of levels up to which the level scheme is
                 complete
        Nc     : number of a level up to which spins and parities are unique
        Sn     : neutron separation energy in MeV
        Sp     : proton separation energy in MeV

        The corresponding FORTRAN format is (a5,6i5,2f12.6)
        '''
        raise NotImplementedError()

    def toTALYSFormat(self):
        '''
        Z, A, total number of lines for this isotope, number of levels for this
        isotope, nuclear symbol with the FORTRAN format
            (2i4,2i5,56x,i4,a2)
        '''
        raise NotImplementedError()

    def __str__(self):
        return self.name


class nuclearLevel:
    '''
    Simple nucleus level implementation, capturing (I hope) all the common elements from the RIPL, TALYS and GNDS formats.
    Member data from nuclearLevel:
        - name : GNDS name of nuclearLevel (usu. something like C12_e3)
        - energy : level energy as PhysicalQuantity
        - index : gets remapped to "label" inside
        - gammas : list
        - attributes : GNDS XML attributes dict.  Two of these keys are 'spin' and 'parity' of JPi
        - groundState : points to xParticle instance
    Other member data added here:
        - lifetime : level lifetime as PhysicalQuantity.  If < 0.0 s, it is stable.  If == 0.0 s, is not stable or is so short it is unknown (the default).
        - decays : list of decays
        - spin : integer or 1/2 integer or None.  We attempt to sync it with the spin entry in the attributes map.
        - parity : either +1, -1 or None. We attempt to sync it with the parity entry in the attributes map.
        - alternateJPiAssignments: free string as one sees in ENSDF e.g "(3/2,5/2,7/2)-"
        - howJPiDetermined: key stating how assignment made (c=chosen, u=uncertain, n=unique, g=guessed, '' or None=unknown)
    '''

    def __init__(self, **kw):
        '''
        Easiest is just to call with GNDS name as sole argument like this:
            nuclearLevel( name="C12_e2", energy=PhysicalQuantity('1.0','MeV'), spin="1/2", parity=-1, index=2, lifetime=PhysicalQuantity('1.0','s') )
        '''
        self.id = kw.get('name')
        self.index = kw.get('index')
        self.energy = kw.get('energy')
        self.lifetime = kw.get('lifetime', PhysicalQuantity(0.0, 's'))
        if self.lifetime.getValue() < 0.0:
            self.lifetime = PhysicalQuantity(2e10, 'yr')  # > age of the universe, so is stable
        self.attributes = kw.get('attributes', {})
        self.spin = self.attributes.get('spin', kw.get('spin', None))
        self.parity = self.attributes.get('parity', kw.get('parity', None))
        self.attributes['spin'] = self.spin
        self.attributes['parity'] = self.parity
        self.howJPiDetermined = quantityDeterminationMap[kw.get('howJPiDetermined', None)]
        self.alternateJPiAssignments = kw.get('alternateJPiAssignments', None)
        self.decayData = kw.get('decayData', []) # also holds gammas I think
        self.__gammas = []

    @property
    def name(self):
        return self.id

    @property
    def gammas(self):
        return self.__gammas

    def __str__(self):
        return self.id

    def renormalizeTransitionProbabilities(self): raise NotImplementedError()

    def getNumberGammas(self): return len(self.gammas)

    def getNumberDecays(self): return len(self.decays)

    def setJPi(self,J,Pi):
        """
        Sets spin and parity

        :param J: is the spin and should be float
        :param Pi: is the parity and should be either +1 or -1
        :return: None
        """
        self.set_spin(J)
        self.set_parity(Pi)

    def set_spin(self,J):
        """
        Sets spin

        :param J: is the spin and should be a float (if 1/2-int or int) or int
        :return: None
        """
        try:
            self.spin=Spin(J)
        except:
            raise TypeError("J must be a str, float or int")

    def set_parity(self,Pi):
        """
        Sets parity

        :param Pi: is the parity and should be either +1 or -1
        :return: None
        """
        if isinstance(Pi,float):Pi=int(Pi)
        if Pi not in [1,-1]: raise ValueError("Pi must be +1 or -1")
        self.parity=Pi

    def addGamma(self, g): self.gammas.append(g)

    def addDecay(self, d): self.decays.append(d)

    def levelReport(self,showGammas=False):
        result = self.name.rjust(12)+' '+str(self.energy).rjust(20)+' '+str(self.spin).rjust(5)+'  '+str(self.parity)
        if showGammas and len(self.gammas)>0:
            result += '\n' \
                      + 50*' ' + ( ''.join( map(\
                            lambda x: str(x).rjust(20),\
                            ( 'Energy (MeV)', 'Final Level', 'Pgamma', "Pnon-gamma", "Multipolarity", "T1/2 (nsec)") ) ) ) + '\n' \
                      + '\n'.join([g.levelReport() for g in self.gammas])
        return result

    def toRIPLFormat(self):
        '''
        The identification record is followed by the level record, which in
        turn is followed by the pertinent gamma records if decay of the level
        is known. This combination is repeated for all levels in a given
        isotope. An example below shows level records for the ground state and
        the first excited state in Mg-22::

            N1  Elv[MeV]  s   p   T1/2    Ng  J  unc  spins   nd  m  percent  mode

            1  0.000000  0.0  1  3.86E+00  0           0+      1  =  100.0000 %EC+%B+
            2  1.246300  2.0  1  2.10E12   1           2+      0


        Each level record may contain the following quantities::

            Nl    :  sequential number of a level
            Elv   :  energy of the level in MeV
            s     :  level spin (unique). Whenever possible unknown spins up to
                     Umax were inferred. Unknown and undetermined spins are entered as -1.0
            p     :  parity (unique). If the parity of the level was unknown, positive or
                     negative was chosen with equal probability. Parities were determined up
                     to Umax as in the case of spins. The method of choice is not coded.
            T1/2  :  half-life of the level (if known). All known half-lives or level widths
                     were converted into seconds. Half-lives of stable nuclei are
                     represented as -1.0E+0.
            Ng    :  number of gamma rays de-exciting the level.
            J     :  flag for spin estimation method.
            unc   :  flag for an uncertain level energy. When impossible to determine, the
                     relative energy of the band, the energy of these band heads were set to
                     0.0 keV or, if the level order is known, to the preceding level energy
                     with a note that an unknown energy X should be added. The notation uses
                     X+, Y+, Z+ etc. for different bands.
            spins :  original spins from the ENSDF file.
            nd    :  number of decay modes of the level (if known). Values from 0 through 10
                     are possible; 0 means that the level may decay via gamma-emission,
                     and other decay modes are not known.
            m     :  decay percentage modifier; informs a user about major uncertainties.
                     The modifiers are copied out of ENSDF with no modification, and can
                     have the following values: =, <, >, ? (unknown, but expected),
                     AP (approximate), GE (greater or equal), LE (less or equal),
                     LT (less then), SY  (value from systematics).
            percent: percentage decay of different decay modes. As a general rule  the
                     various decay modes add up 100%. There are, however, two exceptions:
                     (i) when a small percentage decay is present, the sum may be slightly
                     more then 100% due to rounding error, (ii) when beta -decay is
                     followed by a heavier particle emission, the percentage of the
                     beta-delayed particle emission is given as a portion of the beta-decay
                     and the sum can be substantially larger then 100%. Naturally, when the
                     modifier is "?" the sum is indefinite.
            mode   : short indication of decay modes of a level (see Table below).

        The corresponding FORTRAN format is::

            (i3,1x,f10.6,1x,f5.1,i3,1x,(e10.2),i3,1x,a1,1x,a4,1x,a18,i3,10(1x,a2,1x,f10.4,1x,a7))
        '''
        raise NotImplementedError()

    def toTALYSFormat(self):
        '''
        level number, level energy in MeV, spin, parity, number of gamma-ray branchings,
        lifetime, spin assignment character, parity as- signment character, original
        ENSDF string with the FORTRAN format::

            (i4,f11.6,f6.1,3x,i2,i3,19x,e9.3,1x,2a1,a18)
        '''
        raise NotImplementedError()


class nuclearDecayMode:
    '''
    Simple implementation of a light-weight decay mode indicator.  Unused.
    '''

    def __init__(self, mode=None, probability=0.0):
        self.mode = mode
        self.probability = probability

    def toRIPLFormat(self):
        '''
        Coding of decay modes. Some minor possibilities, such as
        decay through the emission of Ne-20, are neglected.

              -----------------------------------------------------
                     Code    Meaning
              -----------------------------------------------------
                       %B    beta- decay
                      %EC    electron capture
                  %EC+%B+    electron capture and beta+ decay
                       %N    neutron decay
                       %A    alpha  decay
                      %IT    isomeric transition
                       %P    proton decay
                     %3HE    He-3 decay
                     %B+P    beta+ delayed proton decay
                      %BN    beta- delayed neutron decay
                      %SF    spontaneous fission
                     %ECP    electron capture delayed proton decay
                     %ECA    electron capture delayed alpha  decay
                       %G    gamma   decay
                     %B2N    beta- delayed double neutron decay
                    %B+2P    beta+ delayed double proton decay
              ------------------------------------------------------
        '''
        raise NotImplementedError()

    def toTALYSFormat(self): raise NotImplementedError()


class nuclearLevelGamma:
    '''
    Nuclear electromagnetic transition.  This usually is a gamma, hence the name.   It may also include internal conversions.
    Since we are ignoring non-electromagnetic decays, the gamma emission probability + the internal conversion probability = 1.0.

    Member data from nuclearLevelGamma:

        - energy : energy released during transition as a PhysicalQuantity
        - finalLevel : name attribute from a nuclearLevel instance (e.g. "C12_e1")
        - probability : the gamma emission probability ( == gammaEmissionProbability )
        - nonRadiativeProbability : the internal conversion probability ( == internalConversionProbability )
        - angularDistribution : unused

    Other member data added here:

        - gammaEmissionProbability : the gamma emission probability
        - internalConversionProbability : the internal conversion probability
        - emTransitionProbability : the gamma emission probability + the internal conversion probability ( <= 1.0 )
        - internalConversionCoefficient : the internal conversion coefficient ( == internalConversionProbability / gammaEmissionProbability )
        - howBRDetermined : Optional designator to describe how gamma gammaEmissionProbability determined
    '''

    def __init__( self, **kw ):
        '''
        If starting from RIPL, initialize like this:

            >>> nuclearLevelGamma( energy=PhysicalQuantity("1.322 MeV"), finalLevel="Ni58_e1", Pg=0.96, ICC=0.0 )

        If starting from GND, initialize like this:

            >>> nuclearLevelGamma( energy=PhysicalQuantity("1.322 MeV"), finalLevel="Ni58_e1", probability=0.96, nonRadiativeProbability=0.0 )

        If starting from ENDF, initialize like this:

            >>> nuclearLevelGamma( energy=PhysicalQuantity("1.322 MeV"), finalLevel="Ni58_e1", probability=GP*TP, nonRadiativeProbability=(1.0-GP)*TP )
        '''
        self.setEmissionProbabilities(
            Pg=kw.get( 'Pg', kw.get( 'probability' ) ),
            ICC=kw.get( 'ICC' ),
            Pe=kw.get( 'Pe' ),
            nonRadiativeProbability=kw.get( 'nonRadiativeProbability', 0.0 ))
        self.finalLevel=kw.get('finalLevel')
        #                                               angularDistribution=kw.get('angularDistribution', None),
        #                                               probability=self.gammaEmissionProbability,
        #                                               nonRadiativeProbability=self.internalConversionProbability)

        #self.setAncestor(kw.get('startLevel'))

        # This should really be computed from the levels, not specified here
        self.energy = kw.get('energy')

        # Optional designator to describe how gamma Pg determined
        self.howBRDetermined        = quantityDeterminationMap[ kw.get( 'howBRDetermined', None ) ]
        self.multipolarity = kw.get( 'multipolarity', (None,None) )
        self.halflife = kw.get( 'halflife', '?' )

    def setEmissionProbabilities(self, Pg=None, ICC=None, Pe=None, nonRadiativeProbability=None):
        """
        sets the emission probabilities
        """
        self.gammaEmissionProbability = Pg
        self.internalConversionCoefficient = ICC
        self.emTransitionProbability = Pe
        self.internalConversionProbability = nonRadiativeProbability # probably will get overwritten below for consistency
        if self.emTransitionProbability == None and self.internalConversionCoefficient != None:
            self.emTransitionProbability = self.gammaEmissionProbability * ( 1.0 + self.internalConversionCoefficient )
            self.internalConversionProbability = self.gammaEmissionProbability * self.internalConversionCoefficient
        elif self.emTransitionProbability != None and self.internalConversionCoefficient == None:
            self.internalConversionProbability = self.emTransitionProbability - self.gammaEmissionProbability
            self.internalConversionCoefficient = self.internalConversionProbability / self.gammaEmissionProbability
        elif self.emTransitionProbability == None and self.internalConversionCoefficient == None:
            self.emTransitionProbability = self.gammaEmissionProbability
            self.internalConversionProbability = 0.0
            self.internalConversionCoefficient = 0.0
        elif self.internalConversionProbability != 0.0:
            self.emTransitionProbability = self.gammaEmissionProbability + self.internalConversionProbability
            self.internalConversionCoefficient = self.internalConversionProbability / self.gammaEmissionProbability
        else:
            self.internalConversionProbability = self.emTransitionProbability - self.gammaEmissionProbability
            # should we check the internalConversionCoefficient?

        if nonRadiativeProbability is None:
            self.nonRadiativeProbability = 1.0-self.gammaEmissionProbability
        else:
            self.nonRadiativeProbability = nonRadiativeProbability

        # Some checks:
        if self.emTransitionProbability > 1.0:
            raise ValueError( "emTransitionProbability > 1.0: " + str(self.emTransitionProbability) )
        if self.gammaEmissionProbability > 1.0:
            raise ValueError( "gammaEmissionProbability > 1.0: " + str(self.gammaEmissionProbability) )
        if self.internalConversionProbability > 1.0:
            raise ValueError( "internalConversionProbability > 1.0: " + str(self.internalConversionProbability) )

    def setMultipolarity(self, mode, spinChange):
        """
        sets the multipolarity of a transition
        """
        if mode not in ["E","M"]:
            raise ValueError("Mode must be 'E' or 'M'")
        if not isinstance(spinChange, int):
            raise TypeError("spinChange must be an int")
        self.multipolarity=(mode, spinChange)

    def levelReport(self):
        # multipolarity string
        if self.multipolarity[0] is None: multipolarityString='?'
        else: multipolarityString=self.multipolarity[0]
        if self.multipolarity[1] is None: multipolarityString+='?'
        else:multipolarityString+=str(self.multipolarity[1])

        # gamma energy
        try: energy=self.getEnergy('MeV')
        except: energy=self.energy.value

        # halflife
        try: halflifeString=self.halflife.getValueAs('ns')
        except: halflifeString=str(self.halflife)

        result = 50*' ' + ( ''.join( map(
            lambda x: str(x).rjust(20),
            (energy, self.finalLevel, '%15g'%self.gammaEmissionProbability, self.nonRadiativeProbability, multipolarityString, halflifeString) ) ) )
        return result

    def toRIPLFormat(self):
        '''
        Two typical examples of the gamma records are given below (Nb-94)::

                  Nf    Eg[MeV]       Pg          Pe          ICC
                  3      0.055     4.267E02    1.301E01    2.050E+00
                  1      0.113     7.499E01    8.699E01    8.699E01

                Nf    : sequential number of the final state
                Eg    : gamma-ray energy in MeV
                Pg    : Probability that a level decays through the given gamma-ray emission.
                        Pg is the ratio of the total electromagnetic decay of the level to the
                        intensity of the gamma-ray. If no branching ratio is given in the
                        ENSDF file, Pg=0.
                Pe    : Probability that a level decays with the given electromagnetic
                        transition, i.e., ratio of the total electromagnetic decay of the level
                        to the intensity of the given electromagnetic transition. The sum of
                        electromagnetic decays is normalized to 1. If no branching ratio is
                        given in the ENSDF file, Pe=0.
                ICC   : Internal conversion coefficient of a transition.

        The corresponding FORTRAN format is (39x,i4,1x,f10.3,3(1x,e10.3))
        '''
        raise NotImplementedError()

    def toTALYSFormat(self):
        '''
        The number of the level to which the gamma-ray decay takes place and the corresponding
        branching ratio and electron-conversion factor. If a branching ratio has been assigned by
        us, a 'B' is placed in column 58. The FORTRAN format for this line is
            (29x,i3,f10.6,e10.3,5x,a1)
        '''
        raise NotImplementedError()


# ---------------------------------------------------------------------------------
#
#    Level scheme readers for different formats
#
# ---------------------------------------------------------------------------------
class ParseError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def readTALYSLevelScheme(f, verbose=False):
    '''
    Read a level scheme from a string "f" (assumed to have been read from a file)

    TALYS inputs have this structure:

  40  79    1    0                                                          79Zr
   0   0.000000   0.5    1  0                   5.600E-02 JP             HFB
  40  80    5    2                                                          80Zr
   0   0.000000   0.0    1  0                                               0+
   1   0.289900   2.0    1  1                                             (2+)
                               0  1.000000 2.441E-02
   2   0.827900   4.0    1  1                                             (4+)
                               1  1.000000 3.300E-03
    '''
    contentMap = xParticleList()
    lines = f.split('\n')
    while lines:
        line = lines.pop(0)

        if line.strip() == "": continue

        # Read in the nucleus
        name = line[75:81].strip()  # e.g. 92Zr
        A = int(line[5:9])  # e.g. 92
        Z = int(line[1:5])  # e.g. 40
        myNucleus = nucleus(name=nucleusNameFromZA(Z * 1000 + A))
        numLevels = int(line[14:19]) + 1  # e.g. 136
        if verbose:
            try:
                Sn = myNucleus.Sn.inUnitsOf('MeV/c**2')
            except:
                Sn = None
            try:
                Sp = myNucleus.Sp.inUnitsOf('MeV/c**2')
            except:
                Sp = None
            print( "Nucleus:", name, \
                "  ZA =", myNucleus.ZA, \
                "  mass =", myNucleus.mass, \
                "  Sn =", Sn, \
                "  Sp =", Sp, \
                '  numLevels =', numLevels)

        # Save the nucleus
        contentMap.addParticle(myNucleus)

        # Read in the level scheme
        for iLevel in range(numLevels):
            line = lines.pop(0)

            # Read in the level energy & index
            levelIndex = int(line[0:5])
            eLevel = PhysicalQuantity(float(line[7:16]), "MeV")

            # JPi assignment
            J = float(line[17:22])
            if J == -1: J = None
            Pi = int(line[22:27])
            if Pi == 0: Pi = None
            howJPiDetermined = line[59:65].strip()
            if howJPiDetermined in ['JP', 'J', 'P']: howJPiDetermined = 'chosen'
            alternateJPiAssignments = line[65:].strip()

            # Lifetime
            try:
                lifetime = PhysicalQuantity(float(line[41:58]), 's')
            except ValueError:
                lifetime = PhysicalQuantity(0.0, 's')

            # Decay information (we'll keep only the gamma data for now)
            numGammas = int(line[27:30])

            myLevel = nuclearLevel(name=myNucleus.name + '_e' + str(levelIndex), spin=J, parity=Pi,
                                   alternateJPiAssignments=alternateJPiAssignments, index=levelIndex, energy=eLevel)
            if verbose:
                print( "   ", myLevel, ':', \
                    '  energy =', myLevel.energy, \
                    '  JPi =', myLevel.spin, myLevel.parity, \
                    '  lifetime =', myLevel.lifetime, \
                    '  numGammas =', numGammas)

            # Save the level in the contentsMap and inform the nucleus about the level
            myNucleus.addLevel(myLevel)

            # Read in the gammas
            for iGamma in range(numGammas):
                line = lines.pop(0)

                # Read in the gamma
                indexFinal = int(line[25:33])
                finalLevel = myNucleus.name + '_e' + str(indexFinal)
                eGamma = myLevel.energy - myNucleus[indexFinal].energy
                Pe = float(line[33:43])
                ICC = float(line[43:53])
                howBRDetermined = line[53:63].strip()
                if howBRDetermined in ['B']: howJPiDetermined = 'chosen'

                myGamma = nuclearLevelGamma(energy=eGamma, finalLevel=finalLevel, Pg=Pe / (1.0 + ICC), ICC=ICC)
                if verbose:
                    print( "       ", myLevel.name + "->" + myGamma.finalLevel, ':', \
                        '  energy =', myGamma.energy, \
                        '  Pg =', myGamma.gammaEmissionProbability, \
                        '  ICC =', myGamma.internalConversionCoefficient)

                # Save the gamma
                myLevel.addGamma(myGamma)


            # Check if number of gammas found meet our expectations
            if numGammas != myLevel.getNumberGammas():
                raise ValueError("number of gammas stated in nucleus definition != number of gammas found: " + str(
                    numGammas) + " != " + str(myNucleus.getNumberGammas()))

        # Check if number of levels found meet our expectations
        if numLevels != myNucleus.getNumberLevels():
            raise ValueError("number of levels stated in nucleus definition != number of levels found: " + str(
                numLevels) + " != " + str(myNucleus.getNumberLevels()))

    if verbose:
        print( "Processed the following nuclei:", ', '.join(contentMap.keys()))

    return contentMap


def readEMPIRELevelScheme(f, verbose=False, updateMasses=False): return readRIPLLevelScheme(f, verbose=verbose,
                                                                                            updateMasses=updateMasses)


def readRIPLLevelScheme(f, verbose=False, updateMasses=False):
    '''
    Read a level scheme from a string "f" (assumed to have been read from a file)

    EMPIRE and RIPL have the same format (according to Mike Herman).  They have the following structure:

 92Zr   92   40  136  116   55    6    8.634800    9.397201
  1   0.000000   0.0  1  -1.00E+00  0 u                      0+  0
  2   0.934470   2.0  1   5.00E-12  1 u                      2+  0
                                          1      0.934  9.992E-01  1.000E+00  7.693E-04
  3   1.382810   0.0  1   8.80E-11  2 u                      0+  0
                                          2      0.448  9.942E-01  1.000E+00  5.842E-03
                                          1      1.383  0.000E+00  0.000E+00  0.000E+00
  4   1.495460   4.0  1   1.02E-10  1 u                      4+  0
                                          2      0.561  9.971E-01  1.000E+00  2.922E-03
    '''
    contentMap = collections.OrderedDict()
    lines = f.split('\n')
    numLines = len(lines)
    try:
        while lines:
            line = lines.pop(0)

            if line.strip() == "": continue

            # Read in the nucleus
            name = line[0:5].strip()  # e.g. 92Zr
            A = int(line[6:10])  # e.g. 92
            Z = int(line[11:16])  # e.g. 40
            myNucleus = nucleus(name=nucleusNameFromZA(Z * 1000 + A))
            if not updateMasses:
                try:
                    myNucleus.Sn = PhysicalQuantity(float(line[37:48]), 'MeV/c**2') * cSquared  # e.g. 8.634800 MeV/c**2
                except:
                    myNucleus.Sn = None
                try:
                    myNucleus.Sp = PhysicalQuantity(float(line[48:60]), 'MeV/c**2') * cSquared  # e.g. 9.397201 MeV/c**2
                except:
                    myNucleus.Sp = None
            myNucleus.lastLevelInCompleteScheme = int(line[26:31])  # e.g. 55
            myNucleus.lastLevelWithAssignedJPiIndex = int(line[31:36])  # e.g. 6
            numLevels = int(line[16:21])  # e.g. 136
            totNumGammas = int(line[21:26])  # e.g. 116
            if verbose:
                print( "Nucleus:", name, \
                    "  ZA =", myNucleus.ZA, \
                    "  mass =", myNucleus.mass, \
                    "  Sn =", myNucleus.Sn, \
                    "  Sp =", myNucleus.Sp, \
                    '  numLevels =', numLevels)

            # Save the nucleus
            contentMap[str(myNucleus)] = myNucleus

            # Read in the level scheme
            for iLevel in range(numLevels):
                line = lines.pop(0)

                # Read in the level energy & index
                levelIndex = int(line[0:4]) - 1
                eLevel = PhysicalQuantity(float(line[4:15]), "MeV")

                # JPi assignment
                J = float(line[15:21])
                if J == -1:
                    J = None
                else:
                    J = Spin(J)
                Pi = int(line[21:24])
                if Pi == 0:
                    Pi = None
                else:
                    Pi=Parity(Pi)
                howJPiDetermined = line[38:40].strip()
                alternateJPiAssignments = line[40:64].strip()

                # Lifetime
                try:
                    lifetime = PhysicalQuantity(float(line[24:35]), 's')
                except ValueError:
                    lifetime = PhysicalQuantity(0.0, 's')

                # Decay information (we'll keep only the gamma data for now)
                numGammas = int(line[35:38])
                numDecays = int(line[64:67])
                decayModeString = line[67:].strip()

                myLevel = nuclearLevel(name=myNucleus.name + '_e' + str(levelIndex), spin=J, parity=Pi,
                                       howJPiDetermined=howJPiDetermined, alternateJPiAssignments=alternateJPiAssignments,
                                       index=levelIndex, energy=eLevel, lifetime=lifetime, groundState=myNucleus)
                if verbose:
                    print( "   ", myLevel, ':', \
                        '  energy =', myLevel.energy, \
                        '  JPi =', myLevel.spin, myLevel.parity, \
                        '  lifetime =', myLevel.lifetime, \
                        '  numGammas =', numGammas)

                # Read in the gammas
                for iGamma in range(numGammas):
                    line = lines.pop(0)

                    # Read in the gamma
                    indexFinal = int(line[33:44]) - 1
                    eGamma = PhysicalQuantity(float(line[44:55]), "MeV")
                    Pg = float(line[55:66])
                    Pe = float(line[66:77])
                    ICC = float(line[77:88])

                    myGamma = nuclearLevelGamma(energy=eGamma, startLevel=myLevel, finalLevel=myNucleus.levels[indexFinal], Pg=Pg, ICC=ICC, Pe=Pe)

                    if verbose:
                        print( "       ", str(myLevel) + "->" + str(myGamma.finalLevel), ':', \
                            '  energy =', myGamma.energy, \
                            '  Pg =', myGamma.gammaEmissionProbability, \
                            '  ICC =', myGamma.internalConversionCoefficient)

                    # Save the gamma
                    myLevel.addGamma(myGamma)

                # Save the level in the contentsMap and inform the nucleus about the level
                myNucleus.addLevel(myLevel)

                # Check if number of gammas found meet our expectations
                if numGammas != myLevel.getNumberGammas():
                    raise ValueError("number of gammas stated in nucleus definition != number of gammas found: " + str(
                        numGammas) + " != " + str(myNucleus.getNumberGammas()))

            # Check if number of levels found meet our expectations
            if numLevels != myNucleus.getNumberLevels():
                raise ValueError("number of levels stated in nucleus definition != number of levels found: " + str(
                    numLevels) + " != " + str(myNucleus.getNumberLevels()))
    except KeyboardInterrupt: exit()
    except Exception as err:
        badLineNumber = numLines-len(lines)
        print( "\nException raised while parsing levels:\n%s"% '\n'.join(f.split('\n')[max(0,badLineNumber-5):min(numLines,badLineNumber+5)]))
        print( '-'*60)
        traceback.print_exc(file=sys.stdout)
        print( '-'*60)
        raise ParseError("Line %i caused an error %s"% (badLineNumber, str(err)))

    if verbose:
        print ("Processed the following nuclei:", ', '.join(contentMap.keys()))

    return contentMap



# ---------------------------------------------------------------------------------
#
#      A widget to make a "complete" nucleus from RIPL data
#
# ---------------------------------------------------------------------------------


def getNucleus(Z, A, pathToRIPLLevelFiles=None, pathToRIPLResonanceFiles=None, pathToRIPLLevelDensityFiles=None, verbose=False, updateMasses=True):

    RIPLFile = pathToRIPLLevelFiles + os.sep + 'z' + str(Z).zfill(3) + '.dat'

    # Read in the the entire Z chain
    nuclMap = readRIPLLevelScheme(
        ''.join(open(RIPLFile, mode='r').readlines()),
        verbose=verbose, updateMasses=updateMasses)

    # Pick out the nucleus requested
    myNucleus=None
    for nucl in nuclMap:
        if nuclMap[nucl].A==A:
            myNucleus=nuclMap[nucl]
            break
    if myNucleus is None:
        raise RuntimeError( "Z=%i, A=%i not found in file %s"%(Z, A, RIPLFile))

    # Read in the the entire Z chain
    swaveMap = readRIPLMeanLevelSpacing(
        ''.join(open(pathToRIPLResonanceFiles + os.sep + 'resonances0.dat', mode='r').readlines()), verbose=verbose)
    pwaveMap = readRIPLMeanLevelSpacing(
        ''.join(open(pathToRIPLResonanceFiles + os.sep + 'resonances1.dat', mode='r').readlines()), verbose=verbose)

    # Update in the mean level spacings
    # The RIPL resonance files contain data for n + targ.
    # If you want a particular ZA, you have to subtract off the "n" from the ZA to get the actual "targ"
    cnName = myNucleus.name.replace(str(A), str(A - 1))
    try:
        myNucleus.resSpacingSWave = swaveMap.get(cnName)[0]
    except:
        pass
    try:
        myNucleus.uncResSpacingSWave = swaveMap.get(cnName)[1]
    except:
        pass
    try:
        myNucleus.resSpacingPWave = pwaveMap.get(cnName)[0]
    except:
        pass
    try:
        myNucleus.uncResSpacingPWave = pwaveMap.get(cnName)[1]
    except:
        pass

    # Add in the level density.  Assume only HFBM for now.
    hfbFilePrefix = os.sep.join([pathToRIPLLevelDensityFiles, 'total', 'level-densities-hfb', 'z' + str(Z).zfill(3)])
    if not os.path.exists(hfbFilePrefix + '.tab'):
        myNucleus.levelDensity = None
    else:
        try:
            myNucleus.levelDensity = readHFBMLevelDensityTable(''.join(open(hfbFilePrefix + '.tab', mode='r').readlines()),
                                                       myNucleus.Z, myNucleus.A, verbose=verbose)
            if os.path.exists(hfbFilePrefix + '.cor'):
                corTable = readHFBMLevelDensityCorrections(''.join(open(hfbFilePrefix + '.cor', mode='r').readlines()),verbose=verbose)
            else: corTable={}
            try:
                myNucleus.levelDensity.aTildeCorrection = corTable[( myNucleus.Z, myNucleus.A )][0]
            except KeyError:
                myNucleus.levelDensity.aTildeCorrection = 0.0
            try:
                myNucleus.levelDensity.pairingShiftCorrection = corTable[( myNucleus.Z, myNucleus.A )][1]
            except KeyError:
                myNucleus.levelDensity.pairingShiftCorrection = 0.0
        except KeyError:
            myNucleus.levelDensity = None
    return myNucleus


# ---------------------------------------------------------------------------------
#
#    Main !!
#
# ---------------------------------------------------------------------------------


if __name__ == "__main__":
    import unittest

    # ---------------------------------------------------------------------------------
    #
    #   Unit tests
    #
    # ---------------------------------------------------------------------------------

    class TestNucleus(unittest.TestCase):

        def setUp(self):
            self.a = nucleus(name="C12", mass=PhysicalQuantity(12.0, "amu"), attributes={})  # should be same as next version
            self.a = nucleus(name="C12")

        def test_init(self):
            self.assertEqual(str(self.a), "C12")
            self.assertEqual(self.a.Z, 6)
            self.assertEqual(self.a.A, 12)
            self.assertEqual(self.a.ZA, 6012)
            self.assertEqual(self.a.symbol, 'C')
            self.assertEqual(self.a.mass,
                             PhysicalQuantity(12.0, 'amu'))  # better be true, it's the definition of an amu!
            self.assertAlmostEqual(self.a.Sn.getValueAs('keV'), 18721.66,
                                   delta=0.1)  # neutron separation energy must agree w/ Audi-Wapstra 2003 to less than 0.1 eV
            self.assertAlmostEqual(self.a.Sp.getValueAs('keV'), 15956.90,
                                   delta=0.1)  # proton separation energy must agree w/ Audi-Wapstra 2003 to less than 0.1 eV

        def test_addLevel(self): pass

        def test_getNumberLevels(self): pass

        def test_getNumberGammas(self): pass

        def test_getLastLevelInCompleteSchemeIndex(self): pass

        def test_getLastLevelWithAssignedJPiIndex(self): pass

        def test_toRIPLFormat(self): pass

        def test_toTALYSFormat(self): pass


    class TestNuclearLevel(unittest.TestCase):

        def setUp(self):
            self.a = nuclearLevel(name="C12_e1", energy=PhysicalQuantity(4438.91, "keV"), spin=1, index=1, gammas=[],
                                  attributes={}, groundState=None)
            self.b = nuclearLevel(name="C12_e1", energy=PhysicalQuantity(4438.91, "keV"), spin=2.5, parity=-1, index=1,
                                  howJPiDetermined='c', alternateJPiAssignments="(3/2,5/2,7/2)-", gammas=[],
                                  attributes={}, groundState=None)
            self.c = nuclearLevel(name="C12_e1", energy=PhysicalQuantity(4438.91, "keV"), index=1, gammas=[],
                                  attributes={}, groundState=None)

        def test_init(self):
            print(self.a)

        def test_init0(self):
            self.assertEqual(self.a.spin, 1.0)
            self.assertEqual(self.a.attributes['spin'], 1.0)
            self.assertEqual(self.a.parity, None)
            self.assertEqual(self.a.attributes['parity'], None)
            self.assertEqual(self.a.howJPiDetermined, 'unknown')
            self.assertEqual(self.a.alternateJPiAssignments, None)

        def test_init1(self):
            self.assertEqual(self.b.spin, 2.5)
            self.assertEqual(self.b.attributes['spin'], 2.5)
            self.assertEqual(self.b.parity, -1)
            self.assertEqual(self.b.attributes['parity'], -1)
            self.assertEqual(self.b.howJPiDetermined, "chosen")
            self.assertEqual(self.b.alternateJPiAssignments, "(3/2,5/2,7/2)-")

        def test_getSpin(self): pass

        def test_getParity(self): pass

        def test_renormalizeTransitionProbabilities(self): pass

        def test_getNumberGammas(self): pass

        def test_getNumberDecayModes(self): pass

        def test_toRIPLFormat(self): pass

        def test_toTALYSFormat(self): pass

        def test_getName(self): pass

        def test_getToken(self): pass

        def test_getMass(self): pass

        def test_getLevelIndex(self): pass

        def test_getMetaStableIndex(self): pass

        def test_getLevelAsFloat(self): pass

        def test_getZ_A_SuffixAndZA(self): pass

        def test_getJpi(self): pass

        def test_toXMLList(self): pass

        def test_toENDF6(self): pass


    class TestNuclearDecayMode(unittest.TestCase):

        def setUp(self):
            self.a = nuclearDecayMode(mode=None, probability=0.0)

        def test_init(self):
            print (self.a)

        def test_toRIPLFormat(self): pass

        def test_toTALYSFormat(self): pass


    class TestNuclearLevelGamma(unittest.TestCase):

        def setUp(self):
            self.a = nuclearLevelGamma(energy=PhysicalQuantity("1.322 MeV"), finalLevel="Ni58_e1", Pg=0.96, ICC=0.01)
            self.b = nuclearLevelGamma(energy=PhysicalQuantity("1.322 MeV"), finalLevel="Ni58_e1", probability=0.96,
                                       nonRadiativeProbability=0.0)

        def test_init_a(self):
            self.assertAlmostEqual(self.a.gammaEmissionProbability, 0.96)
            self.assertAlmostEqual(self.a.internalConversionProbability, 0.0096)
            self.assertAlmostEqual(self.a.internalConversionCoefficient, 0.01)
            self.assertAlmostEqual(self.a.emTransitionProbability, 0.9696)

        def test_init_b(self):
            self.assertAlmostEqual(self.b.gammaEmissionProbability, 0.96)
            self.assertAlmostEqual(self.b.internalConversionProbability, 0.0)
            self.assertAlmostEqual(self.b.internalConversionCoefficient, 0.0)
            self.assertAlmostEqual(self.b.emTransitionProbability, 0.96)

        def test_toRIPLFormat(self): pass

        def test_toTALYSFormat(self): pass

        def test_toENDF6List(self): pass

        def test_toXMLList(self): pass


    class TestMisc( unittest.TestCase ):

        def test_double_factorial(self):
            self.assertEqual( double_factorial(0), 1 )
            self.assertEqual( double_factorial(1), 1 )
            self.assertEqual( double_factorial(2), 2 )
            self.assertEqual( double_factorial(3), 3 )
            self.assertEqual( double_factorial(4), 8 )
            self.assertEqual( double_factorial(9), 945 )
            self.assertEqual( double_factorial(10), 3840 )
            self.assertEqual( double_factorial(11), 10395 )

        def test_is_transition_allowed( self ):
            self.assertEqual(is_transition_allowed( ("E",1), 2.0, +1, 3.0, -1 ), True  )
            self.assertEqual(is_transition_allowed( ("E",1), 2.0, +1, 3.0, +1 ), False )
            self.assertEqual(is_transition_allowed( ("E",1), 7.0, +1, 3.0, -1 ), False )
            self.assertEqual(is_transition_allowed( ("E",2), 2.0, +1, 0.0, +1 ), True  )
            self.assertEqual(is_transition_allowed( ("E",2), 2.0, +1, 0.0, -1 ), False )
            self.assertEqual(is_transition_allowed( ("E",2), 3.0, +1, 0.0, +1 ), False )
            self.assertEqual(is_transition_allowed( ("M",1), 2.5, +1, 1.5, +1 ), True  )
            self.assertEqual(is_transition_allowed( ("M",1), 2.5, +1, 1.5, -1 ), False )
            self.assertEqual(is_transition_allowed( ("M",1), 3.5, +1, 1.5, +1 ), False )
            self.assertEqual(is_transition_allowed( ("M",2), 2.0, -1, 0.0, +1 ), True  )
            self.assertEqual(is_transition_allowed( ("E",5), 5.5, -1, 0.5, +1 ), True  )

        def test_WeisskopfSingleParticleEstimateBXL(self):
            """
            Answers from NS Memo 1B/1 (82)
            """
            self.assertAlmostEqual(WeisskopfSingleParticleEstimateBXL(('E',1),86),6.446e-4*pow(86.,2./3.),5)
            self.assertAlmostEqual(WeisskopfSingleParticleEstimateBXL(('E',2),86),5.940e-6*pow(86.,4./3.),5)
            self.assertAlmostEqual(WeisskopfSingleParticleEstimateBXL(('M',1),86),1.791,2)
            self.assertAlmostEqual(WeisskopfSingleParticleEstimateBXL(('E',3),86),5.940e-8*86.0*86.0,4)
            self.assertAlmostEqual(WeisskopfSingleParticleEstimateBXL(('M',2),86),1.650e-2*pow(86.,2./3.),4)

        def test_gammaHalflife(self):
            """
            Liz checked these values!
            """
            self.assertAlmostEqual(gammaHalflife(('E',1),1.0,PhysicalQuantity(1.0,'MeV')).getValueAs('ns'),4.3587980198378215e-09)
            self.assertAlmostEqual(gammaHalflife(('M',2),1.0,PhysicalQuantity(1.0,'MeV')).getValueAs('ns'),0.5116100181270946)


    # Actually run the tests
    try:
        import xmlrunner

        unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-results'))
    except ImportError:
        unittest.main()
        print()
