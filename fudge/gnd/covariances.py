# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Module with containers for covariances in several different forms
"""

import sys
import math
import fudge
from fudge.gnd import link
from fudge.core.ancestry import ancestry
from fudge.core.math import matrix as gndMatrix
from fudge.core.math.xData import XYs, axes
from fudge.core.utilities import fudgeExceptions
from fudge.legacy.converting import endfFormats
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

__metaclass__ = type

""" possible forms for the covariance (section.nativeData should list one of the following): """
covarianceFormToken = 'covarianceMatrix'
mixedFormToken = 'mixed'
summedFormToken = 'sum'
energyIntervalFormToken = 'perEnergyInterval'
legendreOrderCovarianceFormToken = 'LegendreOrderCovariance'
legendreLValueFormToken = 'LegendreLValue'
relativeToken = 'relative'
absoluteToken = 'absolute'


class covarianceSuite( ancestry ):
    """
    All covariances for a target/projectile combination are stored in a :py:class:`covarianceSuite` in gnd.
    The :py:class:`covarianceSuite` is stored in a separate file from the reactionSuite.

    Within the :py:class:`covarianceSuite`, data is sorted into :py:class:`sections` (see the section class below), each
    of which contains one section of the full covariance matrix.
    
    """
    def __init__(self, projectile=None, target=None, reactionSums=None,
            externalReactions=None, sections=None, modelParameterCovariances=None, styles=None):

        ancestry.__init__( self, 'covarianceSuite', None )

        self.projectile = projectile            #: The projectile
        self.target = target                    #: The target
        self.reactionSums = reactionSums or []  #: List of lumped sums, etc
        self.externalReactions = externalReactions or [] #: List of cross-material covariances
        self.sections = sections or []          #: List of section instances
        self.modelParameterCovariances = modelParameterCovariances or [] #: List of modelParameterCovariance instances
        self.styles = styles or {}              #: Python dict of styles
        self.format = "gnd version 1.0"         #: For now, it is "gnd version 1.0"
    
    def __getitem__(self, idx):
        return (self.sections+self.modelParameterCovariances)[idx]
    
    def __len__(self):
        return len(self.sections+self.modelParameterCovariances)

    def addSection(self,section):
        section.setParent( self )
        self.sections.append(section)

    def addModelParameterCovariance(self,mpcovar):
        mpcovar.setParent( self )
        self.modelParameterCovariances.append(mpcovar)
    
    def addReactionSum(self, reactionSum):
        reactionSum.setParent( self )
        self.reactionSums.append(reactionSum)

    def addExternalReaction(self, externalReaction):
        externalReaction.setParent( self )
        self.externalReactions.append(externalReaction)

    def saveToOpenedFile( self, fOut, flags=None, verbosityIndent='' ):
        xmlString = self.toXMLList( flags=flags, verbosityIndent=verbosityIndent )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, flags=None, verbosityIndent='' ):
        with open(fileName,"w") as fout:
            self.saveToOpenedFile( fout, flags=flags, verbosityIndent=verbosityIndent )
    
    def check( self, **kwargs ):
        """
        Check all covariance sections, returning a list of warnings. 
        
        :keyword bool checkUncLimits: Should we check the uncertainty limits? (default: True)
        :keyword float minRelUnc: Minimum allowable relative uncertainty (default: 0.0)
        :keyword float maxRelUnc: Maximum allowable relative uncertainty (default: 10.0)
        :keyword theData: A reference to the data for this covariance.  This is useful for converting between relative and absolute covariance (default: None)
        :keyword float negativeEigenTolerance: Ignore negative eigenvalues smaller than this (default: -1e-6) 
        :keyword float eigenvalueRatioTolerance: Warn if smallest eigenvalue < this value * biggest (default: 1e-8)
        :keyword float eigenvalueAbsoluteTolerance: Warn if smallest eigenvalue < this value (default: 1e-14)
        :rtype: warning.context
        """

        from fudge.gnd import warning

        # default input options
        options = {
            'checkUncLimits': True,
            'minRelUnc': 0.0,
            'maxRelUnc': 10.0,
            'theData': None,
            'negativeEigenTolerance': -1e-6,  # ignore smaller negative eigenvalues
            'eigenvalueRatioTolerance': 1e-8, # warn if smallest eival < 1e-8 * biggest
            'eigenvalueAbsoluteTolerance': 1e-14,
        }
        for key in kwargs:
            if key in options: options[key] = kwargs[key]
            else: raise KeyError, "check() received unknown keyword argument '%s'" % key

        # assemble some useful info, to be handed down to children's check() functions:
        info = { 'covarianceSuite': self }
        info.update( options )

        warnings = []

        #: summedReactions are sums of other sections, but that opens the possibility of a cyclic dependency
        #: check for that, using algorithm taken from
        #:    http://neopythonic.blogspot.com/2009/01/detecting-cycles-in-directed-graph.html 
        def find_cycle(NODES, EDGES):
            todo = set(NODES)
            while todo:
                node = todo.pop()
                stack = [node]
                while stack:
                    top = stack[-1]
                    for node in EDGES(top):
                        if node in stack:
                            return stack[stack.index(node):]
                        if node in todo:
                            stack.append(node)
                            todo.remove(node)
                            break
                    else:  # this node is not in a cycle
                        node = stack.pop()
            return None

        def get_edges( section ):   #: return list of all pointers from this section
            natDat = section.getNativeData()
            if isinstance(natDat, summedCovariance): return [v.link for v in natDat.pointerList]
            elif isinstance(natDat, mixedForm):
                edges = []
                for part in natDat:
                    if isinstance(part, summedCovariance): edges += [v.link for v in part.pointerList]
                return edges
        nodes = [sec for sec in self.sections if get_edges( sec )]  # sections that contain pointers

        if nodes:
            cycle = find_cycle(nodes, get_edges)
            if cycle:
                warnings.append( warning.cyclicDependency( cycle ) )

        # check each section
        for section in self.sections:
            sectionWarnings = section.check( info )
            if sectionWarnings:
                warnings.append( warning.context('Section: %s' % section.id, sectionWarnings ) )

        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)
    
    def fix( self, **kwargs ):
        '''
        Apply basic fixes to a covariance
        
        :param bool removeNegativeEVs: Should we remove eigenspaces corresponding to negative eigenvalues? (Default: True) 
        :param bool removeSmallEVs: Should we remove eigenspaces corresponding to small eigenvalues? (Default: False) 
        :param bool fixUncLimits: Should we fix the uncertainties to lie within bounds imposed by minRelUnc and maxRelUnc? (Default: False) 
        :param minRelUnc: Minimum allowable uncertainty for this covariance
        :type minRelUnc: float or None 
        :param maxRelUnc: Maximum allowable uncertainty for this covariance
        :type maxRelUnc: float or None 
        :param theData: Reference to the data that accompanies this covariance so that we may convert between absolute and relative covariance as needed.
        :type theData: instance or None 

        '''
        from fudge.gnd import warning
        # default input options
        options = {
            'removeNegativeEVs':True, 
            'removeSmallEVs':False, 
            'fixUncLimits':False, 
            'minRelUnc':None, 
            'maxRelUnc':None, 
            'theData':None,
        }
        for key in kwargs:
            if key in options: options[key] = kwargs[key]
            else: raise KeyError, "fix() received unknown keyword argument '%s'" % key
        # assemble some useful info, to be handed down to children's check() functions:
        info = { 'covarianceSuite': self }
        info.update( options )
        # do the fixing
        warnings = []
        for section in self.sections: warnings += section.fix( **info )        
        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)
    
    def removeExtraZeros(self):
        """ remove columns/rows of zero from all covariances """
        for section in self.sections:
            for form in sections.forms.values():
                if hasattr(form, "removeExtraZeros"):
                    form.removeExtraZeros()
                if isinstance(form, mixedForm):
                    for covar in form.subsections:
                        if hasattr(covar, "removeExtraZeros"):
                            covar.removeExtraZeros()
    
    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        """Write self out to GND-XML"""
        indent2 = indent+'  '; indent3 = indent2+'  '
        if flags is None: flags = {}
        xmlString = ['<?xml version="1.0" encoding="UTF-8"?>']
        xmlString.append(indent+'<covarianceSuite projectile="%s" target="%s" format="%s"'
                % (self.projectile, self.target, self.format))
        xmlString[-1] += ' xmlns:xlink="http://www.w3.org/1999/xlink">'
        xmlString.append('%s<styles>' % (indent+'  '))
        for style in self.styles: xmlString += self.styles[style].toXMLList( indent=indent+'    ' )
        xmlString[-1] += '</styles>'
        if self.reactionSums:
            xmlString += [indent2+'<reactionSums>',
                indent3+'<!-- Covariances may be given for a sum of several reactions. Define these "summed reactions" here: -->']
            for reactionSum in self.reactionSums:
                xmlString += reactionSum.toXMLList(flags,verbosityIndent,indent=indent3)
            xmlString[-1] += '</reactionSums>'
        if self.externalReactions:
            xmlString += [indent2+'<externalReactions>',
                indent3+"<!-- This target has covariances with reactions on different target(s). List other target/reactions: -->"]
            for externalReaction in self.externalReactions:
                xmlString += externalReaction.toXMLList(flags,verbosityIndent,indent=indent3)
            xmlString[-1] += '</externalReactions>'
        for section in self.sections + self.modelParameterCovariances:
            xmlString += section.toXMLList(flags,verbosityIndent,indent+'  ')
        xmlString.append(indent+'</covarianceSuite>')
        return xmlString

    @staticmethod
    def parseXMLNode( CSelement, linkData={} ):
        '''Parse an XML node and load create a covarianceSuite based on its contents'''
        covariances = covarianceSuite( CSelement.get('projectile'), CSelement.get('target') )
        for child in CSelement:
            if child.tag == 'styles':
                styles = [fudge.gnd.miscellaneous.style.parseXMLNode( style ) for style in child]
                for style in styles: covariances.styles[ style.name ] = style
            elif child.tag == 'reactionSums':
                for reactionSumElement in child:
                    covariances.addReactionSum( reactionSum.parseXMLNode( reactionSumElement, linkData ) )
            elif child.tag == 'externalReactions':
                for external in child:
                    covariances.addExternalReaction( externalReaction.parseXMLNode( external, linkData ) )
            elif child.tag == 'section':
                covariances.addSection( section.parseXMLNode( child, linkData ) )
            elif child.tag in ('resonanceParameterCovariance',):
                covariances.addModelParameterCovariance( resonanceParameterCovariance.parseXMLNode( child, linkData ) )
            else:
                raise Exception("unknown element '%s' found in the covarianceSuite!" % child.tag)

        # fix links:
        for link_ in linkData['unresolvedLinks']:
            path = link_.path
            if path.startswith('/reactionSuite'): root = linkData['reactionSuite']
            elif path.startswith('/covarianceSuite'): root = covariances
            if root is None: continue   # in case user didn't specify the reactionSuite
            res = link.follow(path, root)
            link_.link = res
        return covariances
    
    def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
        """ go back to ENDF format """
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        NIS, ABN = 1,1.0; ZAI=ZAM  # assuming one isotope/file
        MTL = 0 # mtl=1 sections are handled in lumpedCovariance

        for section in self.modelParameterCovariances:
            section.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)

        sections = self.sections[:]
        # sort covariances by MF/MT:
        mfmts = []
        for section in sections:
            mfmts.append(map(int,section.rowData['ENDF_MFMT'].split(',')))
        if len(mfmts)==0: return

        mfs, mts = zip(*mfmts)
        zipList = zip(mfs,mts,sections)
        idx = 0
        while idx<len(zipList):
            mf,mt,covar = zipList[idx]
            thisMFMT = [a[2] for a in zipList if a[:2]==(mf,mt)]
            idx += len(thisMFMT)
            
            if mf in (31,33):
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, MTL, 0, len(thisMFMT) )]
            elif mf==34:
                if thisMFMT[0].nativeData == legendreOrderCovarianceFormToken: LTT = 1
                NMT1 = len(thisMFMT)
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, LTT, 0, NMT1 )]
            elif mf==35:
                NK = 1
                if thisMFMT[0].nativeData == energyIntervalFormToken:
                    NK = len(thisMFMT[0].forms[energyIntervalFormToken])
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, MTL, NK, 0 )]
            elif mf==40:
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, 0, len(thisMFMT), 0 )]
            for section in thisMFMT:
                MAT1 = 0
                if section.columnData and isinstance(section.columnData.link, externalReaction):
                    otherTarget = section.columnData.link.target
                    from fudge.legacy.converting import endf_endl
                    ZA, MAT1 = endf_endl.ZAAndMATFromParticleName( otherTarget )
                if mf==34:
                    form = section.forms[legendreOrderCovarianceFormToken]
                    L1s = [subsec.L1 for subsec in form]
                    L2s = [subsec.L2 for subsec in form]
                    NL = len( set(L1s) );  NL1 = len( set(L2s) )
                    if section.columnData: raise NotImplemented # cross-reaction or cross-material
                    MT1 = mt
                    endf += [ endfFormats.endfHeadLine( 0.0, 0.0, MAT1, MT1, NL, NL1 ) ]
                if mf==40:
                    rowData = section.rowData
                    if type(rowData) is str: raise Exception("Don't string me along!")
                    else:
                        quant = rowData.link
                        from fudge import gnd
                        if isinstance(quant, gnd.covariances.reactionSum):
                            quant = quant.reactions[0].link
                        reaction = quant.parent
                        QI = reaction.getQ('eV')
                        level = reaction.product.getLevelAsFloat('eV')
                        QM = QI + level
                        LFS = reaction.product.getLevelIndex()
                        NL = 1
                        endf += [endfFormats.endfHeadLine( QM, QI, 0, LFS, 0, NL ) ]
                        XMF1, XLFS1, NC, NI = 10,LFS, 0,1
                        endf += [endfFormats.endfHeadLine( XMF1,XLFS1,MAT1,mt,NC,NI )]
                form = section.forms[section.nativeData]
                targetInfo['MAT1'] = MAT1
                targetInfo['dataPointer'] = [section.rowData,section.columnData]
                endf += form.toENDF6(flags, targetInfo)
                targetInfo.dict.pop('dataPointer')
                targetInfo.dict.pop('MAT1')
            endf.append( endfFormats.endfSENDLineNumber() )
            if mt not in endfMFList[mf]:
                endfMFList[mf][mt] = []
            endfMFList[mf][mt] += endf
        # also add ENDF-style pointers for lumped covariance data:
        for reactionSum in self.reactionSums:
            MT1 = reactionSum.ENDF_MFMT[1]
            if MT1 not in range(850,872): continue
            for link in reactionSum.reactions:
                mf,mt = map(int, link['ENDF_MFMT'].split(','))
                endfMFList[mf][mt] = [endfFormats.endfHeadLine(ZAM,AWT,0,MT1,0,0),
                        endfFormats.endfSENDLineNumber() ]

        return

class reactionSum( ancestry ):
    """ 
    A single covariance matrix is often given for a sum (or 'lump') of several reaction channels.
    Define the sum here, then in the covariance <section>, refer to this summed reaction   
    """
    def __init__(self, id=None, reactions=None, ENDF_MFMT=None):
        ancestry.__init__( self, 'reactionSum', None, attribute = 'id' )
        self.id = id #: an identifier str
        self.reactions = reactions or [] #: a list of fudge.gnd.link's that point to the reactions that are lumped together
        self.ENDF_MFMT = ENDF_MFMT #: the ENDF MF & MT values, a tuple of form (MF, MT)

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<reactionSum id="%s" ENDF_MFMT="%i,%i">'%(self.id,self.ENDF_MFMT[0],self.ENDF_MFMT[1])]
        for ch in self.reactions:
            xmlString.append( ch.toXML( indent+'  ' ) )
        xmlString[-1] += '</reactionSum>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        from fudge.gnd import link
        rsum = reactionSum( **dict(element.items()) )
        rsum.ENDF_MFMT = map(int, rsum.ENDF_MFMT.split(','))
        for child in element:
            link_ = link.parseXMLNode( child, linkData )
            linkData['unresolvedLinks'].append( link_ )
            rsum.reactions.append( link_ )
        return rsum

class externalReaction( ancestry ):
    """ 
    Covariance may relate this target with another material ('cross-material covariance'). In this case,
    specify the other material and reaction here 
    """
    def __init__(self, id=None, target=None, ENDF_MFMT=None):
        ancestry.__init__( self, 'externalReaction', None, attribute = 'id' )
        self.id = id #: an identifier str
        self.target=target #: the name of the target isotope/evaluation
        self.ENDF_MFMT = ENDF_MFMT #: the ENDF MF & MT values, a tuple of form (MF, MT) 

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<externalReaction id="%s" target="%s" ENDF_MFMT="%i,%i"/>'%(self.id,
            self.target, self.ENDF_MFMT[0], self.ENDF_MFMT[1])]
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        exReac = externalReaction( **dict(element.items()) )
        exReac.ENDF_MFMT = map(int, exReac.ENDF_MFMT.split(','))
        return exReac

class inputParameter( ancestry ):
    """
    For use within a modelParameterCovariance: the rows/columns of that matrix each correspond to one input
    model parameter. Parameters must each be inputParameter or multipleInputParameter instances.
    """
    def __init__(self, name, path, unit=None):
        ancestry.__init__( self, 'inputParameter', None )
        self.name = name #: just a Python str
        self.path = path #: a fudge.gnd.link 
        self.unit = unit #: an acceptable PhysicalQuantityWithUncertainty unit

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        if self.unit: unit = ' unit="%s"' % self.unit
        else: unit = ''
        return [indent+'<parameter name="%s"%s xlink:href="%s"/>' % (self.name, unit, self.path.toXLink())]

class modelParameterCovariance( ancestry ):
    """ 
    Express covariance between input parameters for a model 
    """
    def __init__(self, label=None, inputParameters=None, matrix=None, type=None, **kwargs):
        ancestry.__init__( self, 'modelParameterCovariance', None, attribute = 'label' )
        self.label = label #: a str
        self.inputParameters = inputParameters #: list of inputParameter instances to relate rows/columns to parameters 
        for param in self.inputParameters: param.setParent( self )
        self.matrix = matrix #: the actual fudge.core.math.matrix instance with the covariance
        #self.matrix.setParent( self )
        self.type = type  #: dunno, got to ask Caleb or Bret
        self.tag = 'modelParameterCovariance' #: usually 'modelParameterCovariance'
        self.attributes = kwargs #: a Python dict
        
    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        matrixWarnings = self.matrix.check( info )
        if matrixWarnings:
            warnings.append( warning.context("Model parameter covariances", matrixWarnings ) )
        return warnings
    
    def fix( self, **kw ): 
        '''assemble some useful info, to be handed down to children's fix() functions'''
        info = {}
        info['rowENDF_MFMT'] = None
        info['columnENDF_MFMT'] = None
        info.update( kw )
        return self.matrix.fix( **info )

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        indent2 = indent+'    '
        attrStr = ''.join( [' %s="%s"' % (key, self.attributes[key]) for key in self.attributes
            if bool(self.attributes[key]) ] )
        xmllist = [indent+( '<%s label="%s" type="%s" %s>' % (self.tag, self.label, self.type, attrStr) )]
        xmllist.extend( [indent+'  <inputParameters>',
            indent2+'<!-- Each row of this matrix corresponds to a model parameter. Parameters may be listed singly,',
            indent2+'  as for scattering radii, or in the case of resonance parameters they may be given all together.',
            indent2+'  In that case, rows of the matrix correspond to a loop over parameters for each resonance,',
            indent2+'  with resonances sorted by energy. -->'] )
        for inputParam in self.inputParameters: xmllist += inputParam.toXMLList(indent=indent+'    ')
        xmllist[-1] += '</inputParameters>'
        xmllist += self.matrix.toXMLList(flags=flags, verbosityIndent=verbosityIndent, indent=indent+'  ')
        xmllist[-1] += ('</%s>' % self.tag)
        return xmllist

class loopOverResonanceParameters( ancestry ):
    """ 
    For resonance region covariances, we need a compact way to express many model inputs.
    Simplest is to specify a loop over the resonances 
    """
    def __init__(self, nResonances, parametersPerResonance, path):
        ancestry.__init__( self, 'loopOverResonanceParameters', None )
        self.nResonances = nResonances #: an int, the number of resonances
        self.parametersPerResonance = parametersPerResonance #: dunno, got to ask Caleb or Bret
        self.path = path #: a fudge.gnd.link pointing to the accompanying reactionSuite's resonance parameters

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        return [ indent +
                '<loopOverResonanceParameters nResonances="%i" parametersPerResonance="%s" xlink:href="%s"/>'
            % (self.nResonances, self.parametersPerResonance, self.path.toXLink()) ]

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        ll = link.parseXMLNode( element )
        path = link.follow( ll.path, linkData['reactionSuite'] )
        return loopOverResonanceParameters( int(ll.attributes["nResonances"]),
                ll.attributes["parametersPerResonance"], path )

class resonanceParameterCovariance( modelParameterCovariance ):
    """
    In the resonance region, covariances are given between resonance parameters (energy and widths).
    Generally, the dimension of the matrix is 3*(number of resonances) for light targets, and 4*(nres)
    for heavy targets (where the fission width must be given).
    
    We also allow including the scattering radius in the covariance, although ENDF files currently only
    have room to list the uncertainty (variance) on the scattering radius. 
    """

    def __init__(self, label=None, inputParameters=None, matrix=None, type=None, **kwargs):
        modelParameterCovariance.__init__(self, label, inputParameters, matrix, type)
        self.tag = 'resonanceParameterCovariance' #: usually set to 'resonanceParameterCovariance'
        self.attributes = kwargs #: a Python dict

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        params = []
        for param in element[0]:
            if param.tag=='loopOverResonanceParameters': param=loopOverResonanceParameters.parseXMLNode(param, linkData)
            elif param.tag=='parameter': pass
            params.append(param)
        Matrix = gndMatrix.parseXMLNode( element[1] )
        return resonanceParameterCovariance( inputParameters=params, matrix=Matrix, **dict(element.items()) )

    def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
        """ go back to ENDF format """
        def swaprows( matrix, i1, i2, nrows ):
            # may need to rearrange parameters: ENDF often sorts first by L rather than by energy
            rows = matrix[i1:i1+nrows].copy()
            matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
            cols = matrix[:,i1:i1+nrows].copy()
            matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

        # need the resonance parameters as well as covariance matrix:
        res = targetInfo['reactionSuite'].resonances
        RPs = res.resolved.nativeData.resonanceParameters
        NRes = self.inputParameters[0].nResonances

        # MF32 header information:
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        NIS, ABN, ZAI = 1, 1.0, ZAM  # assuming only one isotope per file
        endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0 )]
        LFW = RPs.getColumn('fissionWidthA') is not None; NER=1
        endf.append( endfFormats.endfHeadLine( ZAI,ABN,0,LFW,NER,0 ) )
        EL,EH = res.resolved.lowerBound.getValueAs('eV'), res.resolved.upperBound.getValueAs('eV')
        LRU,NRO =1,0
        LRF = {'SingleLevel_BreitWigner':1, 'MultiLevel_BreitWigner':2, 'Reich_Moore':3}[
                res.resolved.nativeData.moniker ]
        NAPS = not res.resolved.nativeData.calculateChannelRadius
        endf.append( endfFormats.endfHeadLine( EL,EH,LRU,LRF,NRO,NAPS ) )
        SPI = targetInfo['spin']
        AP = res.resolved.nativeData.scatteringRadius.getValueAs('10*fm')
        LCOMP=1
        if 'LCOMP=0' in self.attributes.get('endfConversionFlags',''): LCOMP=0
        elif 'LCOMP=2' in self.attributes.get('endfConversionFlags',''): LCOMP=2

        sortByL = ("sortByL" in self.attributes.get('endfConversionFlags',''))
        Ls = RPs.getColumn('L')
        NLS = len(set(Ls))
        if LCOMP==2 or not sortByL: NLS = 0
        ISR = int( isinstance(self.inputParameters[0], inputParameter) and
                ('scatteringRadius' in self.inputParameters[0].name) )
        endf.append( endfFormats.endfHeadLine( SPI,AP,0,LCOMP,NLS,ISR ) )
        MLS = 0
        if ISR:
            MLS = 1 # currently don't handle energy-dependent DAP
            DAP = PhysicalQuantityWithUncertainty( self.matrix.data[0][0], self.inputParameters[0].unit ).getValueAs('10*fm')
            if LRF==3:
                endf.append( endfFormats.endfHeadLine( 0,0,0,0,MLS,1 ) )
            endf.append( endfFormats.endfDataLine( [DAP] ) )

        # MF32 repeats the resonance parameter information.
        # Extract that info from reactionSuite.resonances:
        table = [RPs.getColumn('L'), RPs.getColumn('energy',units='eV'), RPs.getColumn('J'),
                RPs.getColumn('totalWidth',units='eV') or [0]*NRes,
                RPs.getColumn('neutronWidth',units='eV'), RPs.getColumn('captureWidth',units='eV'),
                RPs.getColumn('fissionWidthA') or [0]*NRes,
                RPs.getColumn('fissionWidthB') or [0]*NRes]
        CS = RPs.getColumn('channelSpin')
        if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
            CS = [2*(cs-SPI) for cs in CS]
            Js = [v[0]*v[1] for v in zip(table[2],CS)]
            table[2] = Js
        table = zip(*table)
        matrix = self.matrix.data[MLS:,MLS:].copy()
        MPAR = len(matrix) / len(table)

        if sortByL:
            # reorder resonances, sorting first by L and second by energy:
            table.sort()

            elist1 = [(lis[1],lis[4],lis[5]) for lis in table]
            elist2 = zip( RPs.getColumn('energy',units='eV'),
                    RPs.getColumn('neutronWidth',units='eV'),
                    RPs.getColumn('captureWidth',units='eV') )

            for i in range(len(elist1)):
                i2 = elist2.index( elist1[i] )
                if i2!=i:
                    swaprows( matrix, MPAR*i, MPAR*elist2.index( elist1[i] ), MPAR )
                    val = elist2[i]
                    elist2[i] = elist2[i2]; elist2[i2] = val

        if LCOMP==0:
            for L in set( Ls ):
                NRS = Ls.count(L)
                endf.append( endfFormats.endfHeadLine( AWT, 0, L, 0, 18*NRS, NRS ) )
                for i in range(len(table)):
                    if table[i][0]!=L: break
                    endf.append( endfFormats.endfDataLine( table[i][1:7] ) )
                    block = matrix[MPAR*i:MPAR*(i+1), MPAR*i:MPAR*(i+1)]
                    lis = [block[0,0], block[1,1], block[2,1], block[2,2]]
                    if MPAR==4:
                        lis += [block[3,1],block[3,2],block[3,3],0,0,0,0,0]
                    else:
                        lis += [0,0,0,0,0,0,0,0]
                    endf += endfFormats.endfDataList( lis )


        if LCOMP==1:
            NSRS, NLRS = 1,0    # short-range correlations only
            endf.append( endfFormats.endfHeadLine( AWT, 0, 0, 0, NSRS, NLRS ) )
            MPAR = len( self.inputParameters[0].parametersPerResonance.split(',') )
            NRB = NRes
            NVS = (NRB*MPAR)*(NRB*MPAR+1)/2 # length of the upper diagonal matrix
            endf.append( endfFormats.endfHeadLine( 0,0, MPAR, 0, NVS+6*NRB, NRB ) )

            for res in table:
                if LRF in (1,2):
                    endf.append( endfFormats.endfDataLine( res[1:7] ) )
                elif LRF==3:
                    endf.append( endfFormats.endfDataLine( res[1:3] + res[4:8] ) )

            dataList = []
            for i in range(len(matrix)): dataList.extend( list( matrix[i][i:] ) )
            endf += endfFormats.endfDataList( dataList )

        elif LCOMP==2:
            import numpy
            QX, LRX = 0, 0  # haven't encountered any competitive widths yet
            endf.append( endfFormats.endfHeadLine( AWT,QX,0,LRX, 12*NRes, NRes ) )
            dat = matrix.diagonal()
            for i in range(len(table)):
                if LRF in (1,2):
                    endf += endfFormats.endfDataList( table[i][1:7] )
                    endf += endfFormats.endfDataList( [dat[3*i],0,0,dat[3*i+1],dat[3*i+2],0] )
                elif LRF==3:
                    endf += endfFormats.endfDataList( table[i][1:3] + table[i][4:8] )
                    endf += endfFormats.endfDataList( [dat[3*i],0,dat[3*i+1],dat[3*i+2],0,0] )

            # correlation matrix:
            NDIGIT = [a for a in self.attributes['endfConversionFlags'].split(',') if a.startswith('NDIGIT')]
            NDIGIT = int( NDIGIT[0][-1] )
            nints = 56 // (NDIGIT+1)    # how many numbers fit on each line?
            if NDIGIT==3: nints = 13    # special case
            rsd = numpy.sqrt( matrix.diagonal() )
            rsd[ rsd==0 ] = 1
            corr_mat = matrix / numpy.outer( rsd,rsd )
            corr_mat = numpy.rint( corr_mat * 10**NDIGIT )  # rint: round to nearest int
            # write lower-diagonal as sparse matrix using INTG format:
            endfCorrMat = []
            for i in range(len(corr_mat)):
                vals = corr_mat[i,:i]
                j = 0
                while j < i:
                    if vals[j]!=0:
                        endfCorrMat.append( endfFormats.writeEndfINTG(
                            i+1, j+1, list(vals[j:j+nints]), NDIGIT ) )
                        j += nints
                    else: j+=1
            NNN = NRes * 3
            NM = len(endfCorrMat)
            endf.append( endfFormats.endfHeadLine( 0,0, NDIGIT, NNN, NM, 0 ) )
            endf += endfCorrMat
        endf.append( endfFormats.endfSENDLineNumber() )
        endfMFList[32][151] = endf


class section( ancestry ):
    """
    A covarianceSuite contains sections, where each section represents either a self-covariance for one quantity,
    or a cross-covariance between two quantities

    More generally, the covarianceSuite can be thought of as a single covariance matrix with all covariance data
    for a target/projectile. It is broken into sections, where each section holds a chunk of the full matrix.

    Within each section, covariance data can take multiple forms: :py:class:`covarianceMatrix` is the most common,
    but 'summed', 'mixed' are also possible. 
    
    Valid values in the :py:attr:`forms` dictionary are:
        * mixedForm
        * energyIntervalForm
        * summedCovariance
        * LegendreOrderCovarianceForm
        * covarianceMatrix
    """

    def __init__(self, label=None, id=None, rowData=None, columnData=None, nativeData=None):
        """ each section needs a unique id, pointers to the central values (row and column), 
        and one or more forms """
        ancestry.__init__( self, 'section', None, attribute = 'label' )
        self.label = label #: a str label that gets used on plots, etc.
        self.id = id #: a str identifier, useful for resolving links
        self.rowData = rowData #: a fudge.gnd.link pointing to the corresponding data for the covariance row
        self.columnData = columnData #: a fudge.gnd.link pointing to the corresponding data for the covariance column
        self.nativeData = nativeData #: a str identifying which form is the one the evaluator intended.  This is the key to the forms Python dict.
        self.forms = {} #: a Python dict.  They key is just a string identifying the form.   The value is the actual covariance matrix instance.  

    def addForm( self, form ):
        form.setParent( self )
        self.forms[ form.moniker ] = form

    def getNativeData( self ):

        return self.forms[ self.nativeData ]

    def check( self, info ):
        """ check each section """

        from fudge.gnd import warning
        warnings = []
        for f in self.forms:
            formWarnings = self.forms[f].check( info )
            if formWarnings:
                warnings.append( warning.context( "Form %s:" % f, formWarnings ) )

        return warnings
    
    def fix( self, **kw ): 
        '''assemble some useful info, to be handed down to children's check() functions'''
        info = {}
        warnings = []
        if self.rowData == None:    info['rowENDF_MFMT'] = None
        else:                       info['rowENDF_MFMT'] = self.rowData['ENDF_MFMT']
        if self.columnData == None: info['columnENDF_MFMT'] = None
        else:                       info['columnENDF_MFMT'] = self.columnData['ENDF_MFMT']
        info.update( kw )
        for f in self.forms: warnings += self.forms[f].fix( **info )
        return warnings

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<section label="%s" id="%s" nativeData="%s"' % (self.label, self.id, self.nativeData)]
        if self.columnData is not None: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData','columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append( getattr(self, dataPointer).toXML(indent+'  ') )
        for key in self.forms.keys():
            xmlString += self.forms[key].toXMLList(flags, verbosityIndent, indent+'  ')
        xmlString[-1] += '</section>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        rowData = link.parseXMLNode( element[0], linkData )
        linkData['unresolvedLinks'].append( rowData )
        columnData = None
        if element[1].tag=="columnData":
            columnData = link.parseXMLNode( element[1], linkData )
            linkData['unresolvedLinks'].append( columnData )
        section_ = section( element.get('label'), element.get('id'), rowData, columnData, element.get('nativeData') )
        start = 2 if (columnData is not None) else 1
        for form in element[start:]:
            formClass = {
                    covarianceFormToken: covarianceMatrix,
                    mixedFormToken: mixedForm,
                    summedFormToken: summedCovariance,
                    energyIntervalFormToken: energyIntervalForm,
                    legendreOrderCovarianceFormToken: LegendreOrderCovarianceForm,
                    }.get( form.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.addForm( formClass.parseXMLNode( form, linkData ) )
        return section_
            



class covarianceMatrix( ancestry ):
    """
    Base class for all covariances. covarianceMatrix contains axes, a list of energy boundaries,
    and the matrix data. Some details :

        * Matrix data is stored in an :py:class:`xData.matrix class`. May be diagonal, symmetric, sparse, etc
        * Symmetric matrices only require one set of energy bounds, but asymmetric
          matrices require bounds for both axes.
          
    """

    moniker = covarianceFormToken

    def __init__(self, index=None, type=absoluteToken, axes=None, matrix=None, energyBounds=None,
            ENDFconversionFlag=None):
        ancestry.__init__( self, 'covarianceMatrix', None, attribute = 'index' )
        self.index = index  #: an int, required inside a 'mixed' section
        self.type = type #: 'relative' or 'absolute'
        self.axes = axes or [] #: list of covarianceAxis instances.  There should be 3: the rows, the columns and the value of the covariance itself
        self.matrix = matrix or [] #: a :py:class:`fudge.core.math.matrix` instance containing the actual covariance matrix
        self.energyBounds = energyBounds #: duh, the energy bounds
        self.ENDFconversionFlag = ENDFconversionFlag #: yes, this is a crutch to help when converting back to ENDF

    def convertAxesToUnits( self, units ):
        '''
        Converts all the axes' units.
        The parameter ``units`` should be a list of units with the same length as self.axes
        '''
        if not type( units ) in [ list, tuple ]: raise TypeError()
        if len( units ) != len( self.axes ): raise ValueError()
        for i,a in enumerate( self.axes ): a.convertToUnit( units[i] )

    def toCovarianceMatrix( self ): 
        import copy
        return copy.copy( self )
        
    def toCorrelationMatrix( self ):
        '''
        Returns the correlation matrix generated from self's covariance matrix.  This is
        essentially a copy of self, but renormalized by the uncertainty:
        
            correlation[i,j] = covariance[i,j]/sqrt(covariance[i,i])/sqrt(covariance[j,j])
        
        We reuse the covariance matrix class so that we can do plotting, etc.  If you have
        a correlation matrix, you can safely recover it provided you have the uncertainty 
        vector.
        
        Currently only works for a square covariance matrix and not a off-diagonal part of 
        another covariance.
        '''
        # Check if is full, square covariance matrix
        if not self.matrix.form in (gndMatrix.symmetricFormToken, 
            gndMatrix.diagonalFormToken,gndMatrix.sparse_symmetricFormToken): 
                raise TypeError( "Can only extract correlation matrices from symmetric covariance matrices" )
                
        # Copy to result matrix
        import copy
        correlation = copy.copy( self )
        
        # Rescale result matrix
        import numpy
        theCorrelationMatrix = numpy.matrix( self.matrix.data )
        theUncertainty = numpy.diag( theCorrelationMatrix )
        theUncertainty[ theUncertainty < 0.0 ] = 0.0
        theUncertainty = numpy.sqrt( theUncertainty )
        for i in range( theCorrelationMatrix.shape[0] ):
            for j in range( theCorrelationMatrix.shape[1] ):
                theCorrelationMatrix[i,j] /= ( theUncertainty[i] * theUncertainty[j] )

        # Return the result
        correlation.axes[2].unit = ''
        correlation.matrix.data = theCorrelationMatrix.tolist()
        return correlation
        
    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param XYs rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param XYs colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note:    If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        import copy
        result = copy.copy( self )
        if self.type==absoluteToken: return result

        # Make sure we have usable data to rescale with
        if not isinstance( rowData, XYs.XYs ): raise TypeError( 'rowData must be of type XYs, found '+str(type(rowData)) )
        gRowData = rowData.group( self.axes[0].data, self.axes[0].unit )
        if isinstance( colData, XYs.XYs ): gColData = colData.group( self.axes[1].data, self.axes[1].unit )
        else: gColData = gRowData
        
        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.data):
            newRow = []
            for j, cell in enumerate(row):
                newRow.append( gRowData[i][1]*gColData[j][1]*cell )
            newData.append( newRow )
        result.matrix.data = newData
        result.type=absoluteToken
        return result
        
    def toRelative( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a absolute covariance) using XYs rowData and colData
        to convert self into a relative covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        import copy
        result = copy.copy( self )
        if self.type==relativeToken: return result

        # Make sure we have usable data to rescale with
        if not isinstance( rowData, XYs.XYs ): raise TypeError( 'rowData must be of type XYs, found '+str(type(rowData)) )
        gRowData = rowData.group( self.axes[0].data, self.axes[0].unit )
        if not self.axes[1].mirrorOtherAxis:
            if not isinstance( colData, XYs.XYs ): raise TypeError( 'colData must be of type XYs, found '+str(type(colData)) )
            gColData = colData.group( self.axes[1].data, self.axes[1].unit )
        else: gColData = gRowData

        # Rescale!
        newData = []
        for i,row in enumerate(self.matrix.data):
            newRow = []
            for j, cell in enumerate(row):
                newRow.append( cell/gRowData[i][1]/gColData[j][1] )
            newData.append( newRow )
        result.matrix.data = newData
        result.type=relativeToken
        return result
        
    def check( self, info ): 
        """Check if uncertainty in the bounds passed into the checker.  
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        from fudge.gnd import warning
        warnings = []

        if self.matrix.form in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken) and info['checkUncLimits']: 
            import numpy
            A = numpy.matrix( self.matrix.data )
            relative = self.type == 'relative'
            if relative:
                varMin = info['minRelUnc']*info['minRelUnc']
                varMax = info['maxRelUnc']*info['maxRelUnc']
            for i in range( A.shape[0] ):
                if not relative:
                    if info['theData'] != None:
                        uncMin = info['minRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[0].data[i]+self.axes[0].data[i+1]) )
                        uncMax = info['maxRelUnc'] * info['theData'].getValue(
                                0.5*(self.axes[0].data[i]+self.axes[0].data[i+1]) )
                        varMin = uncMin*uncMin
                        varMax = uncMax*uncMax
                    else:
                        #warnings.append( "WARNING: can't check absolute uncertainties without data to compare to!\n" )
                        break
                if varMin <= A[i,i] and varMax >= A[i,i]: pass # unc is where is should be
                elif varMin >= A[i,i]:                         # uncertainty too small
                    warnings.append( warning.varianceTooSmall( i, A[i,i], self ) )
                else:                                          # uncertainty too big
                    warnings.append( warning.varianceTooLarge( i, A[i,i], self ) )
        return warnings + self.matrix.check( info )
    
    def fix( self, **kw ): 
        """Fix uncertainty using the bounds passed into the fixer.
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        warnings = []
        if self.matrix.form in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken) and kw['fixUncLimits']: 
            import numpy
            A = numpy.matrix( self.matrix.data )
            relative = self.type == 'relative'
            for i in range( A.shape[0] ):
                eMin = self.axes[0].data[i]
                eMax = self.axes[0].data[i+1]
                if relative:
                    uncMin = kw['minRelUnc']
                    uncMax = kw['maxRelUnc']
                else:
                    uncMin = kw['theData'].getValue( 0.5*(eMax+eMin) )
                    uncMax = kw['theData'].getValue( 0.5*(eMax+eMin) )
                uncMin2 = uncMin*uncMin
                uncMax2 = uncMax*uncMax
                #eThresh = threshold.getValueAs( component.axes[0].units )
                if uncMin2 <= A[i,i] and uncMax2 >= A[i,i]: pass   # unc is where is should be
                elif uncMin2 >= A[i,i]: 
                    if i+1 < A.shape[0] and A[i+1,i+1] >= uncMin2: A[i,i] = A[i+1,i+1]
                    else: A[i,i] = uncMin2
                # else:                                            # above threshold and uncertainty out of bounds
                    # continue #skip this fix for now
                    # if uncMin2 > A[i,i]: # uncertainty too small
                        # print '    WARNING: bin', i, 'uncertainty is too small!!!', '(', uncMin2, '>', A[i,i], ')'
                        # if A[i, i] == 0.0:
                            # for j in range( i+1, A.shape[0] ): 
                                # A[i,j] = 0.0
                                # A[j,i] = 0.0
                            # A[i,i] == uncMin2
                        # else:
                            # for j in range( i+1, A.shape[0] ): 
                                # A[i,j] *= uncMin / math.sqrt( A[i,i] )
                                # A[j,i] = A[i,j]
                            # A[i,i] == uncMin2
                    # else:                # uncertainty too big
                        # print '    WARNING: bin', i, 'uncertainty is too big!!!', '(', A[i,i], '>', uncMax2, ')'
                        # for j in range( i+1, A.shape[0] ): 
                            # A[i,j] *= uncMax / math.sqrt( A[i,i] )
                            # A[j,i] = A[i,j]
                        # A[i,i] = uncMax2
            self.data = B.tolist()
        return warnings + self.matrix.fix( **kw )
        
    def group( self, groupBoundaries = ( None, None ), groupUnit = ( None, None ) ):
        '''
        Group the matrix in self
        
        :param original: the original covarianceMatrix instance we intend to regroup
        :param groupBoundaries: a 2 element list containing the group boundaries for the rows 
                                and columns (respectively) of the covariance to be regrouped
        :param groupUnit: a 2 element list containing the units in which group boundaries are 
                          specified for the rows and columns (respectively) of the covariance 
                          to be regrouped
            
        .. note::  We still need to do flux weighting
            
            
        .. rubric:: Regrouping Theory
        
        Given a function :math:`f(E)`, we write the grouped data using fudge's ``flat`` interpolation 
        scheme.  We note that we could write this scheme as an expansion over basis functions:
        
        .. math::    
            f(E) = \sum_{i=0}^{N+1} w_i(E) * f_i
        
        where the weight functions :math:`w_i(E)` are
        
        .. math::
            w_i(E) = 1  \;\\text{for}\; E_i <= E <= E_{i+1}; \;\; 0 \;\\textrm{otherwise}
            
        These weights are an orthogonal (bot not orthonormal) basis, with 
        
        .. math::
            (E_{i+1}-E_i) \delta_{ij} = \int dE w_i(E) * w_j(E)
        
        So, to transform from basis :math:`w_i(E)` to :math:`v_i(E)` (which has group boundaries 
        :math:`[ E'_0, ... ]`), do: 
        
        .. math::
            f'_j = \sum_i m_{ji} f_i
            
        where :math:`f'` is the regrouped function coefficients and :math:`m_{ji}` is the matrix
        
        .. math::
            m_{ij} = (E'_{i+1}-E'_i)^{-1} \int dE v_i(E) w_j(E) 

            
        .. rubric:: Applying regrouping theory to covariance matrices   
        
        When we are given a covariance matrix :math:`c_{ij}` in ENDF, it is meant to be interpreted
        as a grouped covariance in both the direction of the matrix rows and the matrix 
        columns.  Therefore, we must regroup in both the row direction and the column 
        direction.  The ENDF format gives both the group boundaries for the rows and columns.
        In other words, ENDF gives us the following rule for evaluating the continuous row-
        column covariance:
        
        .. math::
            c( E1, E2 ) = \sum_{ij} w_i(E1) w_j(E2) c_{ij}
            
        Computing :math:`m_{ij}` as before, 
            
        .. math::
            cc_{ij} = \sum_{i',j'} m_{ii'} c_{i'j'} m_{j'j}
            
        It is straightforward to generalize to the case where the row and column bases are 
        different.
        
        In the routine below, we abuse :py:class:`fudge.core.math.xData.XYs` to specify the functions 
        :math:`w_i(E)` and use the :py:func:`XYs.groupOneFunction()` method to perform the integrals to get
        the regrouping matrix.  We do this separately for the rows and the columns.
        The matrix multiplication that converts a covariance from one pair of bases (group 
        structures) to another is accomplished using numpy.

        
        .. rubric:: An explanation of fudge's 'flat' interpolation
        
        Suppose we have a function :math:`f(E)` specified using fudge's `'flat'` interpolation.  
        Then we have :math:`N` entries :math:`[f_0, f_1, ..., f_{N-1}]` and a set of group 
        boundaries :math:`[E_0, E_1, ..., E_N]` and the following rule for interpolation:
        
            * Below :math:`E_0`, :math:`f(E)` evaluates to :math:`0.0`
            * From :math:`E_0 \\rightarrow E_1`, :math:`f(E)` evaluates to :math:`f_0`
            * From :math:`E_1 \\rightarrow E_2`, :math:`f(E)` evaluates to :math:`f_1`
            * ...
            * From :math:`E_{i} \\rightarrow E_{i+1}`, :math:`f(E)` evaluates to :math:`f_i`
            * ...
            * From :math:`E_{N-1} \\rightarrow E_N`, :math:`f(E)` evaluates to :math:`f_{N-1}`
            * Above :math:`E_N`, :math:`f(E)` evaluates to :math:`0.0`
        '''
        import numpy, copy

        # determine where to get the settings for the potentially mirrored second axis
        if self.axes[1].mirrorOtherAxis:    axis1index = 0
        else:                               axis1index = 1
        
        # setup the old axes in a form we can (ab)use in the XYs class
        axes0_ = axes.defaultAxes( 
            labelsUnits={ 
                0:( self.axes[0].label, self.axes[0].unit ),
                1:( 'dummy', '' )},
            dependentInterpolation='flat' )        
        axes1_ = axes.defaultAxes( 
            labelsUnits={ 
                0:( self.axes[axis1index].label, self.axes[axis1index].unit ),
                1:( 'dummy', '' )},
            dependentInterpolation='flat' )        
        
        # define basis functions for the rows and columns
        basis0 = XYs.XYs( axes0_, [ ( x, 0.0 ) for x in self.axes[0].data ], 0.0001 )
        basis1 = XYs.XYs( axes1_, [ ( x, 0.0 ) for x in self.axes[axis1index].data ], 0.0001 )
        basis0 = basis0.convertAxisToUnit( 0, groupUnit[0] )
        basis1 = basis1.convertAxisToUnit( 0, groupUnit[1] )
    
        # build the regrouping matrices for the two bases
        w0 = []
        for i in range( self.matrix.nrows ):
            basis0[i] = ( basis0[i][0], 1.0 )
            w0.append( basis0.groupOneFunction( groupBoundaries[0], norm = 'dx' ) )
            basis0[i] = ( basis0[i][0], 0.0 )
        w0 = numpy.mat( w0 )
        w1 = []
        for j in range( self.matrix.ncols ):
            basis1[j] = ( basis1[j][0], 1.0 )
            w1.append( basis1.groupOneFunction( groupBoundaries[1], norm = 'dx' ) )
            basis1[j] = ( basis1[j][0], 0.0 )
        w1 = numpy.mat( w1 )
                
        # set up the regrouped covariance matrix
        grouped = copy.copy( self )
        grouped.axes[0].data = groupBoundaries[0]
        grouped.axes[1].data = groupBoundaries[1]
        grouped.axes[0].unit = groupUnit[0]
        grouped.axes[1].unit = groupUnit[1]
        odata = numpy.mat( self.matrix.data )
        gdata = w0.T * odata * w1
        grouped.matrix.data = gdata.tolist()
        grouped.matrix.nrows = len(grouped.matrix.data)
        grouped.matrix.ncols = len(grouped.matrix.data[0]) 
        return grouped

    def removeExtraZeros(self):
        import numpy
        matrix = numpy.array( self.matrix.data )
        rowStart, colStart = 0,0
        rowEnd, colEnd = matrix.shape
        while numpy.all(matrix[rowStart,:]==0):
            rowStart += 1
        while numpy.all(matrix[:,colStart]==0):
            colStart += 1
        while numpy.all(matrix[rowEnd-1,:]==0):
            rowEnd -= 1
        while numpy.all(matrix[:,colEnd-1]==0):
            colEnd -= 1

        matrix = matrix[rowStart:rowEnd, colStart:colEnd]
        if len(self.axes)==1:
            assert (rowStart,rowEnd) == (colStart,colEnd)
            self.axes[0].data = self.axes[0].data[rowStart:rowEnd]
        else:
            self.axes[0].data = self.axes[0].data[rowStart:rowEnd]
            self.axes[1].data = self.axes[1].data[colStart:colEnd]
        self.matrix.data = matrix.tolist()
        self.matrix.nrows, self.matrix.ncols = matrix.shape

    def getUncertaintyVector( self, theData=None, relative=True ):
        """ 
        Get an XYs object containing uncertainty for this matrix.
        Convert relative/absolute if requested (if so, must also pass central values as theData)

        Examples:

            - if the covariance matrix is relative and we want relative uncertainty vector, just do:
            
                >>> matrix.getUncertaintyVector()
                
            - if we want the absolute matrix instead:
            
                >>> matrix.getUncertaintyVector( theData=<XYs instance>, relative=False ) 
                
        """

        if self.matrix.form not in (gndMatrix.symmetricFormToken, gndMatrix.diagonalFormToken,
                gndMatrix.sparse_symmetricFormToken):
            raise ValueError("getUncertaintyVector only applies to symmetric matrices!")
        import numpy
        energies = self.axes[0].data
        diag = list( numpy.sqrt( numpy.diagonal( numpy.array( self.matrix.data ) ) ) )
        diag.append( diag[-1] )                             # point corresponding to final energy bin
        yunit = self.axes[-1].unit
        if yunit != '': # get square root of the unit
            yunit = PhysicalQuantityWithUncertainty(1,yunit).sqrt().getUnitSymbol()
        axes_ = axes.defaultAxes( labelsUnits={0:('energy_in',self.axes[0].unit),1:('uncertainty',yunit)},
                dependentInterpolation='flat' )
        uncert = XYs.XYs( axes_, zip(energies,diag), 0.0001 )    # what should accuracy be set to?
        uncert = uncert.changeInterpolation('linear','linear',None,1e-8,1e-8)

        # do we need to convert absolute->relative or vice versa?
        if (relative and self.type==absoluteToken) or (not relative and self.type==relativeToken):
            if theData==None:
                raise ValueError("Need central values to convert relative<->absolute uncertainties!")
            try:
                theData = theData.toPointwiseLinear(1e-8,1e-8)
                uncert, theData = uncert.mutualify(1e-8, 1e-8, False, theData, 1e-8, 1e-8, False)
                if relative: #convert absolute to relative
                    uncert /= theData
                else: #relative to absolute
                    uncert *= theData
            except Exception as err:
                print len( uncert ), uncert.copyDataToXYs()[0], uncert.copyDataToXYs()[-1]
                print len( theData ), theData.copyDataToXYs()[0], theData.copyDataToXYs()[-1]
                raise Exception( err.message )
        return uncert
        
    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import QuadMesh
        
        x = self.axes[0].data
        if self.axes[1].mirrorOtherAxis:    y = self.axes[0].data 
        else:                               y = self.axes[1].data 

        X, Y = np.meshgrid( x, y )        
        XY = np.hstack((X.ravel()[:,np.newaxis], Y.ravel()[:,np.newaxis]))
        Z = (np.array(self.matrix.data)).ravel()
        
        ax = plt.subplot(1,1,1)
        if title == None: title = str( self.toXLink() )
        plt.suptitle(title)

        qc = QuadMesh(
            meshWidth=len(x)-1, 
            meshHeight=len(y)-1, 
            coordinates=XY, 
            #showedges=False, 
            antialiased=True, 
            shading='flat',
            transOffset=ax.transData)
            
        qc.set_array(Z) 
        ax.add_collection(qc,autolim=True)

        if xlim == None:    ax.set_xlim( x[0], x[-1] )
        else:               ax.set_xlim( xlim[0], xlim[1] )
        if ylim == None:    ax.set_ylim( y[0], y[-1] )
        else:               ax.set_ylim( ylim[0], ylim[1] )
        if xlog: ax.set_xscale( 'log' )
        if ylog: ax.set_yscale( 'log' )

        xlabel = self.axes[0].label + ' (' + self.axes[0].unit +')'
        if self.axes[1].mirrorOtherAxis:    ylabel = xlabel 
        else:                               ylabel = self.axes[1].label + ' (' + self.axes[1].unit +')'

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        cbar = plt.colorbar(qc)
        if scalelabel != None: cbar.set_label(scalelabel)
        else: cbar.set_label(str(self.type)+' covariance ('+str(self.axes[2].unit)+')')
        plt.show()

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        indent2 = indent+'  '; indent3 = indent2+'  '
        xmlString = [indent+'<%s' % self.moniker]
        if self.index!=None: xmlString[0] += ' index="%s"' % self.index
        xmlString[0] += ' type="%s"' % self.type
        if self.energyBounds: xmlString[0] += (
                ' lowerBound="%s" upperBound="%s"' %
                self.energyBounds )
        if self.ENDFconversionFlag: xmlString[0] += (
                ' ENDFconversionFlag="%s"' % self.ENDFconversionFlag )
        xmlString[0] += '>'
        xmlString.append( '%s<axes>' % indent2 )
        for axis in self.axes: xmlString.append( axis.toXML(indent=indent3) )
        xmlString[-1] += '</axes>'
        xmlString += self.matrix.toXMLList(flags,verbosityIndent,indent2)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode(xmlElement, linkData={}):
        """Read one <covarianceMatrix> element from xml into python class """
        ax, mat = xmlElement[0], xmlElement[1]
        axes_ = covarianceAxis.parseXMLNode( ax )
        matrix_ = gndMatrix.parseXMLNode( mat )
        return covarianceMatrix( type=xmlElement.get('type'), axes=axes_, matrix=matrix_,
                ENDFconversionFlag=xmlElement.get("ENDFconversionFlag") )
    
    def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
        endf = []
        rowdat, coldat = targetInfo['dataPointer']
        MF,MT1 = map(int, rowdat['ENDF_MFMT'].split(','))
        if not inCovarianceGroup:
            # print header for this subsection (contains one NL sub-subsection)
            MAT1 = targetInfo['MAT1']
            XMF1,XLFS1,NC,NI = 0,0,0,1
            if coldat:
                MF1, MT1 = map(int, coldat['ENDF_MFMT'].split(','))
            if MF in (31,33):
                endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        # header for matrix:
        rows,cols = self.matrix.nrows, self.matrix.ncols
        if self.matrix.form==gndMatrix.diagonalFormToken:
            LS = 0; LB = 1; NP = len(self.axes[0]); NT = 2*NP
            if self.type=='absolute': LB = 0
            if self.ENDFconversionFlag:
                LB = int( self.ENDFconversionFlag.split('=')[1] )
            matrixData = []
            for p in range(NP-1):
                matrixData.extend( [self.axes[0][p], self.matrix[p][p] ] )
            matrixData.extend( [self.axes[0][NP-1], 0] )
        elif self.matrix.form==gndMatrix.symmetricFormToken:
            LS = 1; LB = 5; NT = (rows+1) + rows*(rows+1)/2; NP = rows+1
            matrixData = self.axes[0].data + [
                    i for j in range(rows) for i in self.matrix.data[j][j:]]
        elif self.axes[1].mirrorOtherAxis:
            LS = 0; LB = 5; NT = (rows+1) + rows*cols; NP = rows+1
            matrixData = self.axes[0].data + [
                    i for sublist in self.matrix.data for i in sublist]
        else:
            LS = 0; LB = 6; NT = (rows+1) + (cols+1) + rows*cols; NP = rows+1
            matrixData = self.axes[0].data + self.axes[1].data + [
                    i for sublist in self.matrix.data for i in sublist]
        if MF==35:  # header for fission spectra is different:
            E1,E2 = [a.getValueAs('eV') for a in self.energyBounds]
            if LS: LB = 7
            else:
                raise Exception ("Unknown spectrum (MF35) covariance format")
            endf.append( endfFormats.endfHeadLine( E1,E2,LS,LB,NT,NP ) )
        else:
            endf.append( endfFormats.endfHeadLine( 0,0,LS,LB,NT,NP ) )
        endf += endfFormats.endfDataList( matrixData )
        return endf

class covarianceAxis:
    """
    Store energy group boundaries (similar to xData.axes.axis).
    Covariances in ENDF are often symmetric, NxN matrices. The
    axis then must have N+1 points defining upper and lower bounds for each
    bin in the matrix. 
    """
    def __init__(self, index=0, label=None, unit=None, interpolation=None,
            data=None, mirrorOtherAxis=False):
        self.index = index #: an int
        self.label = label #: a str to label the axis
        self.unit = unit #: any acceptable unit for a PhysicalQuantityWithUncertainty
        self.data = data or [] #: a list of floats in unit "unit" denoting the nodes (in ENDF, they are all group boundaries) for interpolating the rows and or columns of a covariance matrix
        self.interpolation = interpolation #: a string describing how the the covariance is meant to be interpreted.  Currently all ENDF covariance is 'linear,flat'
        self.mirrorOtherAxis = mirrorOtherAxis #: if True, then self is a copy of the accompanying row/column axis

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]
        
    def convertToUnit( self, unit ):
        '''
        Converts the units of a covarianceAxis, but only if the covarianceAxis isn't a mirror of another one
        '''
        if self.mirrorOtherAxis: return
        self.data = [ PhysicalQuantityWithUncertainty( x, self.unit ).inUnitsOf( unit ).getValue() for x in self.data ]
        self.unit = unit

    def toXML(self, indent = ''):
        from pqu.physicalQuantityWithUncertainty import toShortestString
        interpolationFlag, axisFlag, lengthFlag = '','',''
        if self.interpolation: interpolationFlag=' interpolation="%s"' % self.interpolation
        if self.mirrorOtherAxis: axisFlag=' mirror_row_energy_bounds="true"'
        if len(self.data)!=0: lengthFlag = ' length="%i"' % len(self.data)
        xmlString = indent+'<axis index="%i" label="%s" unit="%s"%s%s%s>' % (
                self.index, self.label, self.unit, interpolationFlag, axisFlag, lengthFlag)
        # self.data is usually list of floats (energy boundaries). Could also be list of input parameters
        if len(self.data)!=0:
            xmlString += ''.join([' %s' % toShortestString(val) for val in self.data])
            xmlString += '</axis>'
        else: xmlString = xmlString[:-1] + '/>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        axes_ = []
        for axis in element:
            data = None
            if axis.text is not None: data = map(float, axis.text.split())
            mirrorOtherAxis = False
            if axis.get('mirror_row_energy_bounds') in ('true','True'): mirrorOtherAxis=True
            interp = axis.get('interpolation')
            if interp is not None: interp = axes.interpolationXY( *interp.split(',') )
            axes_.append( covarianceAxis( index=int(axis.get('index')), label=axis.get('label'),
                unit=axis.get('unit'), data=data, interpolation=interp, mirrorOtherAxis=mirrorOtherAxis ) )
        return axes_

class mixedForm( ancestry ):
    """
    Covariance for a single quantity, stored as several separate matrices that must be summed together.
    In general, the energy bounds for these matrices can overlap (unlike piecewise cross section data). 
    """

    moniker = mixedFormToken

    def __init__(self, components=None):
        ancestry.__init__( self, 'mixedForm', None )
        self.components = components or [] #: a Python list containing instances of ``mixedForm``, ``summedCovariance``, and ``covarianceMatrix``

    def __getitem__(self, idx):
        return self.components[idx]

    def __len__(self):
        return len(self.components)

    def check( self, info ): 
        from fudge.gnd import warning
        warnings = []

        for component in self.components:
            componentWarnings = component.check( info )
            if componentWarnings:
                warnings.append( warning.context('Component %i' % component.index, componentWarnings) )
        
        #if warnings:
        #    warnings = [warning.context('Section "%s": %s' % (self.id, form), warnings)]
        return warnings
        
    def fix( self, **kw ): 
        warnings = []
        for comp in self.components: warnings += comp.fix( **kw )
        return warnings

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def addComponent(self, covariance):
        ''':param covariance: an instance of covariance (or inherited class)'''
        covariance.setParent(self)
        self.components.append(covariance)

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combines all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs instance
        """

        xys = [ sec.getUncertaintyVector( theData=theData, relative=relative ) for sec in self.components ]
        uncert = xys[0]
        for xy in xys[1:]:
            uncert, xy = uncert.mutualify( 1e-8, 0, 0, xy, 1e-8, 0, 0 )
            uncert += xy
        return uncert
        
    def toCovarianceMatrix( self ): 
        if len( self.components ) == 1: return self.components[0].toCovarianceMatrix()
        from fudge.core.utilities.brb import uniquify
        import numpy, copy
        
        # set up common data using first component
        firstCovMtx = self.components[0].toCovarianceMatrix()
        commonRowAxis = copy.copy( firstCovMtx.axes[0] )
        if firstCovMtx.axes[1].mirrorOtherAxis: commonColAxis = copy.copy( firstCovMtx.axes[0] )
        else:                                   commonColAxis = copy.copy( firstCovMtx.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.axes[2] )
        commonType = firstCovMtx.type
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
            if cc.type != commonType: raise ValueError( "Incompatable types in "+str(self.__class__)+": "+str(commonType)+' vs. '+str(cc.type) )
            cc.convertAxesToUnits( ( commonRowAxis.unit, commonColAxis.unit, commonMatrixAxis.unit ) )
            commonRowAxis.data = commonRowAxis.data + cc.axes[0].data
            if cc.axes[1].mirrorOtherAxis: commonColAxis.data = commonColAxis.data + cc.axes[0].data
            else:                          commonColAxis.data = commonColAxis.data + cc.axes[1].data
        commonRowAxis.data.sort()
        commonColAxis.data.sort()
        commonRowAxis.data = uniquify( commonRowAxis.data )
        commonColAxis.data = uniquify( commonColAxis.data )
        
        # now sum up the components
        commonMatrix = numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for c in self.components[1:]:
            cc = c.toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        return covarianceMatrix( type=commonType, axes=[ commonRowAxis, commonColAxis, commonMatrixAxis ], matrix=gndMatrix.matrix( commonMatrix.tolist(), form=gndMatrix.symmetricFormToken ) )

    def toAbsolute( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a relative covariance) using XYs rowData and colData
        to convert self into an absolute covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    absolute because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to absoluteToken
        '''
        import copy
        result = copy.copy( self )
        result.components = []
        for c in self.components:
            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toAbsolute() call with the correct row or column data
                result.components.append( c.toCovarianceMatrix().toAbsolute( rowData, colData ) )
            else: result.components.append( c.toAbsolute( rowData, colData ) )
        return result
        
    def toRelative( self, rowData=None, colData=None ): 
        '''
        Rescales self (if it is a absolute covariance) using XYs rowData and colData
        to convert self into a relative covariance matrix.
        
        :param rowData: an XYs instance containing data to rescale covariance in the "row direction"
        :param colData: an XYs instance containing data to rescale covariance in the "col direction"
            
        .. note::   If the column axis is set to 'mirrorOtherAxis', only rowData is needed.  
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already 
                    relative because this will throw an error.
            
        :returns: a copy of self, but rescaled and with the type set to relativeToken
        '''
        import copy
        result = copy.copy( self )
        result.components = []
        for c in self.components:
            if isinstance( c, summedCovariance ): 
                # If the covariance is summed, a call to toCovarianceMatrix() should add up
                # the pointed-to covariances (if they are the same type (relative vs. absolute)),
                # allowing us to do a toRelative() call with the correct row or column data
                result.components.append( c.toCovarianceMatrix().toRelative( rowData, colData ) )
            else: result.components.append( c.toRelative( rowData, colData ) )
        return result

    def getSingleMatrix(self):
        """Combine all matrices in a mixed section into a single matrix.
        Beware: resulting matrix may be huge! """
        from fudge.core.math.linalg import expand_grid
        if not all( [isinstance(a, covarianceMatrix) for a in self.components] ):
            raise Exception( "I don't know how to sum these matrices!" )
        supergrid = set()
        for component in self:
            for axis in component.axes:
                supergrid.update( axis.data )
        supergrid = sorted( supergrid )
        matrix = numpy.zeros( (len(supergrid)-1,len(supergrid)-1) )
        for component in self:
            erows, ecols = [a.data for a in component.axes[:2]]
            matrix += expand_grid(supergrid, component.matrix.data, erows, ecols)
        return supergrid, matrix

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<%s>' % self.moniker]
        for covariance in self.components:
            xmlString += covariance.toXMLList(flags,verbosityIndent,indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, linkData={} ):
        mixed_ = cls()
        for child in element:
            formClass = {
                    covarianceFormToken: covarianceMatrix,
                    summedFormToken: summedCovariance,
                    }.get( child.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % child.tag)
            mixed_.addComponent( formClass.parseXMLNode( child, linkData ) )
        for i,form in enumerate( mixed_.components ): form.index = i
        return mixed_
    
    def toENDF6(self, flags, targetInfo):
        endf = []
        XMF1,XLFS1,MAT1 = 0,0,0
        NI = len([cov for cov in self.components if hasattr(cov,'matrix')])
        NC = len(self.components) - NI
        rowdat, coldat = targetInfo['dataPointer']
        MF, MT1 = map(int, rowdat['ENDF_MFMT'].split(','))
        if MF in (31,33):
            endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        for cov in self.components:
            endf += cov.toENDF6(flags, targetInfo, inCovarianceGroup=True)
        return endf

class energyIntervalForm(mixedForm):
    """ 
    For distributions, the covariances may depend on incident energy. Similar to mixed,
    but each component has an associated energy 
    """

    moniker = energyIntervalFormToken

    def __init__(self, components=None):
        self.components = components or [] #: a list of components that are instances of ``mixedForm`` or ``covarianceMatrix``

    @classmethod
    def parseXMLNode( cls, element, linkData={} ):
        form = super(energyIntervalForm, cls).parseXMLNode( element, linkData )
        # add energy bounds to each component:
        for i in range(len(element)):
            component = element[i]
            lower = PhysicalQuantityWithUncertainty( component.get("lowerBound") )
            upper = PhysicalQuantityWithUncertainty( component.get("upperBound") )
            form.components[i].energyBounds = (lower,upper)
        return form

class LegendreOrderCovarianceForm( ancestry ):
    """ 
    Stores covariance between energy-dependent Legendre coefficients for a reaction.
    This class contains one or more LegendreLValue sections, each section containing the matrix
    between a pair of L-values 
    """
    
    moniker = legendreOrderCovarianceFormToken

    def __init__(self, lvalues=None):
        ancestry.__init__( self, legendreOrderCovarianceFormToken, None )
        self.lvalues = lvalues or [] #: the l values of course

    def __getitem__(self, index):   return self.lvalues[index]

    def addLegendreOrder( self, LValue ):
        LValue.setParent( self )
        self.lvalues.append( LValue )

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<%s>' % self.moniker]
        for lvalue in self.lvalues:
            xmlString += lvalue.toXMLList(flags,verbosityIndent,indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def toENDF6(self, flags, targetInfo):
        endf = []
        for lVal in self.lvalues:
            LCT = {'frameOfMF4':0, 'lab':1, 'centerOfMass':2}[ lVal.frame ]
            form = lVal.forms[ lVal.nativeData ]
            NI = 1
            if form.moniker == mixedFormToken: NI = len(form)
            endf.append( endfFormats.endfHeadLine(0,0,lVal.L1, lVal.L2, LCT, NI) )
            endf += form.toENDF6(flags, targetInfo)
        return endf

    @classmethod
    def parseXMLNode( cls, element, linkData={} ):
        form = LegendreOrderCovarianceForm()
        # add L's to each component:
        for lValue in element:
            form.addLegendreOrder( LegendreLValue.parseXMLNode( lValue, linkData ) )
        return form


class LegendreLValue( ancestry ):
    """ 
    Represents one subsection of the Legendre coefficient covariance matrix:
    covariance between coefficients for two Legendre orders at various energies 
    """

    moniker = legendreLValueFormToken

    def __init__(self, L1, L2, frame, nativeData=None):
        ancestry.__init__( self, legendreLValueFormToken, None )
        self.L1 = L1 #:
        self.L2 = L2 #:
        self.frame = frame #:
        self.nativeData = nativeData #:
        self.forms = {} #:

    def addForm( self, form ):
        form.setParent( self )
        self.forms[ form.moniker ] = form

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        xmlString = [indent+'<%s L1="%i" L2="%i" frame="%s" nativeData="%s">' % (self.moniker,
            self.L1, self.L2, self.frame, self.nativeData)]
        for form in self.forms.values():
            xmlString += form.toXMLList(flags,verbosityIndent,indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, linkData={} ):
        component = LegendreLValue( int(element.get("L1")), int(element.get("L2")),
                element.get("frame"), element.get("nativeData") )
        for form in element:
            formCls = {
                    covarianceFormToken: covarianceMatrix,
                    mixedFormToken: mixedForm
                    }.get( form.tag )
            if formCls is None:
                raise BadCovariance
            component.addForm( formCls.parseXMLNode( form, linkData ) )
        return component

class summedCovariance( ancestry ):
    """ 
    Covariance matrix stored as sum/difference of other matrices. 
    """

    moniker = summedFormToken

    def __init__(self, index=None, lowerBound=None, upperBound=None, 
            coefficients=None, pointerList=None):
        ancestry.__init__( self, summedFormToken, None, attribute = 'index' )
        self.index = index #: an int, only needed within 'mixed' section
        self.lowerBound = lowerBound #: lower bound of row/column direction.  (type: PhysicalQuantityWithUncertainty)
        self.upperBound = upperBound #: upper bound of row/column direction.  (type: PhysicalQuantityWithUncertainty)
        self.pointerList = pointerList or [] #a list of pointers (type: fudge.gnd.link).  Each pointer has a dictionary with entry "coefficient" to use to weight the covariance referred to

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def plot( self, title = None, scalelabel = None, xlim=None, ylim=None, xlog=False, ylog=False ):
        self.toCovarianceMatrix().plot( title=title, scalelabel=scalelabel, xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog )

    def toCovarianceMatrix( self ): 
        if len( self.pointerList ) == 1: return self.pointerList[0].link.getNativeData().toCovarianceMatrix()
        from fudge.core.utilities.brb import uniquify
        import numpy, copy
        
        # set up common data using first element in pointerList
        firstCovMtx = self.pointerList[0].link.getNativeData().toCovarianceMatrix()
        commonRowAxis = copy.copy( firstCovMtx.axes[0] )
        if firstCovMtx.axes[1].mirrorOtherAxis: commonColAxis = copy.copy( firstCovMtx.axes[0] )
        else:                                   commonColAxis = copy.copy( firstCovMtx.axes[1] )
        commonMatrixAxis = copy.copy( firstCovMtx.axes[2] )
        commonType = firstCovMtx.type
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        
        # first pass through components is to collect bins to set up the common grid + do assorted checking
        for p in self.pointerList[1:]:
            cc = p.link.getNativeData().toCovarianceMatrix() # a little recursion to take care of nested covariances
            if cc.type != commonType: raise ValueError( "Incompatable types in "+str(self.__class__)+": "+str(commonType)+' vs. '+str(cc.type) )
            cc.convertAxesToUnits( ( commonRowAxis.unit, commonColAxis.unit, commonMatrixAxis.unit ) )
            commonRowAxis.data = commonRowAxis.data + cc.axes[0].data
            if cc.axes[1].mirrorOtherAxis: commonColAxis.data = commonColAxis.data + cc.axes[0].data
            else:                          commonColAxis.data = commonColAxis.data + cc.axes[1].data
        commonRowAxis.data.sort()
        commonColAxis.data.sort()
        commonRowAxis.data = uniquify( commonRowAxis.data )
        commonColAxis.data = uniquify( commonColAxis.data )
        
        # now sum up the components
        commonMatrix = self.pointerList[0]['coefficient'] * numpy.mat( firstCovMtx.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        for p in self.pointerList[1:]:
            cc = p.link.getNativeData().toCovarianceMatrix() # a little recursion to take care of nested covariances
            commonMatrix += p['coefficient'] * numpy.mat( cc.group( ( commonRowAxis.data, commonColAxis.data ), ( commonRowAxis.unit, commonColAxis.unit ) ).matrix.data )
        
        # now create the instance of the resulting covarianceMatrix
        return covarianceMatrix( type=commonType, axes=[ commonRowAxis, commonColAxis, commonMatrixAxis ], matrix=gndMatrix.matrix( commonMatrix.tolist(), form=gndMatrix.symmetricFormToken ) )
        

    def toXMLList(self, flags=None, verbosityIndent='', indent=''):
        indent2 = indent+'  '
        xmlString = [ indent+'<%s' % self.moniker ]
        if self.index!=None: xmlString[0] += ' index="%s"' % self.index
        xmlString[0] += ' lowerBound="%s" upperBound="%s">' % (
                self.lowerBound, self.upperBound)
        xmlString.append(indent2+'<!-- The matrix for this reaction equals the weighted sum of the following matrices: -->')
        for pointer in self.pointerList:
            xmlString.append( pointer.toXML( indent2 ) )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString
        
    def getReferredCovariance( self, pointer ):
        if pointer.link is not None: return pointer.link
        if 'covarianceSuite' in pointer.path:   return link.follow( pointer.path, self.getRootParent() )
        #elif'reactionSuite' in pointer.path:    return link.follow( pointer.path, None )
        else :                                  raise ValueError( "Need reference to root node of "+str( pointer.path ) )

    def getUncertaintyVector( self, theData=None, relative=True ):
        """
        Combine all subsections into single uncertainty vector, converting to relative if requested.
        
        :returns: an XYs instance 
        """
        #xys = [ self.getReferredCovariance( p ).getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        xys = [ p.link.getNativeData().getUncertaintyVector( theData=theData, relative=relative ) for p in self.pointerList ]
        coefficients = [ p['coefficient'] for p in self.pointerList ]
        uncert = coefficients[0] * xys[0]
        for i, xy in enumerate( xys[1:] ):
            uncert, xy = uncert.mutualify( 1e-8, 0, 0, xy, 1e-8, 0, 0 )
            uncert += coefficients[ i ] * xy
        return uncert

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        from fudge.gnd import link
        index = element.get("index")
        if index is not None: index=int(index)
        lower = PhysicalQuantityWithUncertainty( element.get("lowerBound") )
        upper = PhysicalQuantityWithUncertainty( element.get("upperBound") )
        summed_ = summedCovariance( index, lower, upper )
        for summand in element:
            link_ = link.parseXMLNode( summand, linkData )
            link_['coefficient'] = float( link_['coefficient'] )
            linkData['unresolvedLinks'].append( link_ )
            summed_.pointerList.append( link_ )
        return summed_

    def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
        endf = []
        if not inCovarianceGroup:
            # print header for this subsection (contains one NL sub-subsection)
            XMF1,XLFS1,MAT1,NC,NI = 0,0,0,1,0
            rowdat, coldat = targetInfo['dataPointer']
            MT1 = map(int, rowdat['ENDF_MFMT'].split(',')) [1]
            endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
        # header:
        LTY=0
        endf.append( endfFormats.endfHeadLine(0,0,0,LTY,0,0) )
        NCI = len(self.pointerList)
        endf.append( endfFormats.endfHeadLine(self.lowerBound.getValueAs('eV'),
            self.upperBound.getValueAs('eV'),0,0, 2*NCI,NCI) )
        mtList = [ map(int, a['ENDF_MFMT'].split(','))[1] for a in self.pointerList]
        coefficients = [a['coefficient'] for a in self.pointerList]
        endf += endfFormats.endfDataList( [i for j in zip(coefficients,mtList)
            for i in j] )
        return endf

def readXML( gndCovariancesFile, reactionSuite=None ):
    if reactionSuite is None:
        sys.stderr.write("WARNING: without a reactionSuite instance, covariances will have unresolved links!\n")
    from xml.etree import cElementTree
    csElement = cElementTree.parse( gndCovariancesFile ).getroot()
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode
    csElement = xmlNode( csElement, xmlNode.etree )
    linkData = {'reactionSuite': reactionSuite, 'unresolvedLinks':[]}
    return covarianceSuite.parseXMLNode( csElement, linkData )
