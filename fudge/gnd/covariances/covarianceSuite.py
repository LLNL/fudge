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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
Module with containers for covariances in several different forms
"""

import sys
import fudge
from fudge.gnd import link
from fudge.core.ancestry import ancestry
from fudge.legacy.converting import endfFormats
from . import tokens, section, summed, mixed, modelParameters

__metaclass__ = type

class covarianceSuite( ancestry ):
    """
    All covariances for a target/projectile combination are stored in a :py:class:`covarianceSuite` in gnd.
    The :py:class:`covarianceSuite` is stored in a separate file from the reactionSuite.

    Within the :py:class:`covarianceSuite`, data is sorted into :py:class:`sections` (see the section module), each
    of which contains one section of the full covariance matrix.
    
    """
    def __init__(self, projectile=None, target=None, reactionSums=None,
            externalReactions=None, sections=None, modelParameterCovariances=None, styles=None):

        ancestry.__init__( self, tokens.covarianceSuiteToken, None )

        self.projectile = projectile            #: The projectile
        self.target = target                    #: The target
        self.reactionSums = reactionSums or []  #: List of lumped sums, etc
        self.externalReactions = externalReactions or [] #: List of cross-material covariances
        self.sections = sections or []          #: List of section instances
        self.modelParameterCovariances = modelParameterCovariances or [] #: List of modelParameterCovariance instances
        self.styles = styles or {}              #: Python dict of styles
        self.version = "gnd version 1.0"        #: FIXME should be taken from fudge.gnd.version
    
    def __getitem__(self, idx):
        return (self.sections+self.modelParameterCovariances)[idx]
    
    def __len__(self):
        return len(self.sections+self.modelParameterCovariances)

    def addSection(self,section_):
        section_.setParent( self )
        self.sections.append(section_)

    def addModelParameterCovariance(self,mpcovar):
        mpcovar.setParent( self )
        self.modelParameterCovariances.append(mpcovar)
    
    def addReactionSum(self, reactionSum):
        reactionSum.setParent( self )
        self.reactionSums.append(reactionSum)

    def addExternalReaction(self, externalReaction_):
        externalReaction_.setParent( self )
        self.externalReactions.append(externalReaction_)

    def saveToOpenedFile( self, fOut, flags=None, verbosityIndent='' ):
        xmlString = self.toXMLList( flags = flags )
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

        def get_edges( section_ ):   #: return list of all pointers from this section
            natDat = section_.getNativeData()
            if isinstance(natDat, summed.summedCovariance): return [v.link for v in natDat.pointerList]
            elif isinstance(natDat, mixed.mixedForm):
                edges = []
                for part in natDat:
                    if isinstance(part, summed.summedCovariance): edges += [v.link for v in part.pointerList]
                return edges
        nodes = [sec for sec in self.sections if get_edges( sec )]  # sections that contain pointers

        if nodes:
            cycle = find_cycle(nodes, get_edges)
            if cycle:
                warnings.append( warning.cyclicDependency( cycle ) )

        # check each section
        for section_ in self.sections:
            sectionWarnings = section_.check( info )
            if sectionWarnings:
                warnings.append( warning.context('Section: %s' % section_.id, sectionWarnings ) )

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
        for section_ in self.sections: warnings += section_.fix( **info )        
        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)
    
    def removeExtraZeros(self):
        """ remove columns/rows of zero from all covariances """
        for section_ in self.sections:
            for form in section_.forms.values():
                if hasattr(form, "removeExtraZeros"):
                    form.removeExtraZeros()
                if isinstance(form, mixed.mixedForm):
                    for covar in form.subsections:
                        if hasattr(covar, "removeExtraZeros"):
                            covar.removeExtraZeros()
    
    def toXMLList(self, flags=None, indent=''):
        """Write self out to GND-XML"""
        indent2 = indent+'  '; indent3 = indent2+'  '
        if flags is None: flags = {}
        xmlString = ['<?xml version="1.0" encoding="UTF-8"?>']
        xmlString.append(indent+'<%s projectile="%s" target="%s" version="%s"'
                % (tokens.covarianceSuiteToken, self.projectile, self.target, self.version))
        xmlString[-1] += ' xmlns:xlink="http://www.w3.org/1999/xlink">'
        xmlString.append('%s<styles>' % (indent+'  '))
        for style in self.styles: xmlString += self.styles[style].toXMLList( indent=indent+'    ' )
        xmlString[-1] += '</styles>'
        if self.reactionSums:
            xmlString += [indent2+'<reactionSums>',
                indent3+'<!-- Covariances may be given for a sum of several reactions. Define these "summed reactions" here: -->']
            for reactionSum in self.reactionSums:
                xmlString += reactionSum.toXMLList(flags, indent=indent3)
            xmlString[-1] += '</reactionSums>'
        if self.externalReactions:
            xmlString += [indent2+'<externalReactions>',
                indent3+"<!-- This target has covariances with reactions on different target(s). List other target/reactions: -->"]
            for externalReaction_ in self.externalReactions:
                xmlString += externalReaction_.toXMLList(flags, indent=indent3)
            xmlString[-1] += '</externalReactions>'
        for section_ in self.sections + self.modelParameterCovariances:
            xmlString += section_.toXMLList(flags, indent+'  ')
        xmlString.append(indent+'</%s>' % tokens.covarianceSuiteToken)
        return xmlString

    @staticmethod
    def parseXMLNode( CSelement, linkData={} ):
        '''Parse an XML node and load create a covarianceSuite based on its contents'''

        xPath = ['covarianceSuite'] # keep track of location in the tree, in case errors come up
        try:
            covariances = covarianceSuite( CSelement.get('projectile'), CSelement.get('target') )
            for child in CSelement:
                if child.tag == 'styles':
                    styles = [fudge.gnd.miscellaneous.style.parseXMLNode( style, xPath ) for style in child]
                    for style in styles: covariances.styles[ style.name ] = style
                elif child.tag == 'reactionSums':
                    for reactionSumElement in child:
                        covariances.addReactionSum( section.reactionSum.parseXMLNode(
                            reactionSumElement, xPath, linkData ) )
                elif child.tag == 'externalReactions':
                    for external in child:
                        covariances.addExternalReaction( section.externalReaction.parseXMLNode(
                            external, xPath, linkData ) )
                elif child.tag == tokens.sectionToken:
                    covariances.addSection( section.section.parseXMLNode( child, xPath, linkData ) )
                elif child.tag in ('resonanceParameterCovariance',):
                    covariances.addModelParameterCovariance( modelParameters.resonanceParameterCovariance.parseXMLNode(
                        child, xPath, linkData ) )
                else:
                    print("Warning: encountered unexpected element '%s' in covarianceSuite!" % child.tag)
        except Exception:
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise

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
        from fudge.legacy.converting import endfFormats
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        NIS, ABN = 1,1.0; ZAI=ZAM  # assuming one isotope/file
        MTL = 0 # mtl=1 sections are handled in lumpedCovariance

        for section_ in self.modelParameterCovariances:
            section_.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)

        sections = self.sections[:]
        # sort covariances by MF/MT:
        mfmts = []
        for section_ in sections:
            mfmts.append(map(int,section_.rowData['ENDF_MFMT'].split(',')))
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
                if thisMFMT[0].nativeData == tokens.legendreOrderCovarianceFormToken: LTT = 1
                NMT1 = len(thisMFMT)
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, LTT, 0, NMT1 )]
            elif mf==35:
                NK = 1
                if thisMFMT[0].nativeData == tokens.energyIntervalFormToken:
                    NK = len(thisMFMT[0].forms[tokens.energyIntervalFormToken])
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, MTL, NK, 0 )]
            elif mf==40:
                endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, 0, len(thisMFMT), 0 )]
            for section_ in thisMFMT:
                MAT1 = 0
                if section_.columnData and isinstance(section_.columnData.link, section.externalReaction):
                    otherTarget = section_.columnData.link.target
                    from fudge.legacy.converting import endf_endl
                    ZA, MAT1 = endf_endl.ZAAndMATFromParticleName( otherTarget )
                if mf==34:
                    form = section_.forms[tokens.legendreOrderCovarianceFormToken]
                    L1s = [subsec.L1 for subsec in form]
                    L2s = [subsec.L2 for subsec in form]
                    NL = len( set(L1s) );  NL1 = len( set(L2s) )
                    if section_.columnData: raise NotImplemented # cross-reaction or cross-material
                    MT1 = mt
                    endf += [ endfFormats.endfHeadLine( 0.0, 0.0, MAT1, MT1, NL, NL1 ) ]
                if mf==40:
                    rowData = section_.rowData
                    if type(rowData) is str: raise Exception("Don't string me along!") # FIXME
                    else:
                        quant = rowData.link
                        from fudge import gnd
                        if isinstance(quant, section.reactionSum):
                            quant = quant.reactions[0].link
                        reaction = quant.parent
                        product = reaction.outputChannel[0]
                        QI = reaction.getQ('eV')
                        level = product.getLevelAsFloat( 'eV' )
                        QM = QI + level
                        LFS = product.particle.getLevelIndex()
                        NL = 1
                        endf += [endfFormats.endfHeadLine( QM, QI, 0, LFS, 0, NL ) ]
                        XMF1, XLFS1, NC, NI = 10,LFS, 0,1
                        endf += [endfFormats.endfHeadLine( XMF1,XLFS1,MAT1,mt,NC,NI )]
                form = section_.forms[section_.nativeData]
                targetInfo['MAT1'] = MAT1
                targetInfo['dataPointer'] = [section_.rowData,section_.columnData]
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
