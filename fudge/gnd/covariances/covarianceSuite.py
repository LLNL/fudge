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
Module with containers for covariances in several different forms
"""

import sys
from xData.ancestry import ancestry
from fudge.gnd import styles as stylesModule
from . import section, summed, mixed, modelParameters
from ..version import GND_VERSION

__metaclass__ = type

class covarianceSuite( ancestry ):
    """
    All covariances for a target/projectile combination are stored in a :py:class:`covarianceSuite` in gnd.
    The :py:class:`covarianceSuite` is stored in a separate file from the reactionSuite.

    Within the :py:class:`covarianceSuite`, data is sorted into :py:class:`sections` (see the section module), each
    of which contains one section of the full covariance matrix.
    
    """

    moniker = 'covarianceSuite'

    def __init__(self, projectile=None, target=None, reactionSums=None,
            externalReactions=None, sections=None, modelParameterCovariances=None, style=None):

        ancestry.__init__( self )

        self.projectile = projectile            #: The projectile
        self.target = target                    #: The target
        self.reactionSums = reactionSums or []  #: List of lumped sums, etc
        self.externalReactions = externalReactions or [] #: List of cross-material covariances
        self.sections = sections or []          #: List of section instances
        self.modelParameterCovariances = modelParameterCovariances or [] #: List of modelParameterCovariance instances
        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )
        if( style is not None ) : self.styles.add( style )

        self.version = GND_VERSION
    
    def __getitem__(self, idx):
        return (self.sections+self.modelParameterCovariances)[idx]
    
    def __len__(self):
        return len(self.sections+self.modelParameterCovariances)

    @property
    def styles(self):
        return self.__styles

    def addSection(self,section_):
        self.sections.append(section_)
        self.sections[-1].setAncestor( self, 'label' )

    def addModelParameterCovariance(self,mpcovar):
        self.modelParameterCovariances.append(mpcovar)
        self.modelParameterCovariances[-1].setAncestor( self, 'label' )
    
    def addReactionSum(self, reactionSum):
        self.reactionSums.append(reactionSum)
        self.reactionSums[-1].setAncestor( self, 'id' )

    def addExternalReaction(self, externalReaction_):
        self.externalReactions.append(externalReaction_)
        self.externalReactions[-1].setAncestor( self, 'id' )

    def saveToOpenedFile( self, fOut, **kwargs ) :

        xmlString = self.toXMLList( **kwargs )
        fOut.write( '\n'.join( xmlString ) )
        fOut.write( '\n' )

    def saveToFile( self, fileName, **kwargs ):

        with open(fileName,"w") as fout:
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n' )
            self.saveToOpenedFile( fout, **kwargs )
    
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
        info = { 'covarianceSuite': self, 'style': 'eval' }
        info.update( options )

        warnings = []

        #: summedReactions are sums of other sections, but that opens the possibility of a cyclic dependency
        #: check for that, using algorithm taken from
        #:    http://neopythonic.blogspot.com/2009/01/detecting-cycles-in-directed-graph.html 
        def find_cycle(NODES, EDGES, style):
            todo = set(NODES)
            while todo:
                node = todo.pop()
                stack = [node]
                while stack:
                    top = stack[-1]
                    for node in EDGES(top, style):
                        if node in stack:
                            return stack[stack.index(node):]
                        if node in todo:
                            stack.append(node)
                            todo.remove(node)
                            break
                    else:  # this node is not in a cycle
                        node = stack.pop()
            return None

        def get_edges( section_, style ):   #: return list of all pointers from this section
            natDat = section_.forms[style]
            if isinstance(natDat, summed.summedCovariance): return [v.link for v in natDat.pointerList]
            elif isinstance(natDat, mixed.mixedForm):
                edges = []
                for part in natDat:
                    if isinstance(part, summed.summedCovariance): edges += [v.link for v in part.pointerList]
                return edges
        nodes = [sec for sec in self.sections if get_edges( sec, info['style'] )]  # sections that contain pointers

        if nodes:
            cycle = find_cycle(nodes, get_edges, info['style'])
            if cycle:
                warnings.append( warning.cyclicDependency( cycle ) )

        # check each section
        for section_ in self.sections:
            sectionWarnings = section_.check( info )
            if sectionWarnings:
                warnings.append( warning.context('Section %s (%s):' % (section_.label, section_.id), sectionWarnings ) )

        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides ancestry.findEntity. covarianceSuite contains more than one list,
        so may need to descend into those to find desired entity
        """
        if entityName in ('section','reactionSum','externalReaction'):
            for entity in getattr( self, entityName+'s' ):
                if getattr( entity, attribute, None ) == value:
                    return entity
        return ancestry.findEntity( self, entityName, attribute, value )

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
            'negativeEigenTolerance': -1e-6,  # ignore smaller negative eigenvalues
            'eigenvalueRatioTolerance': 1e-8, # warn if smallest eival < 1e-8 * biggest
            'eigenvalueAbsoluteTolerance': 1e-14,
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
    
    def toXMLList( self, indent = '', **kwargs ) :
        """Write self out to GND-XML"""

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xmlString = [ '%s<%s projectile="%s" target="%s" version="%s"'
                % ( indent, self.moniker, self.projectile, self.target, self.version ) ]
        xmlString[-1] += ' xmlns:xlink="http://www.w3.org/1999/xlink">'
        xmlString += self.styles.toXMLList( indent2, **kwargs )
        if self.reactionSums:
            xmlString += [ indent2 + '<reactionSums>',
                indent3 + '<!-- Covariances may be given for a sum of several reactions. Define these "summed reactions" here: -->' ]
            for reactionSum in self.reactionSums:
                xmlString += reactionSum.toXMLList( indent3, **kwargs )
            xmlString[-1] += '</reactionSums>'
        if self.externalReactions:
            xmlString += [ indent2 + '<externalReactions>',
                indent3 + "<!-- This target has covariances with reactions on different target(s). List other target/reactions: -->" ]
            for externalReaction_ in self.externalReactions:
                xmlString += externalReaction_.toXMLList( indent3, **kwargs )
            xmlString[-1] += '</externalReactions>'
        for section_ in self.sections + self.modelParameterCovariances:
            xmlString += section_.toXMLList( indent2, **kwargs )
        xmlString.append( indent + '</%s>' % self.moniker )
        return xmlString

    @staticmethod
    def parseXMLNode( CSelement, xPath, linkData ):
        '''Parse an XML node and load create a covarianceSuite based on its contents'''

        xPath.append( CSelement.tag ) # keep track of location in the tree, in case errors come up
        try:
            covariances = covarianceSuite( CSelement.get('projectile'), CSelement.get('target') )
            for child in CSelement:
                if child.tag == 'styles':
                    covariances.styles.parseXMLNode( child, xPath, linkData )
                elif child.tag == 'reactionSums':
                    for reactionSumElement in child:
                        covariances.addReactionSum( section.reactionSum.parseXMLNode(
                            reactionSumElement, xPath, linkData ) )
                elif child.tag == 'externalReactions':
                    for external in child:
                        covariances.addExternalReaction( section.externalReaction.parseXMLNode(
                            external, xPath, linkData ) )
                elif child.tag == section.section.moniker:
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
            res = link_.follow(root)
            link_.link = res

        xPath.pop()
        return covariances

def readXML( gndCovariancesFile, reactionSuite=None ):
    if reactionSuite is None:
        sys.stderr.write("WARNING: without a reactionSuite instance, covariances will have unresolved links!\n")
    from xml.etree import cElementTree
    csElement = cElementTree.parse( gndCovariancesFile ).getroot()
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode
    csElement = xmlNode( csElement, xmlNode.etree )
    linkData = {'reactionSuite': reactionSuite, 'unresolvedLinks':[]}
    return covarianceSuite.parseXMLNode( csElement, xPath=[], linkData=linkData )
