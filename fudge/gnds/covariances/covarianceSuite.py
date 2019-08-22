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

"""
Module with containers for covariances in several different forms
"""

import sys
from xData.ancestry import ancestry
from fudge.gnds import styles as stylesModule, suites as suitesModule
from . import section as sectionModule, summed as summedModule, mixed as mixedModule, modelParameters as modelParametersModule
from ..version import GNDS_VERSION

__metaclass__ = type


class covarianceSections( suitesModule.suite ):
    """
    Most covariances are stored in 'sections', each representing either the internal covariance
    for a single quantity (i.e. cross section) or the cross-term between two quantities
    """

    moniker = 'covarianceSections'

    def __init__(self):
        suitesModule.suite.__init__( self, (sectionModule.section,) )


class parameterCovariances( suitesModule.suite ):
    """
    Resolved and unresolved resonance parameters are stored inside the parameterCovariances
    """

    moniker = 'parameterCovariances'

    def __init__( self ):
        suitesModule.suite.__init__(self, (modelParametersModule.parameterCovariance,
                                           modelParametersModule.averageParameterCovariance) )


class covarianceSuite( ancestry ):
    """
    All covariances for a target/projectile combination are stored in a :py:class:`covarianceSuite` in gnds.
    The :py:class:`covarianceSuite` is stored in a separate file from the reactionSuite.

    Within the :py:class:`covarianceSuite`, data is sorted into :py:class:`sections` (see the section module), each
    of which contains one section of the full covariance matrix.
    
    """

    moniker = 'covarianceSuite'

    def __init__(self, projectile, target, evaluation, GNDS_version = GNDS_VERSION):
        """
        :param projectile:  particle id
        :param target:      particle id
        :param evaluation:  evaluation id
        """

        ancestry.__init__( self )

        if GNDS_version not in (GNDS_VERSION,): raise Exception("Unsupported GNDS structure '%s'!" % str(GNDS_version))

        self.projectile = projectile            #: The projectile
        self.target = target                    #: The target
        self.evaluation = evaluation

        self.__externalFiles = suitesModule.externalFiles()
        self.__externalFiles.setAncestor( self )

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )

        self.__covarianceSections = covarianceSections()
        self.__covarianceSections.setAncestor( self )

        self.__parameterCovariances = parameterCovariances()
        self.__parameterCovariances.setAncestor( self )

        self.GNDS_version = GNDS_version
    
    @property
    def styles(self):
        return self.__styles

    @property
    def externalFiles(self):
        return self.__externalFiles

    @property
    def covarianceSections(self):
        return self.__covarianceSections

    @property
    def parameterCovariances(self):
        return self.__parameterCovariances

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        for section in self.covarianceSections: section.convertUnits( unitMap )
        for mpsection in self.parameterCovariances: mpsection.convertUnits( unitMap )

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

        from fudge.gnds import warning

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
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError( "check() received unknown keyword argument '%s'" % key )

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
            natDat = section_[style]
            if isinstance(natDat, summedModule.summedCovariance): return [v.link for v in natDat.pointerList]
            elif isinstance(natDat, mixedModule.mixedForm):
                edges = []
                for part in natDat:
                    if isinstance(part, summedModule.summedCovariance): edges += [v.link for v in part.pointerList]
                return edges
        nodes = [sec for sec in self.covarianceSections if get_edges( sec, info['style'] )]  # sections that contain pointers

        if nodes:
            cycle = find_cycle(nodes, get_edges, info['style'])
            if cycle:
                warnings.append( warning.cyclicDependency( cycle ) )

        # check each section
        for section_ in self.covarianceSections:
            sectionWarnings = section_.check( info )
            if sectionWarnings:
                warnings.append( warning.context("Section '%s':" % section_.label, sectionWarnings ) )

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
        from fudge.gnds import warning
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
            if key in options:
                options[key] = kwargs[key]
            else:
                raise KeyError( "fix() received unknown keyword argument '%s'" % key )
        # assemble some useful info, to be handed down to children's check() functions:
        info = { 'covarianceSuite': self }
        info.update( options )
        # do the fixing
        warnings = []
        for section_ in self.covarianceSections: warnings += section_.fix( **info )
        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)
    
    def removeExtraZeros(self):
        """
        Checks all covariance matrices for rows/columns of all zero, removes them if found.
        """
        for section_ in self.covarianceSections:
            for form in section_:
                if hasattr(form, "removeExtraZeros"):
                    form.removeExtraZeros()
                if isinstance(form, mixedModule.mixedForm):
                    for covar in form:
                        if hasattr(covar, "removeExtraZeros"):
                            covar.removeExtraZeros()
    
    def toXMLList( self, indent = '', **kwargs ) :
        """Write self out to GNDS-XML"""

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent

        xmlString = [ '%s<%s projectile="%s" target="%s" evaluation="%s" format="%s">'
                % ( indent, self.moniker, self.projectile, self.target, self.evaluation, self.GNDS_version ) ]
        xmlString += self.styles.toXMLList( indent2, **kwargs )
        xmlString += self.externalFiles.toXMLList( indent2, **kwargs )
        xmlString += self.covarianceSections.toXMLList( indent2, **kwargs )
        xmlString += self.parameterCovariances.toXMLList( indent2, **kwargs )
        xmlString.append( '%s</%s>' % (indent, self.moniker) )
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag) # keep track of location in the tree, in case errors come up
        try:
            CS = covarianceSuite(element.get('projectile'), element.get('target'), element.get('evaluation'),
                    element.get('format'))
            for child in element:
                if child.tag == 'styles':
                    CS.styles.parseXMLNode( child, xPath, linkData )
                elif child.tag == suitesModule.externalFiles.moniker:
                    CS.externalFiles.parseXMLNode( child, xPath, linkData )
                elif child.tag == covarianceSections.moniker:
                    CS.covarianceSections.parseXMLNode( child, xPath, linkData )
                elif child.tag == parameterCovariances.moniker:
                    CS.parameterCovariances.parseXMLNode( child, xPath, linkData )
                else:
                    print("Warning: encountered unexpected element '%s' in covarianceSuite!" % child.tag)
        except Exception:
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise

        # fix links:
        for link_ in linkData['unresolvedLinks']:
            path = link_.path
            if link_.root is None:
                if path.startswith('/covarianceSuite'):
                    root = CS
                elif path.startswith('.'):
                    root = link_    # relative link
            elif link_.root == '$reactions':
                root = linkData['reactionSuite']
            else:
                root = None

            if root is None:
                continue    # FIXME links to other files (i.e. cross-material covariances) currently broken
            res = link_.follow(root)
            link_.link = res

        xPath.pop()
        return CS

def readXML( gndsCovariancesFile, reactionSuite=None ):
    if reactionSuite is None:
        sys.stderr.write("WARNING: without a reactionSuite instance, covariances will have unresolved links!\n")
    from xml.etree import cElementTree
    csElement = cElementTree.parse( gndsCovariancesFile ).getroot()
    # wrapper around the xml parser:
    from fudge.core.utilities.xmlNode import xmlNode
    csElement = xmlNode( csElement, xmlNode.etree )
    linkData = {'reactionSuite': reactionSuite, 'unresolvedLinks':[]}
    covSuite = covarianceSuite.parseXMLNode( csElement, xPath=[], linkData=linkData )
    if reactionSuite is not None:
        for externalLink in reactionSuite._externalLinks:
            externalLink.link = externalLink.follow( covSuite )
    return covSuite
