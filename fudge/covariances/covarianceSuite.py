# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Module with containers for covariances in several different forms
"""

import sys
import os

from xData import formatVersion as formatVersionModule
from xData.ancestry import ancestry

from fudge import styles as stylesModule, suites as suitesModule
from . import covarianceSection as sectionModule, summed as summedModule, mixed as mixedModule, modelParameters as modelParametersModule

__metaclass__ = type


class covarianceSections( suitesModule.suite ):
    """
    Most covariances are stored in 'sections', each representing either the internal covariance
    for a single quantity (i.e. cross section) or the cross-term between two quantities
    """

    moniker = 'covarianceSections'
    legacyMemberNameMapping = { 'section' : sectionModule.covarianceSection.moniker }

    def __init__(self):
        suitesModule.suite.__init__(self, (sectionModule.covarianceSection,))

class parameterCovariances( suitesModule.suite ):
    """
    Resolved and unresolved resonance parameters are stored inside the parameterCovariances
    """

    moniker = 'parameterCovariances'
    legacyMemberNameMapping = { 'section' : modelParametersModule.averageParameterCovariance.moniker }

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
    childNodeOrder = {
            formatVersionModule.version_1_10 :       (  suitesModule.externalFiles.moniker,                 stylesModule.styles.moniker,
                                                        covarianceSections.moniker,                         parameterCovariances.moniker ),
            formatVersionModule.version_2_0_LLNL_4 : (  suitesModule.externalFiles.moniker,                 stylesModule.styles.moniker,
                                                        covarianceSections.moniker,                         parameterCovariances.moniker ) }

    def __init__( self, projectile, target, evaluation, formatVersion = formatVersionModule.default, sourcePath = None ) :
        """
        :param projectile:  particle id
        :param target:      particle id
        :param evaluation:  evaluation id
        """

        ancestry.__init__( self )

        if formatVersion not in formatVersionModule.allowed: raise Exception("Unsupported GNDS structure '%s'!" % str(formatVersion))

        self.projectile = projectile            #: The projectile
        self.target = target                    #: The target
        self.evaluation = evaluation

        self.__sourcePath = sourcePath

        self.__externalFiles = suitesModule.externalFiles()
        self.__externalFiles.setAncestor( self )

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )

        self.__covarianceSections = covarianceSections()
        self.__covarianceSections.setAncestor( self )

        self.__parameterCovariances = parameterCovariances()
        self.__parameterCovariances.setAncestor( self )

        self.formatVersion = formatVersion 

    @property
    def sourcePath( self ) :
        """Returns the sourcePath member which is the path to the covarianceSuite file for self if self is from a file."""

        return( self.__sourcePath )
    
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

        dirname = os.path.dirname( fileName )
        if len(dirname) > 0 and not os.path.exists(dirname): os.makedirs(dirname)

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

        from fudge import warning

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
                        stack.pop()
            return None

        def get_edges( section_, style ):   #: return list of all pointers from this section
            natDat = section_[style]
            if isinstance(natDat, summedModule.summedCovariance): return [v.link for v in natDat]
            elif isinstance(natDat, mixedModule.mixedForm):
                edges = []
                for part in natDat:
                    if isinstance(part, summedModule.summedCovariance): edges += [v.link for v in part]
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
                warnings.append( warning.context("covarianceSection '%s':" % section_.label, sectionWarnings ) )

        for parameterCovariance in self.parameterCovariances:
            pcwarnings = parameterCovariance.check( info )
            if pcwarnings:
                warnings.append(warning.context("parameterCovariance '%s':" % parameterCovariance.label, pcwarnings))

        return warning.context('CovarianceSuite: %s + %s' % (self.projectile, self.target), warnings)

    def fix( self, **kwargs ):
        """
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

        """
        from fudge import warning
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

    def toXML( self, indent = '', **kwargs ) :
        """Returns an GNDS/XML string representation of self."""

        return( '\n'.join( self.toXMLList( **kwargs ) ) )
    
    def toXMLList( self, indent = '', **kwargs ) :
        """Returns a list of GNDS/XML strings representing self."""

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        formatVersion = kwargs.get( 'formatVersion', self.formatVersion )
        kwargs['formatVersion'] = formatVersion
        if formatVersion not in formatVersionModule.allowed: raise Exception("Unsupported GNDS structure '%s'!" % str(formatVersion))

        xmlString = [ '%s<%s projectile="%s" target="%s" evaluation="%s" format="%s">'
                % ( indent, self.moniker, self.projectile, self.target, self.evaluation, formatVersion ) ]
        xmlString += self.externalFiles.toXMLList( indent2, **kwargs )
        xmlString += self.styles.toXMLList( indent2, **kwargs )
        xmlString += self.covarianceSections.toXMLList( indent2, **kwargs )
        xmlString += self.parameterCovariances.toXMLList( indent2, **kwargs )
        xmlString.append( '%s</%s>' % (indent, self.moniker) )
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath, linkData, **kwargs ):

        xPath.append( element.tag) # keep track of location in the tree, in case errors come up

        sourcePath = kwargs.get( 'sourcePath' )

        try:
            formatVersion = element.get('format')
            linkData['formatVersion'] = formatVersion
            CS = covarianceSuite(element.get('projectile'), element.get('target'), element.get('evaluation'),
                    formatVersion, sourcePath=sourcePath)
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

def readXML( gndsCovariancesFile, reactionSuite=None, warningNoReactionSuite=True ):
    if reactionSuite is None:
        if( warningNoReactionSuite ) : sys.stderr.write("WARNING: without a reactionSuite instance, covariances will have unresolved links!\n")
    from xml.etree import cElementTree
    csElement = cElementTree.parse( gndsCovariancesFile ).getroot()
    # wrapper around the xml parser:
    from LUPY.xmlNode import xmlNode
    csElement = xmlNode( csElement, xmlNode.etree )
    linkData = {'reactionSuite': reactionSuite, 'unresolvedLinks':[]}
    covSuite = covarianceSuite.parseXMLNode( csElement, xPath=[], linkData=linkData, sourcePath = gndsCovariancesFile )
    if reactionSuite is not None:
        for externalLink in reactionSuite._externalLinks:
            externalLink.link = externalLink.follow( covSuite )
    return covSuite
