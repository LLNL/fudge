# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.
"""

import datetime

from LUPY import ancestry as ancestryModule

from .. import text as textModule
from .. import suite as suiteModule
from .. import date as dateModule

class CovarianceScript( textModule.Text ) :

    moniker = 'covarianceScript'

class CorrectionScript( textModule.Text ) :

    moniker = 'correctionScript'

class Note( textModule.Text ) :

    moniker = 'note'

class ExforDataSet(ancestryModule.AncestryIO):
    """A class representing a GNDS authors/exforDataSets/exforDataSet node."""

    moniker = 'exforDataSet'
    keyName = 'subentry'
    ancestryMembers = ( 'covarianceScript', 'correctionScript', 'note' )

    def __init__( self, _subentry, _retrievalDate ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__subentry = textModule.raiseIfNotString( _subentry, 'subentry' )
        self.__retrievalDate = dateModule.raiseIfNotDate( _retrievalDate )
        self.__covarianceScript = CovarianceScript( )
        self.__correctionScript = CorrectionScript( )
        self.__note = Note( )

    @property
    def subentry( self ) :

        return( self.__subentry )

    @property
    def covarianceScript( self ) :

        return( self.__covarianceScript )

    @property
    def correctionScript( self ) :

        return( self.__correctionScript )

    @property
    def note( self ) :

        return( self.__note )

    @property
    def retrievalDate( self ):

        return( self.__retrievalDate.__str__ )

    def toXML_strList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        XMLList  = [ '%s<%s subentry="%s" retrievalDate="%s">' % ( indent, self.moniker, self.__subentry, self.retrievalDate() ) ]
        XMLList += self.__covarianceScript.toXML_strList( indent2, **kwargs )
        XMLList += self.__correctionScript.toXML_strList( indent2, **kwargs )
        XMLList += self.__note.toXML_strList( indent2, **kwargs )        
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs ) :

        _subentry = node.get( 'subentry' )
        _retrievalDate = datetime.datetime.strptime( node.get( 'retrievalDate' ), '%Y-%m-%d' )

        return csl(_subentry, _retrievalDate)

class ExforDataSets( suiteModule.Suite ) :
    """A class representing a GNDS experimentalDataSets/exforDataSets node."""

    moniker = 'exforDataSets'
    suiteName = 'subentry'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ ExforDataSet ] )

class ExperimentalDataSets(ancestryModule.AncestryIO_base):

    moniker = 'experimentalDataSets'
    ancestryMembers = ( 'exforDataSets', )

    def __init__(self):

        ancestryModule.AncestryIO_base.__init__(self)

        self.__exforDataSets = ExforDataSets( )

    @property
    def exforDataSets( self ) :

        return ( self.__exforDataSets )

    def toXML_strList( self, indent = '', **kwargs ) :

        if not len(self.__exforDataSets): return []

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        XMLList  = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLList += self.__exforDataSets.toXML_strList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    def parseNode(self, node, xPath, linkData, **kwargs):

        self.__exforDataSets.parseNode(node, xPath, linkData, **kwargs)
