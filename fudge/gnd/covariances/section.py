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

"""A covarianceSuite is organized into sections, where each section contains either
    - a covariance matrix for a single reaction quantity (cross section, multiplicity, etc), or
    - a covariance matrix between two different quantities (off-diagonal block)
    """
from fudge.core.ancestry import ancestry
from fudge.gnd import link
from . import tokens, covarianceMatrix
from .mixed import mixedForm
from .summed import summedCovariance
from .distributions import energyIntervalForm, LegendreOrderCovarianceForm

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

    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<%s label="%s" id="%s" nativeData="%s"' % (tokens.sectionToken, self.label,
            self.id, self.nativeData)]
        if self.columnData is not None: xmlString[0] += ' crossTerm="true"'
        xmlString[0] += '>'
        for dataPointer in ('rowData','columnData'):
            if getattr(self, dataPointer) is not None:
                xmlString.append( getattr(self, dataPointer).toXML(indent+'  ') )
        for key in self.forms.keys():
            xmlString += self.forms[key].toXMLList(flags, indent+'  ')
        xmlString[-1] += '</%s>' % tokens.sectionToken
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate <section> element from xml."""

        xPath.append( '%s[@id="%s"]' % (element.tag, element.get('id') ) )
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
                    tokens.covarianceFormToken: covarianceMatrix,
                    tokens.mixedFormToken: mixedForm,
                    tokens.summedFormToken: summedCovariance,
                    tokens.energyIntervalFormToken: energyIntervalForm,
                    tokens.legendreOrderCovarianceFormToken: LegendreOrderCovarianceForm,
                    }.get( form.tag )
            if formClass is None:
                raise Exception("encountered unknown covariance matrix form '%s'" % form.tag)
            section_.addForm( formClass.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return section_

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

    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<reactionSum id="%s" ENDF_MFMT="%i,%i">'%(self.id,self.ENDF_MFMT[0],self.ENDF_MFMT[1])]
        for ch in self.reactions:
            xmlString.append( ch.toXML( indent+'  ' ) )
        xmlString[-1] += '</reactionSum>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate <reactionSum> element from xml."""

        xPath.append( '%s[@id="%s"]' % (element.tag, element.get('id') )  )
        from fudge.gnd import link
        rsum = reactionSum( **dict(element.items()) )
        rsum.ENDF_MFMT = map(int, rsum.ENDF_MFMT.split(','))
        for child in element:
            link_ = link.parseXMLNode( child, linkData )
            linkData['unresolvedLinks'].append( link_ )
            rsum.reactions.append( link_ )
        xPath.pop()
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

    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<externalReaction id="%s" target="%s" ENDF_MFMT="%i,%i"/>'%(self.id,
            self.target, self.ENDF_MFMT[0], self.ENDF_MFMT[1])]
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate <externalReaction> element from xml."""

        xPath.append( '%s[@id="%s"]' % (element.tag, element.get('id') ) )
        exReac = externalReaction( **dict(element.items()) )
        exReac.ENDF_MFMT = map(int, exReac.ENDF_MFMT.split(','))
        xPath.pop()
        return exReac
