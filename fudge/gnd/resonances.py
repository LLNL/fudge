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
This module contains the resonances, resolved and unresolved classes.

class hierarchy:

    * resonances contains one or more of scatteringRadius, resolved and unresolved

    * resolved and unresolved contain a list of energy regions. In each region,

        - resolved may have SLBW, MLBW, RM, RML formats
        - unresolved has L-dependant or E-dependant format

      each of these subsections has a list of resonances and attributes
"""

import math
from fudge.core.ancestry import ancestry
from pqu import physicalQuantityWithUncertainty
from fudge.core.math.xData import axes, XYs
from fudge.legacy.converting import endfFormats
from fudge.core.math import table as gndTable
from fudge.gnd import xParticle

__metaclass__ = type

PQU = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty

class resonances(ancestry):
    """ 
    top-level class for resonances.
    contains scatteringRadius and (usually) resolved and/or unresolved sections
    
    bool reconstructCrossSection: tells whether we still need 
    to do the reconstruction
    """
    def __init__(self, scatteringRadius=None, resolved=None, unresolved=None, 
            reconstructCrossSection=True):
        ancestry.__init__(self,'resonances',None)
        if scatteringRadius and (resolved or unresolved):
            raise TypeError, ("Resonances should contain scattering radius OR resonance parameters, not both")
        self.scatteringRadius = scatteringRadius
        self.resolved = resolved
        self.unresolved = unresolved
        for child in (self.scatteringRadius, self.resolved, self.unresolved):
            if child is not None: child.parent = self
        self.reconstructCrossSection = reconstructCrossSection
        
    def  __str__( self ) :
        """ string representation """
        return( self.toString( simpleString = False ) )

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for section in ('scatteringRadius','resolved','unresolved'):
            if getattr(self, section) is not None:
                warningList = getattr(self,section).check(info)
                if warningList:
                    warnings.append( warning.context( section, warningList ) )
        return warnings
    
    def toXMLList( self, flags, indent = '' ) :
    
        xmlString = [indent+'<resonances reconstructCrossSection="%s">' % 
                str(self.reconstructCrossSection).lower()]
        if self.scatteringRadius:
            xmlString += self.scatteringRadius.toXMLList( flags, indent = indent + '  ', writeConstant=True )
        if self.resolved:
            xmlString += self.resolved.toXMLList( flags, indent = indent + '  ' )
        if self.unresolved:
            xmlString += self.unresolved.toXMLList( flags, indent = indent + '  ' )
        xmlString[-1] += '</resonances>'
        return xmlString

    def getDomain( self, unitTo = None, asPQU = False ):
        bounds = []
        if self.scatteringRadius:
            bounds.append( (self.scatteringRadius.lowerBound, self.scatteringRadius.upperBound) )
        if self.resolved:
            if self.resolved.multipleRegions:
                bounds += [(reg.lowerBound, reg.upperBound) for reg in self.resolved.regions]
            else: bounds.append( (self.resolved.lowerBound, self.resolved.upperBound) )
        if self.unresolved:
            bounds.append( (self.unresolved.lowerBound, self.unresolved.upperBound) )
        for i in range(len(bounds)-1):
            assert bounds[i][1] == bounds[i+1][0], "Resonance region boundaries don't match!"
        if( asPQU ) : return (bounds[0][0], bounds[-1][1])
        else: return (bounds[0][0].value, bounds[-1][1].value)

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        xPath.append( element.tag )

        scatRadius, RRR, URR = None,None,None
        kReconstruct = (element.get("reconstructCrossSection") == "true")
        for child in element:
            if child.tag=='scatteringRadius': scatRadius = scatteringRadius.parseXMLNode( child, xPath )
            elif child.tag=='resolved': RRR = resolved.parseXMLNode( child, xPath )
            elif child.tag=='unresolved': URR = unresolved.parseXMLNode( child, xPath )
            else: raise Exception("unknown element '%s' encountered in resonances!" % child.tag)
        res = resonances( scatteringRadius = scatRadius, resolved=RRR, unresolved=URR,
                reconstructCrossSection=kReconstruct )
        xPath.pop()
        return res
    
    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = ''):
        """ """
        ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
        NIS, ABN = 1, 1.0; ZAI=ZAM  # assuming only one isotope per file
        
        # get target spin from the particle list:
        target = targetInfo['reactionSuite'].getParticle( targetInfo['reactionSuite'].target.getName() )
        targetInfo['spin'] = target.getSpin().value
        
        endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0)]
        resolvedCount, unresolvedCount = 0, 0
        # resolved may have multiple energy regions:
        if self.resolved is not None: resolvedCount = max(1,len(self.resolved.regions))
        if self.unresolved is not None: unresolvedCount = 1
        
        # resonances may only contain a scattering radius:
        if not (resolvedCount + unresolvedCount) and self.scatteringRadius:
            scatRadius = self.scatteringRadius
            endf.append( endfFormats.endfHeadLine( ZAM, ABN, 0,0,1,0 ) )
            endf.append( endfFormats.endfHeadLine(scatRadius.lowerBound.getValueAs('eV'),
                scatRadius.upperBound.getValueAs('eV'), 0,0,0,0 ) )
            AP = scatRadius.value.getValueAs('10*fm')
            endf.append( endfFormats.endfHeadLine(targetInfo['spin'], AP, 0,0,0,0 ) )
            endf.append( endfFormats.endfSENDLineNumber() )
            endfMFList[2][151] = endf
            return

        # For now I'm storing the LRF/LFW flags in xml, since they are tricky to compute
        # LFW is a pain: only applies to unresolved, but must be written at the start of MF2
        LRFurr, LFW = 0,0
        if unresolvedCount != 0:
            LRF_LFW = self.unresolved.nativeData.ENDFconversionFlag
            LRFurr, LFW = map(int, LRF_LFW.split('=')[-1].split(',') )
        NER = resolvedCount + unresolvedCount
        endf.append( endfFormats.endfHeadLine( ZAI, ABN, 0, LFW, NER, 0 ) )
        for idx in range(resolvedCount):
            if resolvedCount==1: region = self.resolved
            else: region = self.resolved.regions[idx]
            LRU = 1 #resolved
            LRF = {'SingleLevel_BreitWigner':1,'MultiLevel_BreitWigner':2,'Reich_Moore':3,
                    'R_Matrix_Limited':7}[region.nativeData.moniker]
            EL, EH = region.lowerBound.getValueAs('eV'), region.upperBound.getValueAs('eV')
            NRO = region.nativeData.scatteringRadius.isEnergyDependent()
            NAPS = not( region.nativeData.calculateChannelRadius )
            endf.append(endfFormats.endfHeadLine( EL,EH,LRU,LRF,NRO,NAPS ) )
            endf += region.nativeData.toENDF6( flags, targetInfo, verbosityIndent )
        if unresolvedCount != 0:
            LRU = 2 #unresolved
            region = self.unresolved
            EL, EH = region.lowerBound.getValueAs('eV'), region.upperBound.getValueAs('eV')
            NRO, NAPS = 0,0
            endf.append(endfFormats.endfHeadLine( EL,EH,LRU,LRFurr,NRO,NAPS ) )
            # pass LFW/LRF so we don't have to calculate twice:
            targetInfo['unresolved_LFW'] = LFW
            targetInfo['unresolved_LRF'] = LRFurr
            targetInfo['regionEnergyBounds'] = (region.lowerBound, region.upperBound)
            endf += region.nativeData.toENDF6( flags, targetInfo, verbosityIndent )
        endf.append( endfFormats.endfSENDLineNumber() )
        endfMFList[2][151] = endf


    def toString( self, simpleString = False ) :
        """Returns a string representation of self. If simpleString is True, 
        the string contains only an overview without listing resonances"""
        s = 'Resonances:\n'
        if self.scatteringRadius:
            s += self.scatteringRadius.toString( simpleString = simpleString )
        if self.resolved:
            s += self.resolved.toString( simpleString = simpleString )
        if self.unresolved:
            s += self.unresolved.toString( simpleString = simpleString )
        return( s )


# scatteringRadius can be constant or energy-dependent:
class scatteringRadius(ancestry):

    moniker = 'scatteringRadius'

    def __init__( self, value=None, lowerBound=None, upperBound=None ) :

        ancestry.__init__(self, self.moniker, None)
        self.value = value
        self.lowerBound = lowerBound
        self.upperBound = upperBound

    def __eq__(self, other):
        if not isinstance(other, scatteringRadius): return False
        return( self.value==other.value )

    def __str__(self): return str( self.getXMLAttribute() )

    def __nonzero__(self): return bool(self.value)

    def check( self, info ):
        return []

    def isEnergyDependent(self):
        return isinstance( self.value, XYs.XYs )

    def getValueAs( self, unit, energy_grid=None ):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise Exception("Missing: energy grid to evaluate E-dependent scattering radius")
            energy_unit = self.value.axes[0].unit
            return [self.value.getValue_units( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( e, energy_unit )
                ).getValueAs( unit ) for e in energy_grid]
        else:
            return self.value.getValueAs( unit )

    def toString(self, simpleString = False):
        return ("Scattering radius: %s\n" % self.getXMLAttribute())

    def getXMLAttribute( self ):
        if isinstance(self.value, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty): return self.value
        else: return 'energyDependent'

    def toXMLList(self, flags, indent = '', writeConstant=False ):

        xmlString = []
        if( self.isEnergyDependent( ) ) :
            xmlString = self.value.toXMLList( tag = self.moniker, indent = indent )
        elif( writeConstant ) :       # for when the resonances section *only* contains a scattering radius
            xmlString = [ indent + '<%s value="%s" lowerBound="%s" upperBound="%s"/>' %
                    ( self.moniker, self.value, self.lowerBound, self.upperBound ) ]
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        xPath.append( element.tag )
        lower, upper = None, None
        if len(element)==0:
            value = PQU(element.get('value'))
            lower = PQU(element.get('lowerBound')); upper = PQU(element.get('upperBound'))
        else:
            value = XYs.XYs.parseXMLNode( element )
        SR = scatteringRadius( value, lower, upper )
        xPath.pop()
        return SR


# resolved/unresolved regions:
class resolved(ancestry):
    """ class for resolved resonances """
    def __init__(self, nativeData=None, lowerBound='', upperBound='', multipleRegions=False):
        ancestry.__init__(self,'resolved',None)
        if not multipleRegions:
            setattr(self, nativeData.moniker, nativeData)
        self.nativeData = nativeData
        if self.nativeData is not None: self.nativeData.parent = self
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.multipleRegions = multipleRegions
        self.regions = []  # contains multiple energyInterval (deprecated)
    
    def toString(self, simpleString = False):
        if self.regions:
            return ("Resolved region with DEPRECATED multiple regions\n")
        else:
            return ("Resolved resonances in %s form\n" % self.nativeData.moniker )

    def check( self, info ):
        import warning
        warnings = []
        if self.multipleRegions:
            warnings.append( warning.RRmultipleRegions() )
        return warnings
    
    def toXMLList( self, flags, indent = '' ) :
        if self.multipleRegions:
            xmlString = [indent+'<resolved multipleRegions="true">']
            for region in self.regions:
                xmlString += region.toXMLList( flags, indent = indent + '  ' )
        else:
            xmlString = [indent+'<resolved lowerBound="%s" upperBound="%s" nativeData="%s">' % (
                    self.lowerBound, self.upperBound, self.nativeData.moniker ) ]
            xmlString += self.nativeData.toXMLList( flags, indent + '  ')
        xmlString[-1] += '</resolved>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        xPath.append( element.tag )
        def readResolved( child ):
            formClass = {'SingleLevel_BreitWigner': SLBW,
                    'MultiLevel_BreitWigner': MLBW,
                    'Reich_Moore': RM,
                    'R_Matrix_Limited': RMatrix
                    }.get( child.tag )
            if formClass is None: raise Exception("unknown resolved resonance form '%s'!" % child.tag)
            return formClass.parseXMLNode( child )

        regions = getBool( element.get('multipleRegions','false') )
        if regions:
            regions = []
            for region in element:
                thisRegion = energyInterval( **getAttrs( region, exclude=('multipleRegions',) ) )
                thisRegion.nativeData = readResolved( region[0] )
                regions.append( thisRegion )
            RRR = resolved( multipleRegions = True )
            RRR.regions = regions
        else:
            nativeData = readResolved(element[0])
            RRR = resolved( nativeData, **getAttrs( element, exclude=('nativeData',) ) )
        xPath.pop()
        return RRR


class unresolved(ancestry):
    def __init__(self, nativeData=None, lowerBound='', upperBound=''):
        ancestry.__init__(self,'unresolved',None)
        setattr(self, nativeData.moniker, nativeData)
        self.nativeData = nativeData
        if isinstance( nativeData, ancestry ): self.nativeData.parent = self
        self.lowerBound = lowerBound
        self.upperBound = upperBound
    
    def toString( self, simpleString = False ):
        return ("Unresolved resonances in %s form\n" % self.nativeData.moniker )

    def check( self, info ):
        from fudge.gnd import warning
        warnings = []
        for L in self.tabulatedWidths.L_values:
            for J in L.J_values:
                elist = J.energyDependentWidths.getColumn('energy','eV')
                if elist is None:
                    warnings.append( warning.URRmissingEnergyList( L.L, J.J.value, J ) )
                else:
                    if elist[0] > self.lowerBound.getValueAs('eV') or elist[-1] < self.upperBound.getValueAs('eV'):
                        warnings.append( warning.URRdomainMismatch( L.L, J.J.value, J ) )
                    missingPoints = [i for i in range(1,len(elist)) if elist[i] > 3*elist[i-1]]
                    for idx in missingPoints:
                        warnings.append( warning.URRinsufficientEnergyGrid( L.L, J.J.value, PQU(elist[idx-1],'eV'),
                            PQU(elist[idx],'eV'), J ) )
        return warnings
    
    def toXMLList( self, flags, indent = '' ):
        xmlString = [indent+'<unresolved lowerBound="%s" upperBound="%s" nativeData="%s">' % (
            self.lowerBound, self.upperBound, self.nativeData.moniker ) ]
        xmlString += self.nativeData.toXMLList( flags, indent = indent + '  ' )
        xmlString[-1] += '</unresolved>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        xPath.append( element.tag )
        L_values = []
        for lval in element[0]:
            J_values = []
            for jval in lval:
                consWid, edepWid = {}, None
                for width in jval:
                    if width.tag=='constantWidths':
                        consWid = getAttrs( width )
                    else: edepWid = gndTable.parseXMLNode( width[0], conversionTable={'index':int} )
                if edepWid is None: edepWid = gndTable.table()
                for quant in ('energy','levelSpacing','neutronWidth','captureWidth','fissionWidthA','competitiveWidth'):
                    if ( quant not in consWid and quant not in [e.name for e in edepWid.columns] ):
                        consWid[quant] = PQU(0,'eV')
                J_values.append( URR_Jsection( xParticle.spin(jval.get('J')), edepWid, consWid,
                    **getAttrs( jval, exclude=("J"), required=("neutronDOF","gammaDOF","fissionDOF","competitiveDOF") ) ) )
            L_values.append( URR_Lsection( int(lval.get("L")), J_values ) )

        table = unresolvedTabulatedWidths( L_values, **getAttrs( element[0], required=('forSelfShieldingOnly',) ) )
        URR = unresolved( nativeData=table, **getAttrs( element, exclude=('nativeData',) ) )
        xPath.pop()
        return URR

class energyInterval:
    """ resolved region may be made up of multiple energy intervals (deprecated) """
    def __init__(self,index=None,nativeData=None,lowerBound='',upperBound=''):
        self.index = index
        self.nativeData = nativeData
        self.lowerBound = lowerBound
        self.upperBound = upperBound
    
    def setIndex(self, index):
        self.index = index
    
    def toString(self, simpleString = False):
        return ("%s resonances, %s to %s. Contains %i resonances" % 
                (self.lowerBound,self.upperBound, len(self.nativeData) ) )
    
    def toXMLList( self, flags, indent = ''):
        xmlString = [indent+
                '<region index="%s" lowerBound="%s" upperBound="%s" nativeData="%s">' 
                % (self.index,self.lowerBound,self.upperBound,self.nativeData.moniker ) ]
        if self.nativeData:
            xmlString += self.nativeData.toXMLList( flags, indent = indent+'  ' )
        else:
            raise Exception( "Resonance section contains no data!" )
        xmlString[-1] += '</region>'
        return xmlString


################################################################
# now we come to specific implementations for resolved region: #
################################################################
class resonanceFormalismBaseClass( ancestry ) :

    optAttrList = ( 'scatteringRadius', 'calculateChannelRadius', 'computeAngularDistribution', 'LvaluesNeededForConvergence' )
        # L-dependent AP also optional, but shouldn't appear in attributes list:
    LdepOpt = ('LdependentScatteringRadii',)
    
    def __init__(self, resonanceParameters=None, scatteringRadius=None, **kwargs):
        """LdependentScatteringRadii is a list of dictionaries, one for each L with APL != AP
        only used in ENDF for Reich_Moore form
        
        resonanceParameters should be a gnd.core.math.table.table instance."""
        
        index = 0
        for attr in self.optAttrList + self.LdepOpt:
            #if kwargs.get(attr):
            setattr( self, attr, kwargs.get(attr) )
        if self.computeAngularDistribution:
            self.computeAngularDistribution = bool(self.computeAngularDistribution)
        self.resonanceParameters = resonanceParameters or []
        if self.resonanceParameters:
            self.resonanceParameters.parent = self
            self.resonanceParameters.moniker = 'resonanceParameters'
        self.scatteringRadius = scatteringRadius
        if self.scatteringRadius: self.scatteringRadius.parent = self
        ancestry.__init__( self, self.moniker, None )
    
    def __getitem__(self, idx):
        return self.resonanceParameters[idx]
    
    def __len__(self):
        return len(self.resonanceParameters)
    
    def addResonance(self, resonance):
        """ insert a new resonance in the resonance parameter table """
        #resonance = (energy, J, l, ... )
        self.resonanceParameters.addRow( resonance )

    def getScatteringRadius(self, L=None):
        if ( (L is not None) and (getattr(self, 'LdependentScatteringRadii') is not None)
                and (L in self.LdependentScatteringRadii) ):
            AP = self.LdependentScatteringRadii[L]
        else: AP = self.scatteringRadius
        if AP.value == 0:
            print("WARNING: using scattering radius of 0.0 fm!")
        return AP
    
    def toXMLList( self, flags, indent = '' ):
        if not self.moniker :
            raise NotImplementedError ("Please use specific formalisms (MLBW, RM, etc) instead of the BaseClass")
        xmlString = '%s<%s' % ( indent, self.moniker )
        for attr in self.optAttrList:
            if getattr(self,attr,None):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal )
        xmlString = [xmlString+'>']
        xmlString += self.scatteringRadius.toXMLList( flags, indent=indent+'  ')
        if self.LdependentScatteringRadii:
            for Lval in self.LdependentScatteringRadii:
                xmlString.append('%s<LdependentScatteringRadius L="%s" radius="%s"/>' % (
                    indent+'  ', Lval, self.LdependentScatteringRadii[Lval]) )
        if self.resonanceParameters:
            xmlString.extend( ['%s<resonanceParameters>' % (indent+'  ')] +
                    self.resonanceParameters.toXMLList( indent+'    ' ) )
        xmlString[-1] += '</resonanceParameters></%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, linkData={} ):
        resonanceParameters = gndTable.parseXMLNode( element.find('resonanceParameters/table'),
                conversionTable={'index':int, 'L':int} )
        attrs = getAttrs( element )
        if attrs.get('scatteringRadius') is None:
            attrs['scatteringRadius'] = scatteringRadius( PQU(0,'fm') )
        elif attrs.get('scatteringRadius')=='energyDependent':
            attrs['scatteringRadius'] = scatteringRadius.parseXMLNode( element.find('scatteringRadius') )
        radii = [getAttrs(el) for el in element.findall('LdependentScatteringRadius')]
        if radii:
            attrs['LdependentScatteringRadii'] = dict( [(rad['L'],rad['radius']) for rad in radii] )
        return cls( resonanceParameters, **attrs )
    
    def toENDF6( self, flags, targetInfo, verbosityIndent='' ):
        endf = []
        AP = getattr( self, 'scatteringRadius' )
        if AP.isEnergyDependent():
            scatRadius = AP.value
            NR, NP = 1, len(scatRadius)
            endf.append( endfFormats.endfHeadLine( 0,0,0,0, NR,NP ) )
            endf += endfFormats.endfInterpolationList( (NP,
                endfFormats.twoAxesToENDFInterpolation(scatRadius.axes,0)) )
            endf += endfFormats.endfNdDataList( scatRadius.convertAxisToUnit( 1, '10*fm' ) )
            AP = 0
        else: AP = self.scatteringRadius.value.getValueAs('10*fm')
        L_list = self.resonanceParameters.getColumn('L')
        NLS = len( set(L_list) )
        LAD = getattr(self, 'computeAngularDistribution') or 0
        NLSC = getattr(self, 'LvaluesNeededForConvergence') or 0
        endf += [endfFormats.endfHeadLine( targetInfo['spin'], AP, LAD, 0, NLS, NLSC )]
        
        table = [ self.resonanceParameters.getColumn('energy',units='eV'),
                self.resonanceParameters.getColumn('J') ]
        NE = len(table[0])
        if "BreitWigner" in self.moniker :
            attrList = ('totalWidth','neutronWidth','captureWidth','fissionWidthA')
        elif "Reich_Moore" in self.moniker :
            attrList = ('neutronWidth','captureWidth','fissionWidthA','fissionWidthB')
        for attr in attrList:
            column = self.resonanceParameters.getColumn( attr, units='eV' )
            if not column: column = [0]*NE
            table.append( column )
        CS = self.resonanceParameters.getColumn('channelSpin')
        if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
            targetSpin = targetInfo['spin']
            CS = [2*(cs-targetSpin) for cs in CS]
            Js = [v[0]*v[1] for v in zip(table[1],CS)]
            table[1] = Js
        table = zip(*table)

        for L in set(L_list):
            APL = 0
            if self.LdependentScatteringRadii is not None and L in self.LdependentScatteringRadii:
                APL = self.LdependentScatteringRadii[L].getValueAs('10*fm')
            resonances = [table[i] for i in range(NE) if L_list[i] == L]
            NRS = len(resonances)
            endf.append( endfFormats.endfHeadLine( targetInfo['mass'], APL, L, 0, 6*NRS, NRS ) )
            for row in resonances:
                endf.append( endfFormats.endfDataLine( row ) )
        return endf


class SLBW( resonanceFormalismBaseClass ) :
    """Resonance region in Single-Level Breit-Wigner form."""

    moniker = "SingleLevel_BreitWigner"

    def __init__(self, resonanceParameters, **kwargs):

        resonanceFormalismBaseClass.__init__( self, resonanceParameters, **kwargs )

class MLBW( resonanceFormalismBaseClass ) :
    """Resonance region in Multi-Level Breit-Wigner form."""

    moniker = "MultiLevel_BreitWigner"

    def __init__(self, resonanceParameters, **kwargs):

        resonanceFormalismBaseClass.__init__( self, resonanceParameters, **kwargs )

class RM( resonanceFormalismBaseClass ) :
    """Reich_Moore formalism."""

    moniker = "Reich_Moore"

    def __init__(self, resonanceParameters, **kwargs):

        resonanceFormalismBaseClass.__init__( self, resonanceParameters, **kwargs )

#############################################
# R-matrix (LRF=7 in ENDF) is more complex: #
#############################################
class RMatrix:
    """R-Matrix Limited formalism. """

    moniker = "R_Matrix_Limited"

    optAttrList = ('approximation','scatteringRadius','boundaryCondition','relativisticKinematics',
            'reducedWidthAmplitudes','calculatePenetrability','calculateShift','calculateChannelRadius')

    def __init__(self, openChannels=None, spinGroups=None, **kwargs):
        """R-matrix is sorted by 'spin groups', each with same Jpi (conserved)."""

        self.openChannels = openChannels or []
        self.spinGroups = spinGroups or []
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))

    def __getitem__(self, idx):
        return self.spinGroups[idx]

    def __len__(self):
        return len(self.spinGroups)

    def toXMLList( self, flags, indent = '' ):
        xmlString = ['%s<%s' % ( indent, self.moniker ) ]
        for attr in self.optAttrList:
            if getattr(self,attr):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString[0] += ' %s="%s"' % (attr,attrVal)
        xmlString[0] += '>'
        indent += '  '
        for ch in self.openChannels:
            xmlString.append(ch.toXMLList(indent=indent))
        for gr in self.spinGroups:
            xmlString.append(gr.toXMLList(indent=indent))
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        openChannels = [openChannel.parseXMLNode( el ) for el in element if el.tag=='openChannel']
        spinGroups = []
        for sg in [ch for ch in element if ch.tag=='spinGroup']:
            resonanceParameters = gndTable.parseXMLNode( sg[0][0], conversionTable={'index':int, 'L':int,
                'channelSpin':float, 'effectiveRadius':PQU,} )
            spinGroups.append( spinGroup( resonanceParameters=resonanceParameters, **getAttrs(sg,
                required=('background',) ) ) )
        return RMatrix( openChannels, spinGroups, **getAttrs(element, 
            required=RMatrix.optAttrList ) )

    def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):
        """ """
        KRM = {'SLBW':1, 'MLBW':2, 'Reich_Moore':3, 'Full R-Matrix':4} [self.approximation]
        endf = [endfFormats.endfHeadLine(0,0,self.reducedWidthAmplitudes,KRM,
            len(self.spinGroups),self.relativisticKinematics)]
        
        # first describe all the particle pairs (two-body output channels)
        NPP = len(self.openChannels)
        endf.append( endfFormats.endfHeadLine(0,0,NPP,0,12*NPP,2*NPP) )

        def getENDFtuple( spin, parity ):
            # ENDF combines spin & parity UNLESS spin==0. Then it wants (0,+/-1)
            if spin.value: return (spin.value * parity.value, 0)
            else: return (spin.value, parity.value)
        
        def MZIP( p ):  # helper: extract mass, z, spin and parity from particle list
            mass = p.getMass( 'amu' ) / targetInfo['reactionSuite'].getParticle( 'n' ).getMass( 'amu' )
            Z = p.getZ_A_SuffixAndZA()[0]
            I,P = getENDFtuple( p.getSpin(), p.getParity() )
            return mass, Z, I, P
        
        for pp in self.openChannels:
            pA,pB = pp.channel.split(' + ')
            # get the xParticle instances for pA and pB:
            pA,pB = targetInfo['reactionSuite'].getParticle(pA), targetInfo['reactionSuite'].getParticle(pB)
            MA, ZA, IA, PA = MZIP( pA )
            MB, ZB, IB, PB = MZIP( pB )
            MT = pp.ENDF_MT
            PNT = pp.calculatePenetrability
            if PNT is None: PNT = self.calculatePenetrability
            SHF = pp.calculateShift
            if SHF is None: SHF = self.calculateShift
            if MT in (19,102): PNT = 0  # special case
            try: Q = pp.Qvalue.inUnitsOf('eV').value
            except:
                Q = targetInfo['reactionSuite'].getReaction( pp.channel ).getQ('eV')
                # getQ doesn't account for residual left in excited state:
                for particle in targetInfo['reactionSuite'].getReaction( pp.channel ).outputChannel.particles:
                    Q -= particle.getLevelAsFloat('eV')
            if MT==102: Q = 0
            endf.append( endfFormats.endfDataLine( [MA,MB,ZA,ZB,IA,IB] ) )
            endf.append( endfFormats.endfDataLine( [Q,PNT,SHF,MT,PA,PB] ) )
        
        for spingrp in self.spinGroups:
            AJ,PJ = getENDFtuple( spingrp.spin, spingrp.parity )
            KBK = spingrp.background
            KPS = spingrp.applyPhaseShift
            NCH = len(spingrp.resonanceParameters.columns)-1    # skip the 'energy' column
            endf.append( endfFormats.endfHeadLine( AJ,PJ,KBK,KPS,6*NCH,NCH ) )
            for chan in spingrp.resonanceParameters.columns:
                # which open channel does it correspond to?
                if chan.name == 'energy': continue
                name = chan.name.split(' width')[0]
                openChannel = [ch for ch in self.openChannels if ch.channel == name][0]
                PPI = self.openChannels.index( openChannel )+1  # 1-based index in ENDF
                attr = chan.attributes
                L = attr['L']
                SCH = attr['channelSpin']
                # some data may have been moved up to the channel list:
                BND = attr.get('boundaryCondition') or self.boundaryCondition
                APT = attr.get('scatteringRadius') or self.scatteringRadius
                APE = attr.get('effectiveRadius') or openChannel.effectiveRadius or APT
                APT = APT.getValueAs('10*fm')
                APE = APE.getValueAs('10*fm')
                if openChannel.ENDF_MT==102:
                    APT, APE = 0,0
                endf.append( endfFormats.endfDataLine( [PPI,L,SCH,BND,APE,APT] ) )
            # resonances:
            NRS = len(spingrp.resonanceParameters)
            NX = (NCH//6 + 1)*NRS
            if NRS==0: NX=1 # special case
            endf.append( endfFormats.endfHeadLine( 0,0,0,NRS,6*NX,NX ) )
            for res in spingrp.resonanceParameters:
                #resproperties = [res.energy.inUnitsOf('eV').value] + [
                #        w.inUnitsOf('eV').value for w in res.widths]
                for jidx in range(NCH//6+1):
                    endfLine = res[jidx*6:jidx*6+6]
                    while len(endfLine)<6: endfLine.append(0)
                    endf.append( endfFormats.endfDataLine( endfLine ) )
            if NRS==0:
                endf.append( endfFormats.endfDataLine( [0,0,0,0,0,0] ) )
        return endf


class openChannel:
    """ describes one reaction channel that opens up in the resonance region. In an R-Matrix section,
    all open reaction channels should be described in the list of openChannel elements """
    reqAttrList = ('channel', 'ENDF_MT')
    optAttrList = ('Qvalue', 'scatteringRadius', 'effectiveRadius','boundaryCondition',
            'calculatePenetrability','calculateShift')
    def __init__(self, index, **kwargs):
        self.index = index
        for attr in self.reqAttrList:
            setattr(self, attr, kwargs[attr])
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr,None))
    
    def toXMLList( self, indent = ''):
        xmlString = ('%s<openChannel index="%s"') % (indent, self.index)
        for attr in self.reqAttrList:
            xmlString += ' %s="%s"' % (attr,getattr(self,attr))
        for attr in self.optAttrList:
            if getattr(self,attr):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal)
        xmlString += '/>'
        return xmlString

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        attrs = getAttrs( element )
        if 'Qvalue' not in attrs:
            attrs['Qvalue'] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(0,'eV')
        return openChannel( attrs.pop('index'), **attrs )


class spinGroup:
    """ single group with same Jpi (conserved). spin group contains AP
    (scattering radius), 1 or more resonance widths """
    def __init__(self, index=None, spin=None, parity=None, background=None, applyPhaseShift=False,
            resonanceParameters=None):
        self.index = index
        self.spin = spin
        self.parity = parity
        self.background = background
        self.applyPhaseShift = applyPhaseShift
        self.resonanceParameters = resonanceParameters
    
    def __getitem__(self, idx):
        return self.resonanceParameters[idx]
    
    def __len__(self):
        return len(self.resonanceParameters)
    
    def __lt__(self, other):
        # for sorting spin groups by Jpi. group J values together
        return (self.Jpi.spin,self.Jpi.parity) < (
                other.Jpi.spin,other.Jpi.parity)
    
    def toXMLList( self, indent='' ):
        xmlString = '%s<spinGroup index="%i" spin="%s" parity="%s"' % (indent, self.index, self.spin, self.parity)
        for opt in 'background','applyPhaseShift':
            if getattr(self, opt):
                attrVal = getattr(self,opt)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (opt, attrVal)
        xmlString = [xmlString + '>']
        #for ch in self.channelInfo:
        #    xmlString.append( ch.toXMLList(indent=indent+'  ') )
        if self.resonanceParameters.columns:    # need not contain any data
            xmlString.extend( ['%s<resonanceParameters>' % (indent+'  ')] +
                    self.resonanceParameters.toXMLList( indent+'    ' ) )
        xmlString[-1] += '</resonanceParameters></spinGroup>'
        return '\n'.join(xmlString)

#############################################
#    end of R-Matrix (LRF=7) classes        #
#############################################


####################################################
# remaining classes are for unresolved resonances: #
####################################################
class unresolvedTabulatedWidths:

    moniker = 'tabulatedWidths'
    optAttrList = ('interpolation','scatteringRadius','forSelfShieldingOnly','ENDFconversionFlag')

    def __init__(self, L_values=None, **kwargs):
        """ L_values is a list of URR_Lsections, which in turn contains a list of J_sections. 
        The average widths are given in the J_sections """
        index = 0
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))
        self.L_values = L_values or []
    
    def __len__(self):
        return len(self.L_values)
    
    def toXMLList( self, flags, indent = '' ):
        xmlString = '%s<%s' % (indent, self.moniker)
        for attr in self.optAttrList:
            if getattr(self,attr,None):
                attrVal = getattr(self,attr)
                if type(attrVal) is bool: attrVal = str(attrVal).lower()
                xmlString += ' %s="%s"' % (attr, attrVal)
        xmlString = [xmlString+'>']
        for L in self.L_values:
            xmlString += L.toXMLList( flags, indent = indent+'  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString
    
    def toENDF6( self, flags, targetInfo, verbosityIndent = ''):
        AP = self.scatteringRadius.value.inUnitsOf('10*fm').value
        NLS = len(self.L_values)
        LFW = targetInfo['unresolved_LFW']; LRF = targetInfo['unresolved_LRF']
        
        def v(val): # get value back from PhysicalQuantityWithUncertainty
            if type(val)==type(None): return
            return val.getValueAs('eV')
        
        if LFW==0 and LRF==1:   # 'Case A' from ENDF 2010 manual pg 70
            endf = [endfFormats.endfHeadLine( targetInfo['spin'], AP, 
                self.forSelfShieldingOnly,0,NLS,0) ]
            for Lval in self.L_values:
                NJS = len(Lval.J_values)
                endf.append(endfFormats.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, 6*NJS, NJS ))
                for Jval in Lval.J_values:
                    # here we only have one width per J:
                    ave = Jval.constantWidths
                    endf.append( endfFormats.endfDataLine([v(ave.levelSpacing),Jval.J.value,
                        Jval.neutronDOF,v(ave.neutronWidth),v(ave.captureWidth),0]) )
        
        elif LFW==1 and LRF==1: # 'Case B'
            energies = self.L_values[0].J_values[0].energyDependentWidths.getColumn('energy',units='eV')
            NE = len(energies)
            endf = [endfFormats.endfHeadLine( targetInfo['spin'], AP,
                self.forSelfShieldingOnly, 0, NE, NLS )]
            nlines = int(math.ceil(NE/6.0))
            for line in range(nlines):
                endf.append( endfFormats.endfDataLine( energies[line*6:line*6+6] ) )
            for Lval in self.L_values:
                NJS = len(Lval.J_values)
                endf.append( endfFormats.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ) )
                for Jval in Lval.J_values:
                    cw = Jval.constantWidths
                    endf.append( endfFormats.endfHeadLine(0,0,Lval.L,Jval.fissionDOF,NE+6,0) )
                    endf.append( endfFormats.endfDataLine([v(cw.levelSpacing),Jval.J.value,Jval.neutronDOF,
                        v(cw.neutronWidth),v(cw.captureWidth),0]) )
                    fissWidths = Jval.energyDependentWidths.getColumn('fissionWidthA',units='eV')
                    for line in range(nlines):
                        endf.append( endfFormats.endfDataLine( fissWidths[line*6:line*6+6] ) )
        
        elif LRF==2:            # 'Case C', most common in ENDF
            endf = [endfFormats.endfHeadLine( targetInfo['spin'], AP, 
                self.forSelfShieldingOnly,0,NLS,0) ]
            INT = endfFormats.XYStringToENDFInterpolation( self.interpolation )
            for Lval in self.L_values:
                NJS = len(Lval.J_values)
                endf.append( endfFormats.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ))
                for Jval in Lval.J_values:
                    NE = len( Jval.energyDependentWidths )
                    useConstant = not NE
                    if useConstant: NE = 2
                    endf.append( endfFormats.endfHeadLine( Jval.J.value, 0, INT, 0, 6*NE+6, NE ) )
                    endf.append( endfFormats.endfDataLine([0,0,Jval.competitiveDOF,
                        Jval.neutronDOF,Jval.gammaDOF,Jval.fissionDOF]) )
                    cws = Jval.constantWidths
                    if useConstant:
                        # widths are all stored in 'constantWidths' instead. get energies from parent class
                        NE = 2; useConstant=True
                        eLow,eHigh = targetInfo['regionEnergyBounds']
                        for e in (eLow,eHigh):
                            endf.append(endfFormats.endfDataLine([v(e),v(cws.levelSpacing),
                                v(cws.competitiveWidth),v(cws.neutronWidth),v(cws.captureWidth),
                                v(cws.fissionWidthA)]) )

                    else:
                        table = [ Jval.energyDependentWidths.getColumn('energy',units='eV') ]
                        for attr in ('levelSpacing','competitiveWidth','neutronWidth','captureWidth',
                                'fissionWidthA'):
                            # find each attribute, in energy-dependent or constant width section
                            column = ( Jval.energyDependentWidths.getColumn( attr, units='eV' ) or
                                    [v(getattr( cws, attr ))]*NE )
                            if not any(column): column = [0]*NE
                            table.append( column )
                        for row in zip(*table):
                            endf.append( endfFormats.endfDataLine( row ) )
        return endf


# each of the above sections contains one or more 'L' sections, which each contain one or more 'J's:
class URR_Lsection:
    """ unresolved average widths, grouped by L. Contains list of J-values: """
    def __init__(self, L, J_values):
        self.L = L
        self.J_values = J_values
    
    def toXMLList( self, flags, indent = '' ):
        xmlString = ['%s<L_section L="%s">' % (indent, self.L)]
        for J in self.J_values:
            xmlString += J.toXMLList( flags, indent+'  ' )
        xmlString[-1] += '</L_section>'
        return xmlString


class URR_Jsection:
    """ unresolved average widths, grouped by J. Should contain *either*
     - energyDependentWidths, with a table of widths / level densities as function of energy, or
     - constantWidths + energyGrid (the grid is necessary to tell where to do the reconstruction) """
    optAttrList = ('competitiveDOF','neutronDOF','gammaDOF','fissionDOF')
    def __init__(self, J, energyDependentWidths=None, constantWidths={}, **kwargs):
        self.J = J
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr))
        self.constantWidths = constantWidthSection( **constantWidths )
        self.energyDependentWidths = energyDependentWidths
    
    def eliminateRedundantInfo(self):
        """ ENDF often lists more info than necessary. 
        If we encounter section where everything (except energy) is constant, 
        move data from 'energyDependentWidths' to 'constantWidths' """

        allEliminated = False
        edep = self.energyDependentWidths
        for colId in range(edep.nColumns)[::-1]:
            column = edep.columns[colId]
            columnData = edep.getColumn( column.name, units='eV' )
            if len(set( columnData ) ) == 1:
                setattr( self.constantWidths, column.name, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(
                    columnData[0], column.units ) )
                [d.pop(colId) for d in edep.data]
                edep.columns.pop(colId)
        for idx, col in enumerate( edep.columns ): col.index = idx  #re-number
        #if edep.nColumns == 1 and edep.columns[0].name == 'energy':
        #    edep.columns, edep.data = [],[] # all widths are constant
        #    allEliminated = True
        return allEliminated
    
    def toXMLList( self, flags, indent='' ):
        indent2 = indent+'  '; indent3 = indent+'    '
        xmlString = ['%s<J_section J="%s"' % (indent,self.J) ]
        for attr in URR_Jsection.optAttrList:
            if getattr(self, attr):
                xmlString[-1] += ' %s="%s"' % (attr, getattr(self,attr))
        xmlString[-1] += '>'
        constwidth = self.constantWidths.toXMLList( indent=indent2 )
        if constwidth: xmlString.append(constwidth)
        if self.energyDependentWidths:
            xmlString.extend( ['%s<energyDependentWidths>' % (indent2)] +
                    self.energyDependentWidths.toXMLList( indent3 ) )
            xmlString[-1] += '</energyDependentWidths>'
        xmlString[-1] += '</J_section>'
        return xmlString


# single region (specified by L, J and optionally energy) for which we provide average widths in URR:
class constantWidthSection:
    """
    for widths or level spacings that depend only on L and J (not on energy)
    """
    optAttrList = ('levelSpacing','neutronWidth','captureWidth','competitiveWidth',
            'fissionWidthA','fissionWidthB')
    def __init__(self, **kwargs):
        for attr in self.optAttrList:
            setattr(self, attr, kwargs.get(attr, None))
    
    def toXMLList( self, indent = '' ):
        if not any( [ getattr(self, attr) for attr in self.optAttrList ] ):
            return ''
        xmlString = '%s<constantWidths' % indent
        for attr in self.optAttrList:
            if getattr(self, attr):
                xmlString += ' %s="%s"' % (attr, getattr(self, attr))
        xmlString += '/>'
        return xmlString


# helper functions for reading in from xml:
def getBool( value ):
    return {'true':True, '1':True, 'false':False, '0':False}[value]
def floatOrint( value ):
    if float( value ).is_integer(): return int( value )
    return float( value )
def getAttrs(element, exclude=(), required=()):
    """ read attributes from one xml element. Convert to PhysicalQuantityWithUncertainty (or bool, int, etc)
    where appropriate, and skip anything in the 'exclude' list: """
    conversionTable = {'lowerBound':PQU, 'upperBound':PQU, 'value':PQU, 'energy':PQU,
            'neutronWidth':PQU, 'captureWidth':PQU, 'fissionWidthA':PQU, 'fissionWidthB':PQU, 'competitiveWidth':PQU,
            'levelSpacing':PQU, 'Qvalue':PQU, 'radius':PQU, 'effectiveRadius':PQU,
            'reconstructCrossSection':getBool, 'multipleRegions': getBool, 'LdependentScatteringRadii': getBool,
            'calculateChannelRadius':getBool, 'computeAngularDistribution':getBool, 'forSelfShieldingOnly': getBool,
            'calculateShift':getBool,'calculatePenetrability':getBool,
            'LvaluesNeededForConvergence':int, 'ENDF_MT':int, 'index':int, 'L':int,
            'neutronDOF':floatOrint, 'gammaDOF':floatOrint, 'competitiveDOF':floatOrint, 'fissionDOF':floatOrint,
            'spin':xParticle.spin, 'parity':xParticle.parity,
            'scatteringRadius':(lambda foo: scatteringRadius(PQU(foo)) if foo!='energyDependent' else foo),
            }
    attrs = dict( element.items() )
    for key in attrs.keys():
        if key in exclude: attrs.pop(key)
        elif key in conversionTable: attrs[key] = conversionTable[key]( attrs[key] )
    for val in required:
        if val not in attrs: attrs[val] = False
    return attrs

