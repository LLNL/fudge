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
from fudge.gnd import link
from fudge.core.ancestry import ancestry
from fudge.core.math import matrix as gndMatrix
from pqu import PQU

class inputParameter( ancestry ):
    """
    For use within a modelParameterCovariance: the rows/columns of that matrix each correspond to one input
    model parameter. Parameters must each be inputParameter or multipleInputParameter instances.
    """
    def __init__(self, name, path, unit=None):
        ancestry.__init__( self, 'inputParameter', None )
        self.name = name #: just a Python str
        self.path = path #: a fudge.gnd.link 
        self.unit = unit #: an acceptable PQU unit

    def toXMLList(self, flags=None, indent=''):
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

    def toXMLList(self, flags=None, indent=''):
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
        xmllist += self.matrix.toXMLList(flags=flags, indent=indent+'  ')
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

    def toXMLList(self, flags=None, indent=''):
        return [ indent +
                '<loopOverResonanceParameters nResonances="%i" parametersPerResonance="%s" xlink:href="%s"/>'
            % (self.nResonances, self.parametersPerResonance, self.path.toXLink()) ]

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        ll = link.parseXMLNode( element, linkData )
        path = link.follow( ll.path, linkData['reactionSuite'] )
        LORPs = loopOverResonanceParameters( int(ll.attributes["nResonances"]),
                ll.attributes["parametersPerResonance"], path )
        xPath.pop()
        return LORPs

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
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate <resonanceParameterCovariance> element from xml."""

        xPath.append( element.tag )
        params = []
        for param in element[0]:
            if param.tag=='loopOverResonanceParameters':
                param=loopOverResonanceParameters.parseXMLNode(param, xPath, linkData)
            elif param.tag=='parameter': pass
            params.append(param)
        Matrix = gndMatrix.parseXMLNode( element[1], xPath )
        RPCs = resonanceParameterCovariance( inputParameters=params, matrix=Matrix, **dict(element.items()) )
        xPath.pop()
        return RPCs

    def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
        """ go back to ENDF format """
        from fudge.legacy.converting import endfFormats
        def swaprows( matrix, i1, i2, nrows ):
            # may need to rearrange parameters: ENDF often sorts first by L rather than by energy
            rows = matrix[i1:i1+nrows].copy()
            matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
            cols = matrix[:,i1:i1+nrows].copy()
            matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

        # need the resonance parameters as well as covariance matrix:
        res = targetInfo['reactionSuite'].resonances
        RPs = res.resolved.nativeData.resonanceParameters
        NRes = self.inputParameters[-1].nResonances

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
            DAP = PQU( self.matrix.data[0][0], self.inputParameters[0].unit ).getValueAs('10*fm')
            if LRF in (1,2):
                endf.append( endfFormats.endfDataLine( [0,DAP] ) )
            elif LRF==3:
                endf.append( endfFormats.endfHeadLine( 0,0,0,0,MLS,1 ) )
                endf.append( endfFormats.endfDataLine( [DAP] ) )
            else:
                raise Exception("ISR>0 not yet supported for LRF=%i!" % LRF)

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
            tableIndex = 0
            for L in set( Ls ):
                NRS = Ls.count(L)
                endf.append( endfFormats.endfHeadLine( AWT, 0, L, 0, 18*NRS, NRS ) )
                for i in range(tableIndex, len(table)):
                    if table[i][0]!=L: break
                    endf.append( endfFormats.endfDataLine( table[i][1:7] ) )
                    block = matrix[MPAR*i:MPAR*(i+1), MPAR*i:MPAR*(i+1)]
                    lis = [block[0,0], block[1,1], block[2,1], block[2,2]]
                    if MPAR==4:
                        lis += [block[3,1],block[3,2],block[3,3],0,0,0,0,0]
                    else:
                        lis += [0,0,0,0,0,0,0,0]
                    endf += endfFormats.endfDataList( lis )
                tableIndex += NRS


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
                    params = table[i][1:7]
                    uncerts = [dat[MPAR*i],0,0,dat[MPAR*i+1],dat[MPAR*i+2],0]
                    if MPAR==4: uncerts[-1] = dat[MPAR*i+3]
                elif LRF==3:
                    params = table[i][1:3] + table[i][4:8]
                    uncerts = [dat[MPAR*i],0,dat[MPAR*i+1],dat[MPAR*i+2],0,0]
                    if MPAR==5: uncerts[-2:] = [dat[MPAR*i+3], dat[MPAR*i+4]]
                endf += endfFormats.endfDataList( params )
                endf += endfFormats.endfDataList( uncerts )

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
            NNN = NRes * MPAR
            NM = len(endfCorrMat)
            endf.append( endfFormats.endfHeadLine( 0,0, NDIGIT, NNN, NM, 0 ) )
            endf += endfCorrMat
        endf.append( endfFormats.endfSENDLineNumber() )
        endfMFList[32][151] = endf

