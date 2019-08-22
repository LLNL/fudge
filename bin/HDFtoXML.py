#!/usr/bin/env python

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
translate hdf5 to xml
cmattoon, 1/11/2014
"""
import sys
import os
import numpy

from xml.etree import cElementTree as etree
import h5py

from pqu import PQU

class ElementTreeCDATA( etree.ElementTree ):
    # must wrap text data in <![CDATA[ ... ]]>, or special characters like <,>,/,& will cause problems
    def _write(self, file, node, encoding, namespaces):
        text = node.text.encode(encoding)
        file.write("\n<![CDATA[%s]]>\n" % text)

def addNode( parent, node ):
    """
    recursively copy HDF data (one node at a time) to xml
    parent is an xml element,
    node is an HDF5 group to be added to xml
    """
    name = os.path.split( node.name )[-1]
    attrs = dict( node.attrs )  # make a copy
    # some names must be modified since hdf5 doesn't allow duplicates
    if 'tag' in attrs:
        name = attrs.pop('tag')
    name = name.replace('_','/')    # "1/2" not allowed in hdf5 group name

    try:
        xmlnode = etree.Element( name )
        for key,val in attrs.items():
            xmlnode.set( str(key), str(val) )
        if isinstance( parent, etree.ElementTree ): parent._setroot( xmlnode )
        else: parent.append( xmlnode )
        newParent = xmlnode

        if node.attrs.get('xData') is not None:
            addXData( newParent, node )
        elif name=='table':
            addTable( newParent, node )
        elif name=='documentation':
            addDocumentation( newParent, node )
        else:
            for child in node.values():
                addNode( newParent, child )
    except:
        print ("addNode backtrace: node=%s, name=%s" % (node,name))
        raise

def addXData( parent, node ):
    """
    called when elements with 'xData' attributes are encountered
    """
    def addDataSet( parent, dataNode ):
        name = os.path.split( dataNode.name )[-1]
        attrs = dict( dataNode.attrs )
        if 'tag' in attrs:
            name = attrs.pop('tag')
        xmlnode = etree.Element( name )
        for key,val in attrs.items():
            xmlnode.set( str(key), str(val) )

        data = dataNode.value[:]
        data.shape = numpy.product( data.shape )
        xmlnode.text = ' '.join( map( PQU.toShortestString, data ) )
        parent.append( xmlnode )

    # xData always needs axes info:
    addNode( parent, node['axes'] )

    xData = node.attrs.get('xData')
    if xData=='XYs': addDataSet( parent, node['data'] )
    elif xData=='polynomial': addDataSet( parent, node['data'] )
    elif xData in ('W_XYs', 'W_XYs_LegendreSeries'):
        idx = 0
        while True:
            energy_in = node.get('energy_in: %i' % idx)
            if energy_in is None: break
            addDataSet( parent, energy_in )
            idx += 1
    elif xData in ('V_W_XYs', 'V_W_XYs_LegendreSeries'):
        idx = 0
        while True:
            energy_in = node.get('energy_in: %i' % idx)
            if energy_in is None: break
            xmlenergy_in = etree.Element( 'energy_in' )
            for key,val in energy_in.attrs.items():
                xmlenergy_in.set( str(key), str(val) )
            jdx = 0
            while True:
                W = energy_in.get('energy_out: %i' % jdx) or energy_in.get('mu: %i' % jdx)
                if W is None: break
                addDataSet( xmlenergy_in, W )
                jdx += 1
            parent.append( xmlenergy_in )
            idx += 1
    elif xData in ('regionsXYs', 'regionsW_XYs_LegendreSeries'):
        idx = 0
        while True:
            region = node.get('region: %i' % idx)
            if region is None: break
            xmlregion = etree.Element( 'region' )
            for key,val in region.attrs.items():
                xmlregion.set( str(key), str(val) )
            addNode( xmlregion, region['interpolationAxes'] )
            if xData=='regionsXYs':
                addDataSet( xmlregion, region['data'] )
            else: # piecewise Legendre data
                jdx = 0
                while True:
                    energy_in = region.get('energy_in: %i' % jdx)
                    if energy_in is None: break
                    addDataSet( xmlregion, energy_in )
                    jdx += 1
            parent.append( xmlregion )
            idx += 1
    else:
        raise Exception("Encountered unknown form of xData: %s" % xData)

def addTable( parent, node ):
    addNode( parent, node['columnHeaders'] )
    dimensions = map(int, [node.attrs.get(val) for val in ('rows','columns')])
    xmldata = etree.Element('data')
    xmldata.text = '\n'.join( [ ' '.join( map( PQU.toShortestString, line ) ) for line in node['data'].value ] )
    parent.append( xmldata )

def addDocumentation( parent, docNode ):
    parent.text = docNode['text'].value
    parent.__class__ = ElementTreeCDATA


if __name__ == '__main__':
    if len(sys.argv)<2:
        print ("Usage: python HDFtoXML.py <HDF_filename>")
        sys.exit(1)

    h5file = sys.argv[1]
    h5 = h5py.File( h5file )

    xdoc = etree.ElementTree()

    root = h5.values()[0]
    addNode( xdoc, root )

    xdoc.write( os.path.splitext(h5file)[0] + ".xml" )
