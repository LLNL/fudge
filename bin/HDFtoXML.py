#!/usr/bin/env python

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
translate hdf5 to xml
cmattoon, 1/11/2014
"""
import sys
import os
import numpy

from xml.etree import cElementTree as etree
import h5py

from pqu.physicalQuantityWithUncertainty import toShortestString

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
        xmlnode.text = ' '.join( map( toShortestString, data ) )
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
    xmldata.text = '\n'.join( [ ' '.join( map( toShortestString, line ) ) for line in node['data'].value ] )
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
