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
translate xml to hdf5
cmattoon, 10/21/2010

updates:
October 2011, check the 'xData' attribute to determine how to store data
June 2011: use etree instead of DOM
"""
import sys
import os
import numpy

from xml.etree import cElementTree as parser
import h5py

keydir = {
        'alias': 'key',
        'particle': 'name',
        'level': 'label',
        'gamma': 'finalLevel',
        'openChannel': 'index',
        'spinGroup': 'index',
        'J_section': 'J',
        'L_section': 'L',
        'LdependentScatteringRadius': 'L',
        'column': 'index',
        'reaction': 'label',
        'summedReaction': 'label',
        'fissionComponent': 'label',
        'production': 'label',
        'product': 'label',
        'region': 'index',
        'weighted': 'index',
        'key': 'label',
        'axis': 'index',
        'energy_in': 'index',
        'energy_out': 'index',
        'mu': 'index',
        'summand': '{http://www.w3.org/1999/xlink}href',
        # covariances:
        'covariance': 'id'
        }

def addNode( parent, node ):
    """
    recursively copy xml data (one <node> at a time) to hdf5
    parent is a hdf5 element,
    node is an xml.dom element to be added to hdf5
    """
    name = node.tag
    # some names must be modified since hdf5 doesn't allow duplicates
    if name in keydir.keys():
        node.set('tag', name)
        name = "%s: %s" % (name, node.get( keydir[node.tag] ))
    name = name.replace('/','_')    # "1/2" not allowed in hdf5 group name

    try:
        h5node = parent.create_group( name )
        for key,val in node.items():
            h5node.attrs.create(str(key),str(val))
        newParent = h5node

        if node.get('xData') is not None:
            addXData( newParent, node )
        elif name=='table':
            addTable( newParent, node )
        elif name=='documentation':
            addDocumentation( newParent, node )
        else:
            for child in node.getchildren():
                addNode( newParent, child )
    except:
        print ("addNode backtrace: node=%s, name=%s" % (node,name))
        raise

def addXData( parent, node ):
    """
    called when elements with 'xData' attributes are encountered
    Data should be first stored in numpy array, then dumped into hdf5 dataset
    """
    def addDataSet( parent, dataNode, twoDimensional=True ):
        data = numpy.array( map(float, dataNode.text.split()) )
        if twoDimensional: data.shape = (len(data)//2, 2)
        name = dataNode.tag
        if name in keydir:
            dataNode.set('tag', name)
            name="%s: %s" % (name, dataNode.get( keydir[dataNode.tag] ))
        h5data = parent.create_dataset( name, data=data )
        for key,val in dataNode.items(): h5data.attrs.create(str(key),str(val))

    # xData always needs axes info:
    addNode( parent, node.find('axes') )

    xData = node.get('xData')
    if xData=='XYs': addDataSet( parent, node.find('data') )
    elif xData=='polynomial': addDataSet( parent, node.find('data') )
    elif xData in ('W_XYs', 'W_XYs_LegendreSeries'):
        twoDimensional = (xData=='W_XYs')   # otherwise it's a 1-d list of Legendre coefs
        for energy_in in node[1:]:
            addDataSet( parent, energy_in, twoDimensional )
    elif xData in ('V_W_XYs', 'V_W_XYs_LegendreSeries'):
        twoDimensional = (xData=='V_W_XYs')
        for energy_in in node[1:]:
            h5energy_in = parent.create_group( "energy_in: %s" % energy_in.get("index") )
            for key,val in energy_in.items(): h5energy_in.attrs.create(str(key),str(val))
            for W in energy_in: # could be 'mu' or 'energy_out'
                addDataSet( h5energy_in, W, twoDimensional )
    elif xData in ('regionsXYs', 'regionsW_XYs_LegendreSeries'):
        for region in node[1:]:
            h5region = parent.create_group( "region: %s" % region.get('index') )
            for key,val in region.items():
                h5region.attrs.create(str(key),str(val))
            addNode( h5region, region.find('interpolationAxes') )
            if xData=='regionsXYs':
                addDataSet( h5region, region.find('data'), twoDimensional=True )
            else: # piecewise Legendre data
                for energy_in in region[1:]:
                    addDataSet( h5region, energy_in, twoDimensional=False )
    else:
        raise Exception("Encountered unknown form of xData: %s" % xData)

def addTable( parent, node ):
    addNode( parent, node.find('columnHeaders') )
    dimensions = map(int, [node.get(val) for val in ('rows','columns')])
    data = numpy.array( map(float, node[1].text.split()) )
    data.shape = dimensions
    parent.create_dataset( 'data', data=data )

def addDocumentation( parent, docNode ):
    text = docNode.text
    h5node = parent.create_dataset( 'text', data=numpy.array(text) )


if __name__ == '__main__':
    if len(sys.argv)<2:
        print ("Usage: python toHDF5.py <xml_filename>")
        sys.exit(1)

    xmlfile = sys.argv[1]
    xdoc = parser.parse( xmlfile )

    h5file = os.path.splitext(xmlfile)[0] + ".hdf5"
    if h5py.version.version_tuple > (1,4,0,''): # use new, more efficient HDF5 format:
        h5 = h5py.File( h5file, "w", libver='latest' )
    else:
        h5 = h5py.File( h5file, "w" )

    root = xdoc.getroot()
    addNode( h5, root )

