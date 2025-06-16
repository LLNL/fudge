#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Translate GNDS/XML to HDF5. Each xml element becomes a group in HDF, except a few elements that become
HDF datasets. Attributes are translated to HDF5 metadata
"""
import os
import argparse
import numpy
from collections import Counter

from xml.etree import ElementTree as ET
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("xml", help="XML file to be converted")
parser.add_argument("output", nargs='?',
        help="output file name. By default use the input file with extension changed to 'h5'")
parser.add_argument("-l", "--latestHDF", action='store_true',
        help="Use newer HDF format specification: smaller files but requires HDF 1.10 or later")
parser.add_argument("-c", "--compress", action='store_true',
        help="Use gzip compression (with shuffle) in HDF5 datasets")

h5Opts = {}

def fixName( node, suffix=None ):
    name = node.tag
    if suffix is not None:
        name += '_%s' % str(suffix)
    return name

def addAttributes( h5node, node, index ):
    for key,val in node.items():
        h5node.attrs.create(str(key),numpy.string_(val))
    h5node.attrs.create("_xmltag", numpy.string_(node.tag)) # save original tag for translating back to XML
    h5node.attrs.create("_xmlindex", numpy.uint32(index))   # hdf5 doesn't preserve group order, so preserve it here

def addNode( parent, node, index, suffix=None ):
    """
    recursively copy xml data (one <node> at a time) to hdf5

    :param parent: hdf5 file or group
    :param node:   xml node to be added to hdf5
    :param suffix: suffix to append to xml tag name (since hdf5 doesn't allow duplicate group names)
    :return:
    """
    name = fixName(node, suffix)

    try:
        h5node = parent.create_group( name )
        addAttributes(h5node, node, index)
        newParent = h5node

        child_names = Counter( [ child.tag for child in node ] )
        for key,value in child_names.items():
            if value == 1: child_names[key] = None
            else: child_names[key] = 0

        for idx,child in enumerate(node):
            suffix = child_names[ child.tag ]
            if suffix is not None: child_names[ child.tag ] += 1
            if child.tag == 'values':
                addValues( newParent, child, idx, suffix, isTwoDimensional=(node.tag=='XYs1d') )
            elif child.tag == 'data':
                addTableData( newParent, child, idx, int(node.get('rows')), int(node.get('columns')), suffix )
            elif child.tag=='documentation':
                addDocumentation( newParent, child, idx, suffix )
            else:
                addNode( newParent, child, idx, suffix )
    except:
        print( "addNode backtrace: node=%s, name=%s" % ( node, name ) )
        raise

def addValues( parent, node, index, suffix=None, isTwoDimensional=False ):
    """
    <values> elements in GNDS store a list of numbers. Convert to HDF5 dataset
    """
    name = fixName(node, suffix)

    if node.text is None:
        data = numpy.array([])
    else:
        data = numpy.array( list( map( float, node.text.split( ) ) ) )
    if isTwoDimensional:
        data.shape = (len(data)//2,2)
    h5data = parent.create_dataset( name, data=data, **h5Opts )
    addAttributes(h5data, node, index)

def addTableData( parent, node, index, rows, columns, suffix=None ):
    """
    <data> elements in GNDS store a table. Convert to HDF5 dataset.
    Note that tables can be empty
    """
    name = fixName(node, suffix)

    if node.text is None:
        data = numpy.array([])
    else:
        data = numpy.array( list( map(float, node.text.split( ) ) ) )
    data.shape = (rows, columns)
    h5data = parent.create_dataset( name, data=data, **h5Opts )
    addAttributes(h5data, node, index)

def addDocumentation( parent, node, index, suffix=None ):
    """
    <documentation> element stores CDATA. Convert to HDF5 text array
    """
    name = fixName( node, suffix )

    h5node = parent.create_dataset( name, data=numpy.string_(node.text) )
    addAttributes(h5node,node,index)


if __name__ == '__main__':
    args = parser.parse_args()

    xdoc = ET.parse( args.xml )

    if args.output is not None:
        h5file = args.output
    else:
        h5file = os.path.splitext(os.path.basename(args.xml))[0] + ".h5"

    if args.latestHDF and h5py.version.version_tuple > (1,4,0,''):
        # use new, more efficient HDF5 format (requires hdf v1.10 or later):
        h5 = h5py.File( h5file, "w", libver='latest' )
    else:
        h5 = h5py.File( h5file, "w" )

    if args.compress:
        h5Opts['compression'] = 'gzip'
        h5Opts['shuffle'] = True

    root = xdoc.getroot()
    addNode( h5, root, index=0 )
