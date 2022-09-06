#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Translate GNDS/HDF5 to GNDS/XML
"""
import sys
import os

from xml.etree import cElementTree as etree
import numpy
import h5py

from pqu import PQU

def addNode( parent, node ):
    """
    recursively copy HDF data (one node at a time) to xml

    :param parent: xml element
    :param node:   HDF group or dataset to be added to xml
    """
    attrs = dict( node.attrs )
    for key in attrs:
        if type(attrs[key]) is numpy.bytes_:
            attrs[key] = attrs[key].decode("utf8")
    name = attrs.pop( "_xmltag" )

    try:
        xmlnode = etree.Element( name )
        for key,val in attrs.items():
            if key.startswith("_xml"): continue
            xmlnode.set( str(key), str(val) )
        if isinstance( parent, etree.ElementTree ): parent._setroot( xmlnode )
        else: parent.append( xmlnode )
        newParent = xmlnode

        if isinstance( node, h5py.Dataset ):
            addDataset( xmlnode, node )
        else:
            for child in sorted( node.values(), key=lambda tmp: tmp.attrs["_xmlindex"] ):
                addNode( newParent, child )
    except:
        print( "addNode backtrace: node=%s, name=%s" % ( node, name ) )
        raise

def addDataset( parent, dataset ):
    """
    Convert HDF5 dataset into text

    :param parent:  xml element corresponding to the dataset
    :param dataset: HDF5 dataset to be converted
    """
    if str(dataset.dtype).startswith("|S"):
        parent.text = dataset[()].decode("utf8")
    else:
        parent.text = ' '.join( map( PQU.toShortestString, dataset[()].flatten() ) )


if __name__ == '__main__':
    if len(sys.argv)<2:
        print( "Usage: python HDFtoXML.py <HDF_filename>" )
        sys.exit(1)

    h5file = sys.argv[1]
    h5 = h5py.File( h5file, "r" )

    xdoc = etree.ElementTree()

    root = list(h5.values())[0]
    addNode( xdoc, root )

    xdoc.write( os.path.splitext(h5file)[0] + ".hdf5.xml" )
