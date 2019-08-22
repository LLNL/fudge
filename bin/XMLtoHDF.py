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
Translate GND/XML to HDF5. Each xml element becomes a group in HDF, except a few elements that become
HDF datasets. Attributes are translated to HDF5 metadata
"""
import sys
import os
import numpy
from collections import Counter

from xml.etree import cElementTree as parser
import h5py

def fixName( node, suffix=None ):
    name = node.tag
    if suffix is not None:
        name += '_%s' % str(suffix)
    return name

def addAttributes( h5node, node, index ):
    for key,val in node.items():
        h5node.attrs.create(str(key),str(val))
    h5node.attrs.create("_xmltag", node.tag)    # save original tag for translating back to XML
    h5node.attrs.create("_xmlindex", index)  # hdf5 doesn't preserve group order, so preserve it here

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

        child_names = Counter( [ child.tag for child in node.getchildren() ] )
        for key,value in child_names.items():
            if value == 1: child_names[key] = None
            else: child_names[key] = 0

        for idx,child in enumerate(node.getchildren()):
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
        print ("addNode backtrace: node=%s, name=%s" % (node,name))
        raise

def addValues( parent, node, index, suffix=None, isTwoDimensional=False ):
    """
    <values> elements in GND store a list of numbers. Convert to HDF5 dataset
    """
    name = fixName(node, suffix)

    data = numpy.array( map(float, node.text.split()) )
    if isTwoDimensional:
        data.shape = (len(data)//2,2)
    h5data = parent.create_dataset( name, data=data )
    addAttributes(h5data, node, index)

def addTableData( parent, node, index, rows, columns, suffix=None ):
    """
    <data> elements in GND store a table. Convert to HDF5 dataset.
    Note that tables can be empty
    """
    name = fixName(node, suffix)

    data = node.text
    if data is None:
        data = numpy.array([])
    else:
        data = numpy.array( map(float, node.text.split()) )
    data.shape = (rows, columns)
    h5data = parent.create_dataset( name, data=data )
    addAttributes(h5data, node, index)

def addDocumentation( parent, node, index, suffix=None ):
    """
    <documentation> element stores CDATA. Convert to HDF5 text array
    """
    name = fixName( node, suffix )

    h5node = parent.create_dataset( name, data=numpy.array(node.text) )
    addAttributes(h5node,node,index)


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
    addNode( h5, root, index=0 )
