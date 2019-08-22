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

try   :
    from xml.etree.cElementTree import ElementTree, tostring, SubElement, Element
except ImportError   :
    from ElementTree import ElementTree, tostring, SubElement, Element

class XML2PY( object ) :
    def __init__( self ) :
        self.etree = ElementTree( )
    def __repr__( self ): return repr( self.__dict__ )
    def __setattr__( self, attr, val ) :
        if attr not in ['xml', 'etree'] :
            if attr not in [e.tag.lower( ) for e in self.xml.getchildren( )] :
                self.xml.set( attr, str( val ) )
        return object.__setattr__( self, attr, val )
    def newXMLElement( self, attr, val={}, index=0 ) :
        element = Element( attr.capitalize( ), val )
        element.text = None
        element.tail = '\n    '
        self.xml.insert( index, element )
        setattr( self, element.tag.lower( ), self._get_py( element ) )
    def read( self, xmlfile ) :
        self.etree.parse( xmlfile )
        update( self, self.etree.getroot( ) )
        for element in self.etree.getroot( ) :
            setattr( self, element.tag.lower( ), self._get_py( element ) )
    @staticmethod
    def _get_py( xml ) :
        x = XMLElement( )
        update( x, xml )
        for element in xml :
            x.append( XML2PY._get_py( element ) )
        return x
    def toxml( self ) :
        return tostring( self.xml )
    def write( self, xmlfile ) :
        self.etree.write( xmlfile )
    def fixIndentation( self ) :
        root = self.etree.getroot( )
        root.text = '\n    '
        root.tail = None
        self.textTail( root.getchildren( ) )
    @staticmethod
    def textTail( elements, indentation=0 ) :
        def space( indentation ) : return ''.join( [' '] * ( indentation * 4 ) )
        for element in elements :
            if len( element ) :
                element.text = '\n' + space( 2 + indentation )
                XML2PY.textTail( element, indentation + 1 )
            else :
                element.text = None
            element.tail = '\n' + space( 1 + indentation )
        if len( element ) :
            element.text = '\n' + space( 2 + indentation )
        else :
            element.text = None
        element.tail = '\n' + space( indentation )

class XMLElement( list ) :
    def __init__( self, l=[] ) :
        list.__init__( self, l )
    def __repr__( self )  :
        return repr( dict( [( key, val ) for key, val in self.__dict__.items( ) if key != 'xml'] ) )
    def __setattr__( self, attr, val ) :
        if attr != 'xml' :
            if self.xml.tag == 'Level' and attr == 'energy' :
                self.xml.set( attr, '%.6f' % val )
            else :
                self.xml.set( attr, str( val ) )
        return list.__setattr__( self, attr, val )
    def __delitem__( self, index ) :
        self.xml.remove( self.xml.getchildren( )[index] )
        return list.__delitem__( self, index )
    def remove( self, item ) :
        index = self.index( item )
        return self.__delitem__( index )

def evaluate( item ) :
    try: return eval( item )
    except: return str( item )

def update( obj, element ) :
    if hasattr( element, 'keys' ) :
        obj.xml = element
        for attr in element.keys( ) :
            setattr( obj, attr, evaluate( element.get( attr ) ) )

if __name__=='__main__':
    import fudge,os
    x=XML2PY()
    x.read(os.path.join(fudge.fudgeDefaults.XENSL_DATABASE_DIR,'za036072.xml'))
    x.newXMLElement('isomers')
    x.isomers.m0='e0'
    print x.toxml()
    x.fixIndentation()
    print x.toxml()

