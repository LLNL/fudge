# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains miscellaneous routines for supporting the endl1dmath class from module endl1dmathClasses.

Module's global variables::
    endl1d_repr_xFormat : Default format string used for converting an element of the data of
                          an endl1dmath object to a string.
"""

import types
import endl1dmathClasses
import fudgemath

endl1d_repr_xFormat = "%14.7e"

def get1dmathData( object, callingRoutine, msg ) :
    """If object is a python list, then it is returned. Else, if object is a sub-class
    of endl1dmath its data is return otherwise a raise is executed."""

    if( type( object ) == type( [] ) ) :
        data = object
    else :
        data = valid1dClassType( object, callingRoutine, msg ).data
    return( data )

def valid1dClassType( object, callingRoutine = "", msg = "" ) :
    """Returns the first argument, object, if it is a subclass of the endl1dmath class; else, triggers a raise."""

    if( isinstance( object, endl1dmathClasses.endl1dmath ) ) : return( object )
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check1dData( data ) :
    """Checks that data is a python list of numbers. Triggers a raise if it is not.  Data may be a
python list or an instance that is a subclass of endl1dmath."""

    data = get1dmathData( data, "check1dData", "" )
    i = 0
    for p in data :
        if( not( fudgemath.isNumber( p ) ) ) : raise Exception( "\nError in check1dData: data at index %d is not a number." % i )
        i += 1

def list1dToXMLPointwiseString( data, varString = '', indent = '', floatFormat = '%16.9e' ) :
    iscomplex = any( [isinstance(a,complex) for a in data] )
    s = [ '%s<xData type="1d.x" length="%d"' % ( varString, len( data ) ) ]
    if iscomplex: s[0] += ' complex="True"'
    s[0] += '>'
    if iscomplex:
        d = [ ' %s+%sj' % (floatFormat,floatFormat) % (x.real,x.imag) for x in data ]
    else:
        d = [ ' %s' % ( floatFormat % x ) for x in data ]
    s += d
    s.append( '</xData>' )
    return( ''.join( s ) )

def list1dToXMLEqualProbableBins1dString( data, indent = '', floatFormat = '%16.9e' ) :
    
    s = [ '<xData type="1d.x" length="%d">' % len( data ) ]
    s += [ ' %s' % ( floatFormat % x ) for x in data ]
    s.append( '</xData>' ) 
    return( ''.join( s ) )

def endl1dToXMLGroup1dStringList( self, indent = '', floatFormat = '%16.9e' ) :

# ?????????? This should no long be used.

    from fudge import gnd
    form = self.form
    data = self.data
    return( [ "%s<%s>%s%s</%s>" % ( indent, form, datasVariableInfoToXML( data, indent+'  ' ), 
        list1dToXMLEqualProbableBins1dString( data ), form ) ] )

def datasVariableInfoToXML( data, indent='' ) :

    if( not hasattr( data, 'variables' ) ) : return( '' )
    return( variableInfoToXML( data.variables, indent ) )

def variableInfoToXML( variables, indent='' ) :

    varList = []
    for index in variables :
        fStr = ''
        for name in variables[index] : fStr += ' %s="%s"' % ( name, variables[index][name] )
        varList.append( '<variable index="%d"%s/>' % ( index, fStr ) )
    if len(varList) > 2:
        return ''.join(['\n'+indent+var for var in varList])
    return( ''.join(varList) )
