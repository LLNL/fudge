# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module defines 'pointwiseXY' in pure python. Users will generally prefer to use 'pointwiseXY_C',
which is a faster and more capable C implementation.
If the C version is not available, fudge will revert to this version.

cmattoon, July 2011
"""
from endl2dmathClasses import endl2dmath

__metaclass__ = type

class pointwiseXY( endl2dmath ):
    def __init__(self, initialSize=0, overflowSize=0, data=None, dataForm="xys", accuracy=None, 
            biSectionMax=4, interpolation=None, infill=True, safeDivide=False):

        if( dataForm.lower() == 'xsandys' ) :
            data = zip(*data)
        if( type( interpolation ) is str ) :
            if( 'flat' in interpolation ) : interpolation = 'linear,linear'     # This is a kludge
            interpolation = { 'linear,linear' : 0, 'log,linear' : 1, 'linear,log' : 2, 'log,log' : 3 }[interpolation]
            endl2dmath.__init__(self, data=data, interpolation=interpolation)
        else:
            endl2dmath.__init__(self, data=data)

        self.accuracy = accuracy
        self.biSectionMax = biSectionMax
        self.infill = infill
        self.safeDivide = safeDivide

    # override some endl2dmath functions for compatibility with pointwiseXY_C:
    def __add__(self, other):
        return pointwiseXY( data = endl2dmath.__add__(self,other) )

    def __mul__(self, other):
        return pointwiseXY( data = endl2dmath.__mul__(self,other) )

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        if type(value) in (int,float):
            self.data[index][1] = value
        else:
            self.data[index] = value

    def getAccuracy(self):
        return self.accuracy

    def getBiSectionMax(self):
        return self.biSectionMax

    def getInfill(self):
        return self.infill

    def getSafeDivide(self):
        return self.safeDivide

    def thin(self, accuracy):
        return pointwiseXY( data = endl2dmath.thin( self, accuracy ) )

    def union(self, other, fillWithSelf, trim):
        d = endl2dmath.union( self, other )
        return pointwiseXY( data=d.data, interpolation=self.interpolation )

    def xSlice(self, xMin, xMax, fill=None, dullEps=None):
        return self
        d = endl2dmath.slicex( self, xMin, xMax )
        return pointwiseXY( data=d.data, interpolation=self.interpolation )

    def copyDataToXYs(self, xScale=None, yScale=None):
        if xScale or yScale:
            xScale = xScale or 1; yScale = yScale or 1
            return [[xScale*x,yScale*y] for (x,y) in self.data]
        import copy
        return copy.deepcopy( self.data )

    def copy(self):
        return pointwiseXY( data=self.data, accuracy=self.accuracy, biSectionMax=self.biSectionMax, 
                interpolation=self.interpolation, infill=self.infill, safeDivide=self.safeDivide )

    def changeInterpolation( self, interpolation, accuracy = None, lowerEps = 0, upperEps = 0 ) :
        if interpolation != "linear,linear":
            raise Exception("Only linear interpolation currently supported")
        return self.toInterpolation(0, accuracy)

    def __str__( self ) :

        return( endl2dmath.toString( self ) )
