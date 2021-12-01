# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import datetime

__metaclass__ = type

class Resolution :

    date = 'date'
    time = 'time'
    undefined = 'undefined'

    allowed = ( date, time, undefined )

class Date :
    """
    The class supports a date or datetime in one of 3 forms. The forms are

        -) date only (e.g., "2021-01-21"),
        -) date and time (e.g., "2021-01-21T08:42:04") or
        -) date and time with UTC offset (e.g., "2021-01-21T08:42:04-08:00").
    """

    def __init__( self, _datetime = None, resolution = Resolution.date, UTC = False ) :

        if( resolution not in Resolution.allowed ) : raise TypeError( 'Invalid resolution instance.' )
        self.__resolution = resolution

        self.__datetime = None
        if( resolution == Resolution.undefined ) : return

        if( _datetime is None ) :
            _datetime = datetime.datetime.now( ).replace( microsecond = 0 )
            if( resolution == Resolution.date ) :
                _datetime = _datetime.replace( hour = 0, minute = 0, second = 0, microsecond = 0 )
            elif( ( sys.version_info.major == 3 ) and UTC ) :
                _datetime = _datetime.astimezone( )

        if( not( isinstance( _datetime, datetime.datetime ) ) ) : raise TypeError( 'Invalid _datetime instance.' )

        self.__datetime = _datetime

    def __str__( self ) :

        if( self.__resolution == Resolution.undefined ) : return( '' )

        stringDate = self.__datetime.isoformat( sep = 'T', timespec = 'seconds' )
        if( self.__resolution == Resolution.date ) : stringDate = stringDate.split( 'T' )[0]

        return( stringDate )

    @property
    def date( self ) :

        return( self.__datetime )

    @property
    def resolution( self ) :

        return( self.__resolution )

    def asXML_attribute( self, name = 'date' ) :

        if( self.resolution == Resolution.undefined ) : return( '' )
        return( ' %s="%s"' % ( name, self ) )

    @staticmethod
    def parse( stringDate ) :

        resolution = Resolution.date
        datetime1 = None

        if( len( stringDate ) == 0 ) :
            resolution = Resolution.undefined
        elif( len( stringDate ) == 10 ) :
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%d" )
        elif( len( stringDate ) == 19 ) :
            resolution = Resolution.time
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%dT%H:%M:%S" )
        else :
            resolution = Resolution.time
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%dT%H:%M:%S%z" )

        return( Date( datetime1, resolution = resolution ) )

    @staticmethod
    def fromYearMonthDay( year, month, day ) :

        return( Date( datetime.datetime( year, month, day ) ) )

def raiseIfNotDate(date):
    
    if not isinstance(date, Date): raise TypeError('Invalid date instance.')

    return date

def raiseIfNotPythonDateTime(_datetime):

    if not isinstance(_datetime,datetime.datetime): raise TypeError('Invalid python datetime instance.')

    return( Date(_datetime) )
