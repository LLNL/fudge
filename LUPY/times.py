# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the :py:class:`Times` class that is useful for timing sections of a python code.
"""

import os, time

timeIndicesNames = { 0 : 'user', 1 : 'sys', 2 : 'children', 3 : 'childrenSys', 4 : 'wall' }

class Times :
    """
    This class is useful for timing sections of a python code. This class uses the python class :py:func:`os.times` 
    to measure time information. Five components of time are measured. They are:

        * user:         user cpu time.
        * sys:          user system time.
        * children:     total children cpu time.
        * children sys: total children system time.
        * wall:         wall time.

    All times are since :py:func:`Times.reset` was called with the :py:func:`Times.reset` method calls.

    Example usage is::

        from LUPY import times as timesModule

        time = timesModule.Tiems()
        # Some python code to time.
        print(time.toString(includeChildren=False)

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | times     | The instance return by :py:func:`os.times`.               |
    +-----------+-----------------------------------------------------------+
    """

    def __init__( self ) :

        self.reset( )

    def __repr__( self ) :
        """Calls :py:func:`self.toString` with no arguments."""

        return( self.toString( ) )

    def delta( self, reset = False ) :
        """
        Returns time information since the last time :py:func:`Times.reset` was called as a dictionary.
        The dictionary as the keys 'user', 'sys', 'children', 'childrenSys' and 'wall'.

        :param reset:       If **True**, :py:func:`Times.reset` is called after the delta time is determined.

        :returns:           A directory containing time components and their times.
        """

        current, ts = os.times( ), {}
        for index in timeIndicesNames : ts[timeIndicesNames[index]] = current[index] - self.times[index]
        if( reset ) : self.reset( )
        return( ts )

    def _delta_cpu( self, times, includeChildren = True, includeSystem = True ) :
        """
        Return the cpu component of the time form *times*. If *includeChildren* is **True** the cpu 
        time of the child processes are include in the cpu time. If *includeSystem* is **True** the 
        system time is include in the cpu time.

        This should probably be a static method since it does not use *self*.

        :param times:               The return instance from :py:func:`Times.delta`.
        :param includeChildren:     If **True** the cpu time of the child processes are include.
        :param includeSystem:       If **True** the system time is include.

        :returns:                   The requested time as a float.
        """

        indices = [ 0 ]
        if( includeChildren ) : indices.append( 2 )
        if( includeSystem ) : indices += [ index + 1 for index in indices ]
        t = 0
        for index in indices : t += times[timeIndicesNames[index]]
        return( t )

    def delta_cpu( self, reset = False, includeChildren = True, includeSystem = True ) :
        """
        Return the cpu component of the time for *self*. If *includeChildren* is **True** the cpu time of the child processes
        are include in the cpu time. If *includeSystem* is **True** the system time is include in the cpu time.

        :param reset:               If **True**, :py:func:`Times.reset` is called after the delta time is determined.
        :param includeChildren:     If **True** the cpu time of the child processes are include.
        :param includeSystem:       If **True** the system time is include.

        :returns:           The requested time as a float.
        """

        times = self.delta( reset )
        return( self._delta_cpu( times, includeChildren = includeChildren, includeSystem = includeSystem ) )

    def delta_wall( self, reset = False ) :
        """
        Return the wall component of the time for *self*.

        :param reset:       If **True**, :py:func:`Times.reset` is called after the delta time is determined.

        :returns:           The wall time as a float.
        """

        return( self.delta( reset )['wall'] )

    def reset( self ) :
        """Reset time to the current time."""

        self.times = os.times( )

    def toString( self, prefix = '', reset = False, includeChildren = True, includeSystem = True, current = True, names = True ) :
        """
        Returns a string representation of *self*.

        :param prefix:              A prefix string to add to the returned string.
        :param reset:               If **True**, :py:func:`Times.reset` is called after the delta time is determined.
        :param includeChildren:     If **True** the cpu time of the child processes are include.
        :param includeSystem:       If **True** the system time is include.
        :param current:             If **True**, the current date/time information is added to the returned string.
        :param names:               If **True**, the names of each component is added to the returned string.

        :returns:                   A string representing the time information.
        """

        deltas = self.delta( reset = reset )
        cpu = self._delta_cpu( deltas, includeChildren = includeChildren, includeSystem = includeSystem )
        currentTime = ''
        if( current ) : currentTime = ' on %s' % time.ctime( )
        if( names ) : return( '%sdelta times: cpu = %.3f s, wall = %.5g s%s' % ( prefix, cpu, deltas['wall'], currentTime ) )
        return( '%s%.3f s, %.5g s%s' % ( prefix, cpu, deltas['wall'], currentTime ) )

def timeCode(command, globals, locals):
    """
    This function uses a :py:class:`Times` instance to time how long *command* takes to execute.

    Example usage::

        import math
        timesModule.timeCode("z = math.sqrt( 3**2 + 4**2 ); print(z)", globals(), locals())

    :param command:         A string of the python code to execute.
    :param globals:         A global disctionary to pass to the exec function, normally should be globals().
    :param locals:          A local disctionary to pass to the exec function, normally should be locals().
    """

    times = Times()
    exec(command, globals, locals)
    print(times)

if( __name__ == '__main__' ) :

    t = Times( )
    print(t)
    time.sleep( 1.4 )
    print(t.toString( ))
    print(t.toString( current = False ))
    for i in range( 1000 * 1000 ) : y = i**2
    print(t.toString( prefix = '# ', current = False ))
    print(t.toString( prefix = '# ', current = True ))
    print(t.delta_cpu( ))
    print(t.delta_wall( ))
    print(t.toString( prefix = '# ', current = True ))
    print(t.toString( prefix = '# ', current = True, includeChildren = False ))
    print(t.toString( prefix = '# ', current = True, includeChildren = False, includeSystem = False ))
