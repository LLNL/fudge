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

import os, time
timeIndicesNames = { 0 : 'user', 1 : 'sys', 2 : 'children', 3 : 'childrenSys', 4 : 'wall' }

class times :

    def __init__( self ) :

        self.reset( )

    def __repr__( self ) :

        return( self.toString( ) )

    def delta( self, reset = False ) :

        current, ts = os.times( ), {}
        for index in timeIndicesNames : ts[timeIndicesNames[index]] = current[index] - self.times[index]
        if( reset ) : self.reset( )
        return( ts )

    def _delta_cpu( self, times, includeChildren = True, includeSystem = True ) :

        indices = [ 0 ]
        if( includeChildren ) : indices.append( 2 )
        if( includeSystem ) : indices += [ index + 1 for index in indices ]
        t = 0
        for index in indices : t += times[timeIndicesNames[index]]
        return( t )

    def delta_cpu( self, reset = False, includeChildren = True, includeSystem = True ) :

        times = self.delta( reset )
        return( self._delta_cpu( times, includeChildren = includeChildren, includeSystem = includeSystem ) )

    def delta_wall( self, reset = False ) :

        return( self.delta( reset )['wall'] )

    def reset( self ) :

        self.times = os.times( )

    def toString( self, prefix = '', reset = False, includeChildren = True, includeSystem = True, current = True ) :

        deltas = self.delta( reset = reset )
        cpu = self._delta_cpu( deltas, includeChildren = includeChildren, includeSystem = includeSystem )
        currentTime = ''
        if( current ) : currentTime = ' on %s' % time.ctime( )
        return( '%sdelta times: cpu = %s s, wall = %.5g s%s' % ( prefix, cpu, deltas['wall'], currentTime ) )

if( __name__ == '__main__' ) :

    t = times( )
    print t
    time.sleep( 1.4 )
    print t.toString( )
    print t.toString( current = False )
    for i in xrange( 1000 * 1000 ) : y = i**2
    print t.toString( prefix = '# ', current = False )
    print t.toString( prefix = '# ', current = True )
    print t.delta_cpu( )
    print t.delta_wall( )
    os.system( 'python timesTest.py' )
    print t.toString( prefix = '# ', current = True )
    print t.toString( prefix = '# ', current = True, includeChildren = False )
    print t.toString( prefix = '# ', current = True, includeChildren = False, includeSystem = False )
