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

import unittest
import os
from fudge.legacy.endl import endlProject, __path__

class testEndlProject( unittest.TestCase ): 
    
    def testInit( self ): 
        e = endlProject( os.sep.join( __path__ + [ 'test', 'testdb' ] ) )
        self.assertEqual( e.ZAList(), [ 'za001001' ] )

    def testReadZA( self ): 
        e = endlProject( os.sep.join( __path__ + [ 'test', 'testdb' ] ) )
        self.assertEqual( e.ZAList(), [ 'za001001' ] )
        za = e.readZA( 1001 )
        za.read()
        self.assertEqual( [ x.C for x in za.findDatas( I=0 ) ], [ 1, 10, 46 ] )

if __name__== "__main__":
    unittest.main()
