# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import parsers as PR
from pqu.PQU import PQU

def a( s ) :

    print()
    print("<%s>" % s.strip( ))
    try :
        value, unit, uncertainty = PR.parsePQUString( s )
#        print('<%s> <%s> <%s>' % ( value.toString( ), uncertainty.toString( prefix = ' ' ), unit ))
        unitStr = str( unit )
        if( len( unitStr ) > 0 ) : unitStr = ' ' + unitStr
#        print('%s%s%s' % ( value.toString( ), uncertainty.toString( prefix = ' ' ), unitStr ))
        pqu = PQU( s )
        print('<%s>' % pqu)
    except :
        print('ERROR')
        raise

a( '0.0 +/- 0.05 Ang' )
a( '0.0 +/- 0.050 Ang' )
a( '0.0 +/- 0.05012 Ang' )
a( '0.0 +/- 0.05000 Ang' )
a( '-0.0 +/- 0.05000 Ang' )
a( '  0.135    +/- 1.2% eV/mm  ' )
a( '  135    +/- 1.2% eV/mm  ' )
a( '  11.e3 ' )
a( '  11.009e3 ' )
a( '  .009e3 ' )
a( '  .135 \t ' )
a( '  .135 \t ' )
a( '  0.135  ' )
a( '  .135 \t\t eV  ' )
a( '  0.135    +/- 1.2% eV/mm  ' )
a( '  -0.135+/- 1.2% eV/mm  ' )
a( '  +0.135    +/- 1.2% eV/mm  ' )
a( '  0.135+/-1.2 ' )
a( '  0.135 +/- 1.2e-2 eV/mm  ' )
a( '  0.135(12) eV/mm  ' )
a( '  0.135 eV / mm  ' )
a( '  -0.135e-12 +/- 1.2% eV/mm  ' )
a( '  -0.135e-12 +/- 1.2e-14 eV/mm  ' )
a( '  -0.135e-13 +/- 1.2e-14 eV/mm  ' )
a( '  -0.135e-12(12) eV/mm  ' )
a( '  -0.135e-12(4) eV/mm  ' )
a( '  -0.135e+12(12) eV/mm  ' )
a( '  -0.135e+12(4) eV/mm  ' )
a( '  13.5  ' )
a( '  13.5%  ' )
a( '0.0 +/- 0.0555 Ang' )
a( '  -0.135e-12(12)eV/mm  ' )
