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

import sys
sys.path.insert( 0, '../../' )

from pqu.PQU import pqu_float, PQU

f1 = pqu_float( 1.2, 2, True )
f2 = pqu_float( 11, 2, True )
f = f1 * 4
print f1, f
f = 5 * f1
print f1, f
print f1 * f2, f2 * f1

print f1.info( significantDigits = 15 )
f = f1 / 5
print f, f.info( significantDigits = 15 )
f = 5 / f1
print f, f.info( significantDigits = 15 )

print
f1 = pqu_float(  1.234567e10, 7, False )
f2 = pqu_float( -1.234558e10, 6, False )
f3 = pqu_float( -7123e4, 6, False )
print f1.info( significantDigits = 15 )
print f2.info( significantDigits = 15 )
f = f1 + f2
print f, f.info( significantDigits = 15 )
f = f2 + f1
print f, f.info( significantDigits = 15 )
f = f + f3
print f, f.info( significantDigits = 15 )
f = 1.2345e10 - f1
print f, f.info( significantDigits = 15 )

print f1 == 1.234567e10
print f1 == 1.23456e10
print f1 != 1.23456e10
print f1 <= 1.23456e10
print f1 < 1.23456e10
print f1 >= 1.23456e10
print f1 > 1.23456e10


print
p1 = PQU( '0.0 +/- 0.05 Ang' )
p2 = PQU( '  0.135    +/- 1.2% eV/mm  ' )
p3 = PQU( '  11.e4 ' )
p4 = PQU( '13.5 +/- 3.2% eV/mm  ' )
p5 = PQU( '12.1 +/- 12.21% eV/mm  ' )
p6 = PQU( '12.1 +/- 1.2% eV/mm  ' )

print p1
print p1.info( significantDigits = 15 )
print p2
print p2.info( significantDigits = 15 )
print p3
print p3.info( significantDigits = 15 )
print p4
print p4.info( significantDigits = 15 )
print p5
print p5.info( significantDigits = 15 )
print p6
print p6.info( significantDigits = 15 )

def operator( oper, p1, p2 ) :

    def operator2( oper, p1, p2 ) :

        print 
        print '===== ( %s ) %s ( %s ) =====' % ( p1, oper, p2 )
        try :
            p = eval( 'p1 %s p2' % oper, { 'p1' : p1, 'p2' : p2 } )
            print p
            print p.info( significantDigits = 15 )
        except ZeroDivisionError : 
            print 'ERROR 1/0: ( %s ) %s ( %s ) ' % ( p1, oper, p2 )
        except :
            raise

    operator2( oper, p1, p2 )
    operator2( oper, p2, p1 )

operator( '*', p1, p2 )
operator( '*', p3, p2 )
operator( '*', p2, '2.2 +/- 0.05 Ang' )

operator( '/', p2, '2.2 +/- 0.05 Ang' )
operator( '/', p1, p2 )
operator( '/', p3, p2 )
operator( '/', p1, f1 )
operator( '/', p2, f1 )
operator( '/', p3, f1 )

operator( '+', p2, p4 )
operator( '+', p5, p4 )
operator( '+', p6, p4 )
operator( '-', p2, p4 )
operator( '-', p5, p4 )
operator( '-', p6, p4 )
