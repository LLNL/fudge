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
