from pqu.PQU import PQU

pqu1 = PQU( "1.23", "m", "12%" )
ref = pqu1
pqu1 += "3.21 m"
print pqu1
print ref

pqu1 -= "3.21 m"
print pqu1
print ref

pqu1 *= "2.2 kg"
print pqu1
print ref

pqu1 /= "2.2 kg"
print pqu1
print ref
