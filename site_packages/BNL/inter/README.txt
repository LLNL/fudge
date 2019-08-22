INTERrogate an evaluation!
==========================

20 Jun 2017
David Brown

Super fancy replacement for legacy INTER code. This does everything that legacy INTER does 
(use the "--report legacy" option), but so much more.  Invoke thusly to get help:

$ fudge/site_packages/bin/inter.py -h

To add your own spectra for flux averaging, see fudge/site_packages/BNL/inter/example_spectra.json
for examples.  Pay special attention to the format field.  The energy is denoted in the “x” column 
or in “xmin” and “xmax” columns.  The flux in the “y” column.  I assume the flux is grouped, so there  
must be one more energy point than flux points.  If “xmin” and “xmax” are given, you are set.  If 
only “x” is given, you’ll be short by one group boundary.  Don’t worry, I’ll make something up ;).

Note, not all of the resonance reports behave yet.






First results!
==============

Nov 13, 2015
David Brown

Here are first real results from inter.py, compared to the Atlas, legacy INTER and 
from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.  Overall all we are in 
good-excellent agreement.

These calcs were done on the ENDF/B-VII.1 n-092_U_233.endf file

From 0.000010 to 600.000000 eV, reconstructing using Reich_Moore (LRF=3)
12939 points were added (for total of 188215) to achieve tolerance of 0.1%
Unresolved from 600.000000 to 40000.000000 eV
Skipping unresolved: for self shielding only
*************** summary reactions report ***************
name                                                 MT          Q Threshold  
unit                                                            eV        eV   
n + U233                                              2          0         0   
n[multiplicity:'energyDependent', emissionMode:...   18  191040000         0   
U234 + gamma                                        102    6844200         0   
total                                                 1          0         0   

name                                                               RI      Boris                Atlas          INTER (INTER cuts off
unit                                                                b          b                    b              b  integration at 
n + U233                                            169.6691534332043   1.696E+2                         1.43696E+02  1e5 eV)
n[multiplicity:'energyDependent', emissionMode:...  775.5378882051494   7.755E+2    7.750E+2±1.700E+1    7.64818E+02
U234 + gamma                                        141.0467349894915   1.411E+2    1.380E+2±6.000E+0    1.40571E+02
total                                               1091.577791339657                                    1.04923E+03

name                                                         14.2 MeV   INTER(14.0 MeV)
unit                                                                b                 b
n + U233                                                    2.8361696       2.81722E+00
n[multiplicity:'energyDependent', emissionMode:...          2.3500125       2.33500E+00
U234 + gamma                                              1.245229e-3       1.25780E-03
total                                               5.778268363636363       5.75726E+00

name                                                         2200 m/s      Boris               Atlas       Standards         INTER 
unit                                                                b          b                   b               b             b
n + U233                                            12.14812890310712   1.220E+1                                       1.21512E+01
n[multiplicity:'energyDependent', emissionMode:...  532.0075437829548   5.314E+2    5.291E+2±1.200E0    531.22±1.328   5.31215E+02
U234 + gamma                                        45.27169009237033   4.526E+1    4.550E+1±7.000E-1   45.56±0.6834   4.52376E+01
total                                               589.4610140892717                                                  5.88604E+02

name                                                   Westcott factor      Boris                Atlas     INTER (INTER cuts off 
unit                                                                                                              integration at 
n + U233                                             1.125691948727054   1.126E+0                        1.12583  10 eV and 
n[multiplicity:'energyDependent', emissionMode:...  0.9937414875835585   9.966E-1    9.955E-1±1.500E-3   0.99652  ignores recoil)
U234 + gamma                                         1.030197591211457   1.032E+0                        1.03203
total                                               0.9992209816705242                                   1.00191

name                                                     MACS(30. keV)               Boris
unit                                                                 b                   b
n + U233                                             11.46808147500756   
n[multiplicity:'energyDependent', emissionMode:...   3.000088363913941   3.001E+0±4.987E-2
U234 + gamma                                        0.4239508275461996   4.240E-1±6.106E-2
total                                                15.04574318099927   

name                                                    MACS(1420. keV)               Boris
unit                                                                  b                   b
n + U233                                              4.751480983374386  
n[multiplicity:'energyDependent', emissionMode:...    2.096723739994021   2.097E+0±6.232E-3
U234 + gamma                                        0.04900796580060433   4.903E-2±8.460E-3
total                                                 8.383778009408696  

name                                                     ARR(30. keV)  
unit                                                      cm**3/s/mol  
n + U233                                            1658106431.893352  
n[multiplicity:'energyDependent', emissionMode:...  433766172.9466244    4.339E+8±7.210E-6
U234 + gamma                                        61296703.85520167    6.131E+7±8.822E+6
total                                               2175380737.867837  

*************** addendum ***************
ALF = 0.08509595516344783 <-- Compare with EXFOR SUBENT #10013002's value of 0.0899 +/- 4.0E-4
ETA = 2.300994661457297   <-- Compare with EXFOR SUBENT #40426003's value of 2.297