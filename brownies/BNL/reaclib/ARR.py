#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt



#read the endf file in here
#from brownies.legacy.converting.endfFileToGNDS import endfFileToGNDS
#input_file='/Users/laurenosojnak/Desktop/ENDF-B-VIII.0/neutrons/n-079_Au_197.endf'
#resultsDict = endfFileToGNDS( input_file, toStdOut=True, skipBadData=True )

#def check_is_cross_section( x ):
#    if not isinstance( x, component ):
#        raise TypeError( "Not instance of fudge.reactionData.crossSection.component" )


#def computeReaclibARR(xs, T):
#    computeAstrophysicalReactionRate(xs, T)
#    Lambda=np.exp(a_zero+ thesum + a_six*np.log(T))
#    thesum=0
#    while i in [1,2,3,4,5]:
#        thesum+=(a_i*(T**(((2*i)-5)/3)))
#        return thesum

#convert Kelvn to keV
Kelvin=0.00000008617328149741
Temp_inkeV=[]
Temperature=[.1, .2, .3, .4, .5, 1, 1.5, 2, 10]
for iT in Temperature:
    r=(Kelvin*10**9)*(iT)
    Temp_inkeV.append(r)
b_matrix=[]
# for n+Fe55 ->Fe56
# for T9=.1 -> 8.617328149741 keV -> ARR= 20324269.946729 cm**3/s/mol
# for T9=.2 -> 17.234656299482 keV, ARR= 16076263.1081873 cm**3/s/mol
# for T9=.3 -> 25.851984449223 keV, ARR= 13908572.9977194 cm**3/s/mol
# for T9=.4 -> 34.469312598964 keV, ARR= 12652186.6350402 cm**3/s/mol
# for T9=.5 -> 43.086640748705 keV, ARR=  11838149.1880259 cm**3/s/mol
# for T9=1 -> 86.17328149741 keV, ARR= 10019033.7916955 cm**3/s/mol
# for T9=1.5 -> 129.259922246115 keV, ARR=  9311391.12427268 cm**3/s/mol
# for T9=2 -> 172.34656299482 keV, ARR=  8925449.81071796 cm**3/s/mol
# for T9=10 -> 861.7328149741 keV, ARR=  6715355.09181261 cm**3/s/mol
Fe55_ARR=[20324269.946729,16076263.1081873 , 13908572.9977194, 12652186.6350402,  11838149.1880259,10019033.7916955 , 9311391.12427268, 8925449.81071796, 6715355.09181261  ]
# for n+Au197 -> Au198
# for T9=.1 -> 8.617328149741 keV -> ARR= 108359277.292297 cm**3/s/mol
# for T9=.2 -> 17.234656299482 keV, ARR= 95124164.2944291 cm**3/s/mol
# for T9=.3 -> 25.851984449223 keV, ARR= 90524017.6799109 cm**3/s/mol
# for T9=.4 -> 34.469312598964 keV, ARR= 88587427.1385287 cm**3/s/mol
# for T9=.5 -> 43.086640748705 keV, ARR= 87814745.0164856 cm**3/s/mol
# for T9=1 -> 86.17328149741 keV, ARR= 88060941.0447786 cm**3/s/mol
# for T9=1.5 -> 129.259922246115 keV, ARR= 87912244.8950413 cm**3/s/mol
# for T9=2 -> 172.34656299482 keV, ARR= 86778030.3983979 cm**3/s/mol
# for T9=10 -> 861.7328149741 keV, ARR= 67611130.5968464 cm**3/s/mol
Au197_ARR=[108359277.292297,95124164.2944291,90524017.6799109,88587427.1385287,87814745.0164856,88060941.0447786,87912244.8950413, 86778030.3983979, 67611130.5968464]

# for n+Fe56 -> Fe57
# for T9=.1 -> 8.617328149741 keV -> ARR=   998726.266197593 cm**3/s/mol
# for T9=.2 -> 17.234656299482 keV, ARR=  1544758.74041983 cm**3/s/mol
# for T9=.3 -> 25.851984449223 keV, ARR=  1871147.6025848 cm**3/s/mol
# for T9=.4 -> 34.469312598964 keV, ARR=  2053634.66089578 cm**3/s/mol
# for T9=.5 -> 43.086640748705 keV, ARR=  2157733.52889866 cm**3/s/mol
# for T9=1 -> 86.17328149741 keV, ARR= 2329157.77560464 cm**3/s/mol
# for T9=1.5 -> 129.259922246115 keV, ARR= 2410355.77382783 cm**3/s/mol
# for T9=2 -> 172.34656299482 keV, ARR= 2490303.81666388 cm**3/s/mol
# for T9=10 -> 861.7328149741 keV, ARR= 2526008.85828584 cm**3/s/mol
ARR_format=[]
Fe56_ARR=[998726.266197593,1544758.74041983,1871147.6025848,2053634.66089578,2157733.52889866,2329157.77560464, 2410355.77382783,2490303.81666388,2526008.85828584]
# for n +p -> d
# for T9=.1 -> 8.617328149741 keV -> ARR= 38131.7916485332 cm**3/s/mol
# for T9=.2 -> 17.234656299482 keV, ARR= 34394.1219690182  cm**3/s/mol
# for T9=.3 -> 25.851984449223 keV, ARR=  31799.3968665715 cm**3/s/mol
# for T9=.4 -> 34.469312598964 keV, ARR=  29905.031665433 cm**3/s/mol
# for T9=.5 -> 43.086640748705 keV, ARR= 28482.0059060774  cm**3/s/mol
# for T9=1 -> 86.17328149741 keV, ARR=  25070.3690553475 cm**3/s/mol
# for T9=1.5 -> 129.259922246115 keV, ARR= 24457.5461388397 cm**3/s/mol
# for T9=2 -> 172.34656299482 keV, ARR=  24961.4400818657 cm**3/s/mol
# for T9=10 -> 861.7328149741 keV, ARR= 44601.1275895236 cm**3/s/mol
H_ARR=[38131.7916485332,34394.1219690182, 31799.3968665715, 29905.031665433, 28482.0059060774, 25070.3690553475,24457.5461388397, 24961.4400818657, 44601.1275895236]
for iR in H_ARR:
    number=iR*pow(10,-6)
    ARR_format.append(number)
plt.plot(Temperature, ARR_format, label='ENDF/B-VIII.0')

for iARR in Fe56_ARR:
    b_matrix_row=[np.log(iARR)]
    b_matrix.append(b_matrix_row)
b=np.array(b_matrix)
a_matrix=[]
for iT in Temperature:
    a_matrix_row=[1, pow(iT,-1), pow(iT,-1/3), pow(iT,1/3), pow(iT,1), pow(iT,5/3), np.log(iT)]
    a_matrix.append(a_matrix_row)
a=np.array(a_matrix)

x=np.linalg.lstsq(a,b, rcond=None)
print( x[0] )

#myplot
#myplot=plt.plot(Temperature, total_ARR, label='FUDGE')


# constants for Au197 +n ->Au198
#a_zero= 19.55230
#a_one= 0.000000
#a_two= 0.000000
#a_three= -2.702790
#a_four= 1.959220
#a_five=  -0.6369410
#a_six= 0.000
# constants for Fe56+ n ->Fe57
#a_zero= 13.09190
#a_one= 0.000000
#a_two= 0.000000
#a_three= 1.018590
#a_four= 0000
#a_five=  00000
#a_six= 0.000
#constants for n+p->d
#Kadonis 2002
K_a_zero=3.746377
K_a_one=0.01391141
K_a_two=-2.672281
K_a_three=10.48889
K_a_four=-0.7198877
K_a_five=0.04782074
K_a_six=-3.543090
K_ARR_reaclib_pre=[]
K_ARR_reaclib=[]
for iT in Temperature:
    K_inside=(K_a_zero+ (K_a_one*pow(iT,-1))+ (K_a_two*pow(iT,-.33333333333))+ (K_a_three*pow(iT,.3333333333)) + (K_a_four*pow(iT,1))+ (K_a_five*pow(iT,1.66666666666)) + (K_a_six*np.log(iT)))
    K_ARR_reaclib_section=np.exp(K_inside)
    K_ARR_reaclib_pre.append(K_ARR_reaclib_section)
for iR in K_ARR_reaclib_pre:
    number=iR*pow(10,-6)
    K_ARR_reaclib.append(number)
reaclibplot=plt.plot(Temperature, K_ARR_reaclib, label='Kadonis')
#Bruslib 2006 first one
B1_a_zero=10.75480
B1_a_one=0
B1_a_two=0
B1_a_three=-2.304720
B1_a_four=-0.8878620
B1_a_five=0.1376630
B1_a_six=0
B1_ARR_reaclib_pre=[]
B1_ARR_reaclib=[]
for iT in Temperature:
    B1_inside=(B1_a_zero+ (B1_a_one*pow(iT,-1))+ (B1_a_two*pow(iT,-.33333333333))+ (B1_a_three*pow(iT,.3333333333)) + (B1_a_four*pow(iT,1))+ (B1_a_five*pow(iT,1.66666666666)) + (B1_a_six*np.log(iT)))
    B1_ARR_reaclib_section=np.exp(B1_inside)
    B1_ARR_reaclib_pre.append(B1_ARR_reaclib_section)
for iR in B1_ARR_reaclib_pre:
    number=iR*pow(10,-6)
    B1_ARR_reaclib.append(number)
reaclibplot=plt.plot(Temperature, B1_ARR_reaclib, label='BRUSLIB 1')
#Bruslib 2006 second one
B2_a_zero=8.846880
B2_a_one=0
B2_a_two=0
B2_a_three=-0.01020820
B2_a_four=-0.08939590
B2_a_five=0.006967040
B2_a_six=1.0
B2_ARR_reaclib_pre=[]
B2_ARR_reaclib=[]
for iT in Temperature:
    B2_inside=(B2_a_zero+ (B2_a_one*pow(iT,-1))+ (B2_a_two*pow(iT,-.33333333333))+ (B2_a_three*pow(iT,.3333333333)) + (B2_a_four*pow(iT,1))+ (B2_a_five*pow(iT,1.66666666666)) + (B2_a_six*np.log(iT)))
    B2_ARR_reaclib_section=np.exp(B2_inside)
    B2_ARR_reaclib_pre.append(B2_ARR_reaclib_section)
for iR in B2_ARR_reaclib_pre:
    number=iR*pow(10,-6)
    B2_ARR_reaclib.append(number)
#reaclibplot=plt.plot(Temperature, B2_ARR_reaclib, label='BRUSLIB 2')
#Bruslib 2006 third one
B3_a_zero=12.36870
B3_a_one=0
B3_a_two=0
B3_a_three=-2.706180
B3_a_four=11.71800
B3_a_five=-.003127880
B3_a_six=0.4691270
ARR_reaclib_pre=[]
ARR_reaclib=[]
for iT in Temperature:
    inside=(B3_a_zero+ (B3_a_one*pow(iT,-1))+ (B3_a_two*pow(iT,-.33333333333))+ (B3_a_three*pow(iT,.3333333333)) + (B3_a_four*pow(iT,1))+ (B3_a_five*pow(iT,1.66666666666)) + (B3_a_six*np.log(iT)))
    ARR_reaclib_section=np.exp(inside)
    ARR_reaclib_pre.append(ARR_reaclib_section)
for iR in ARR_reaclib_pre:
    number=iR*pow(10,-6)
    ARR_reaclib.append(number)
#reaclibplot=plt.plot(Temperature, ARR_reaclib, label='BRUSLIB 3')
plt.xlabel('$Temperature (GK)$')
plt.ylabel('$ Astrophysical Reaction Rate(1e6)(cm^{3} s^{-1} mol^{-1})$')
plt.title('H neutron capture Astrophysical Reaction Rate')
plt.legend(loc='upper center', shadow=True, fontsize='x-large')
plt.show()





