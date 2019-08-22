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

from fudge.legacy.endl.structure import masses as massModule
import collections

class massTracker:
    """
    ENDF has complex rules for particle masses, so the massTracker serves as a single interface for masses.
    It stores ground state masses in amu, and has methods to read in from ENDF AWR
    (converting from nuclear to atomic where appropriate).
    It can return either atomic mass in amu or as AWR (converting back to nuclear where appropriate)

    When reading in from AWR it also compares computed atomic mass to AME2003, and warns if values differ
    by more than electronMass / 10
    """
    eVc2 = 1.07354416620656e-9   # in amu
    electronMass = 0.0005485801
    neutronMass = massModule.getMassFromZA(1)
    electronBindingEnergiesAmu = {
        1001: -13.8885758218 * eVc2,
        1002: -13.9351505263 * eVc2,
        1003: -14.3729527487 * eVc2,
        2003: -70.5606714777 * eVc2,
        2004: -75.2181415144 * eVc2,
        2005: -78.    * eVc2,
        2006: -79.94  * eVc2,
        2007: -79.99  * eVc2
    }
    elementalMass = {
         1000: 9.992414e-01, 2000: 3.968215e+00, 3000: 6.881371e+00, 4000: 8.934758e+00, 5000: 1.071713e+01,
         6000: 1.190782e+01, 7000: 1.388637e+01, 8000: 1.586195e+01, 9000: 1.883519e+01, 10000: 2.000565e+01,
         11000: 2.279230e+01, 12000: 2.409620e+01, 13000: 2.674971e+01, 14000: 2.784422e+01, 15000: 3.070771e+01,
         16000: 3.178458e+01, 17000: 3.514843e+01, 18000: 3.960482e+01, 19000: 3.876242e+01, 20000: 3.973568e+01,
         21000: 4.456969e+01, 22000: 4.748850e+01, 23000: 5.050387e+01, 24000: 5.154932e+01, 25000: 5.446604e+01,
         26000: 5.536723e+01, 27000: 5.842692e+01, 28000: 5.819572e+01, 29000: 6.300009e+01, 30000: 6.481834e+01,
         31000: 6.912106e+01, 32000: 7.196640e+01, 33000: 7.427797e+01, 34000: 7.828168e+01, 35000: 7.921757e+01,
         36000: 8.308009e+01, 37000: 8.473357e+01, 38000: 8.686728e+01, 39000: 8.814213e+01, 40000: 9.043635e+01,
         41000: 9.210826e+01, 42000: 9.511580e+01, 43000: 9.715810e+01, 44000: 1.002017e+02, 45000: 1.020210e+02,
         46000: 1.054859e+02, 47000: 1.069413e+02, 48000: 1.114443e+02, 49000: 1.138336e+02, 50000: 1.176704e+02,
         51000: 1.207041e+02, 52000: 1.265038e+02, 53000: 1.258138e+02, 54000: 1.301720e+02, 55000: 1.317632e+02,
         56000: 1.361502e+02, 57000: 1.377117e+02, 58000: 1.389163e+02, 59000: 1.396975e+02, 60000: 1.430009e+02,
         61000: 1.437543e+02, 62000: 1.491080e+02, 63000: 1.506545e+02, 64000: 1.558991e+02, 65000: 1.575597e+02,
         66000: 1.611040e+02, 67000: 1.635131e+02, 68000: 1.658231e+02, 69000: 1.674827e+02, 70000: 1.715535e+02,
         71000: 1.734639e+02, 72000: 1.769566e+02, 73000: 1.793935e+02, 74000: 1.822706e+02, 75000: 1.846073e+02,
         76000: 1.885660e+02, 77000: 1.905687e+02, 78000: 1.934140e+02, 79000: 1.952739e+02, 80000: 1.988668e+02,
         81000: 2.026143e+02, 82000: 2.054200e+02, 83000: 2.071847e+02, 84000: 2.072045e+02, 85000: 2.081959e+02,
         86000: 2.200928e+02, 87000: 2.210843e+02, 88000: 2.240585e+02, 89000: 2.250499e+02, 90000: 2.300446e+02,
         91000: 2.290155e+02, 92000: 2.359841e+02, 93000: 2.349640e+02, 94000: 2.419039e+02, 95000: 2.409124e+02,
         96000: 2.448781e+02, 97000: 2.448781e+02, 98000: 2.488437e+02, 99000: 2.498351e+02, 100000: 2.547922e+02}

    def __init__(self):
        self.amuMasses = { 1: self.neutronMass }
        self.ZA_AWRMasses = {}
        self.ZA_AWRMasses_nuclear = {}

    def addMassAMU(self, ZA, massInAmu):
        self.amuMasses[ ZA ] = massInAmu

    def addMassAWR(self, ZA, AWR, asTarget=True):

        # what's this all about? Found in ITYPE_0_Misc readMF6
        """
        if( ZAP not in info.ZAMasses ) :
            info.ZAMasses[ZAP] = AWP * info.ZAMasses[1]
        elif( info.ZAMasses[ZAP] is None ) :
            info.ZAMasses[ZAP] = -AWP * info.ZAMasses[1]
        """
        if ZA==0: return    # exclude gammas

        warning = ''
        amuMass = AWR * self.neutronMass
        Z, A = divmod(ZA, 1000)
        if not asTarget and Z in (1,2): # AWR should be a nuclear mass rather than atomic
            self.ZA_AWRMasses_nuclear.setdefault( ZA, collections.Counter() ).update( [AWR] )
            amuMass += self.electronMass * Z + self.electronBindingEnergiesAmu[ZA]
        else:
            self.ZA_AWRMasses.setdefault( ZA, collections.Counter() ).update( [AWR] )

        AMEmass = massModule.getMassFromZA(ZA)
        if AMEmass is not None:
            ratioToAME = amuMass / AMEmass
            #if abs(AMEmass - amuMass) > self.electronMass / 10:
            #    warning = " WARNING: ENDF differs from AME by more than electron mass/10: %g vs %g" % (amuMass, AMEmass)
            if ratioToAME < 0.9 or ratioToAME > 1.1:
                warning = " WARNING: ENDF/AME masses differ: %g vs %g" % (amuMass, AMEmass)

        self.amuMasses[ZA] = amuMass
        return warning


    def getMassAMU(self, ZA):
        """
        First checks internal table. If requested ZA not present there, query the AME database using massModule.
        If still not found, attempt to compute it from nearest neighbors.
        """

        mass = self.amuMasses.get(ZA)
        if mass is not None:
            if mass < 0:
                raise Exception("Negative mass encountered?")
                mass *= -1
                self.addMassAMU( ZA, mass )
        else :
            mass = massModule.getMassFromZA( ZA )
            if( mass is None ) :               # If mass is not found, let's attempt to make something up that is close.
                A = ZA % 1000
                for idx in xrange( 1, A ) :
                    massLower = massModule.getMassFromZA( ZA - idx )
                    massHigher = massModule.getMassFromZA( ZA + idx )
                    if( massLower is not None ) :
                        mass = massLower + idx * massModule.getMassFromZA( 1 )
                        break
                    elif( massHigher is not None ) :
                        mass = massHigher - idx * massModule.getMassFromZA( 1 )
                        break
                if( mass is None ) : raise Exception( 'Could not obtain mass for ZA = %s' % ZA )
                mass = round(mass,1)    # since we're extrapolating
        return mass

    def getMassAWR(self, ZA, asTarget=True, levelEnergyInEv=0):
        """For going back to ENDF. Convert light charged particles back to nuclear masses (except for targets)."""

        try :           # BRBFIXME - fails for ENDF-B-VII.1/photoat/photoat-001_H_000.endf because ZA = 1000 not in self.amuMasses.
            massAMU = self.amuMasses[ZA]
        except :
            if( ( ZA % 1000 ) != 0 ) : raise
            massAMU = massModule.getMassFromZA( ZA )

        if not asTarget:
            Z,A = divmod(ZA,1000)
            if Z in (1,2):
                massAMU -= self.electronMass * Z - self.electronBindingEnergiesAmu[ZA]
        if levelEnergyInEv:
            massAMU += levelEnergyInEv * self.eVc2
        return massAMU / self.neutronMass

    def getElementalMassAMU(self, ZA):
        return self.elementalMass[ZA]

    def getElementalMassAWR(self, ZA):
        return self.getElementalMassAMU(ZA) / self.neutronMass

    def useMostCommonAMUmasses(self):
        """
        ENDF evaluations may have internal discrepancies, including for particle masses.
        self.ZA_AWRMasses keeps track of how many times each value appears. This method updates
        self.amuMasses to reflect the value that appears most often in ZA_AWRMasses.
        """
        for ZA in self.ZA_AWRMasses:
            AWR, count = self.ZA_AWRMasses[ZA].most_common(1)[0]
            self.addMassAMU( ZA, AWR * self.neutronMass )
