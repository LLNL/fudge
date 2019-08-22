/*
 * # <<BEGIN-copyright>>
  Copyright (c) 2017, Lawrence Livermore National Security, LLC.
  Produced at the Lawrence Livermore National Laboratory.
  Written by the LLNL Nuclear Data and Theory group
          (email: mattoon1@llnl.gov)
  LLNL-CODE-725546.
  All rights reserved.
  
  This file is part of the Merced package, used to generate nuclear reaction
  transfer matrices for deterministic radiation transport.
  
  
      Please also read this link - Our Notice and Modified BSD License
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the disclaimer below.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the disclaimer (as noted below) in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of LLNS/LLNL nor the names of its contributors may be used
        to endorse or promote products derived from this software without specific
        prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
  THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  
  Additional BSD Notice
  
  1. This notice is required to be provided under our contract with the U.S.
  Department of Energy (DOE). This work was produced at Lawrence Livermore
  National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
  
  2. Neither the United States Government nor Lawrence Livermore National Security,
  LLC nor any of their employees, makes any warranty, express or implied, or assumes
  any liability or responsibility for the accuracy, completeness, or usefulness of any
  information, apparatus, product, or process disclosed, or represents that its use
  would not infringe privately-owned rights.
  
  3. Also, reference herein to any specific commercial products, process, or services
  by trade name, trademark, manufacturer or otherwise does not necessarily constitute
  or imply its endorsement, recommendation, or favoring by the United States Government
  or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
  herein do not necessarily state or reflect those of the United States Government or
  Lawrence Livermore National Security, LLC, and shall not be used for advertising or
  product endorsement purposes.
  
 * # <<END-copyright>>
 */
// code to calculate hypoerbolas \eta = const for center-of-mass coordinates
// Fig. 3 in the xndfgen document
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  double myi = 1.0;  // incident particle
  double myo = 1.0;  // emitted particle
  double mtarg = 6.0;  // target
  double mres = mtarg + ( myi - myo );  // approximate residual
  double gamma = myi * myo /( ( myi + mtarg ) * ( myi + mtarg ) );
  double beta = mres / ( myo + mres );
  double alpha = beta * mtarg / ( myi + mtarg );
  double Q = -4.0;
  double threshold = -beta * Q / alpha;

  int num_eta = 2;
  int num_E = 10;
  double range = 5.0;
  double dE = range/( num_E * num_E );

  for( int i_eta = -num_eta; i_eta <= num_eta; ++i_eta )
  {
    double eta = i_eta * 1.0 / num_eta;
    cout << "% eta: " << eta << endl;
    cout << "\\spline(" << threshold << ", " << gamma * threshold <<
      ")%" << endl;
    for( int i_E = 1; i_E <= num_E; ++i_E )
    {
      // grid to resolve the sqrt at the threshold
      double E = threshold + i_E * i_E * dE;
      double E_lab = ( alpha + gamma )*E + beta * Q +
        2.0 * eta * sqrt( gamma * E * (alpha * E + beta * Q));
      cout << "(" << E << ", " << E_lab << ")%" << endl;
    }
  }
}
