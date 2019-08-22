/*
 * # <<BEGIN-copyright>>
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
