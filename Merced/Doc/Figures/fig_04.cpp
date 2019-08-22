/*
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// code to calculate hyperbolas E = const in the (E', \eta) plane
// Fig. 4 in the xndfgen document
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

  int num_Ein = 10;
  double Ein_min = 7.0;
  double Ein_max = 8.0;
  double dEin = ( Ein_max - Ein_min )/num_Ein;
  int num_Eout = 2;
  double Eout_min = 2.0;
  double Eout_max = 3.0;
  double dEout = ( Eout_max - Eout_min )/(num_Eout-1);

  for( int i_Eout = 0; i_Eout < num_Eout; ++i_Eout )
  {
    double Eout = Eout_min + i_Eout*dEout;
    cout << "% Eout: " << Eout << endl;
    for( int i_Ein = 0; i_Ein <= num_Ein; ++i_Ein )
    {
      double Ein = Ein_min + i_Ein * dEin;
      double eta = ( Eout - ( alpha + gamma )*Ein - beta * Q )/
	(2.0 * sqrt( gamma * Ein * (alpha * Ein + beta * Q)));
      cout << "(" << Ein << ", " << eta << ")%" << endl;
    }
  }
}
