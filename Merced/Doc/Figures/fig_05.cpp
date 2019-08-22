/*
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// code to graph curves for the quadrature box in Fig. 5
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  double Ein0 = 1.0;
  double Ein1 = 3.0;
  double dE0 = 1.0;  // range of outgoing energies at Ein0
  double dE1 = 4.0;  // range of outgoing energies at Ein4
  int num_Ein = 10;
  double Ein_min = 0.8;
  double Ein_max = 2.4;
  double dEin = ( Ein_max - Ein_min )/num_Ein;
  int num_Eout = 2;
  double Eout_min = 0.8;
  double Eout_max = 1.8;
  double dEout = ( Eout_max - Eout_min )/(num_Eout-1);

  for( int i_Eout = 0; i_Eout < num_Eout; ++i_Eout )
  {
    double Eout = Eout_min + i_Eout*dEout;
    cout << "% Eout: " << Eout << endl;
    for( int i_Ein = 0; i_Ein <= num_Ein; ++i_Ein )
    {
      double Ein = Ein_min + i_Ein * dEin;
      double alpha = (Ein - Ein0)/(Ein1 - Ein0 );
      double dE = (1.0 - alpha)*dE0 + alpha*dE1;  // interpolated outgoing energy range
      double Ehat = Eout/dE;
      // scale by 3 for the plot and shift by 6
      cout << "(" << 6 + Ein << ", " << 3*Ehat << ")%" << endl;
    }
  }
}
