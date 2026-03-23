#include "Basis.h"
#include <cmath>

double Basis::basis1( double x )
{
   return 2 * ( x - 0.5 ) * ( x - 1 );
}

double Basis::basis1Grad( double x )
{
   return 4 * x - 3;
}

double Basis::basis2( double x )
{
   return -4 * x * ( x - 1 );
}

double Basis::basis2Grad( double x )
{
   return -4 * ( 2 * x - 1 );
}

double Basis::basis3( double x )
{
   return 2 * x * ( x - 0.5 );
}

double Basis::basis3Grad( double x )
{
   return 4 * x - 1;
}

double **Basis::Quadrature( )
{
   double **res = new double *[3];
   double t1 = -std::sqrt( 15.0 ) / 5.0;
   double t2 = 0.0;
   double t3 = +std::sqrt( 15.0 ) / 5.0;
   double w1 = 5.0 / 9.0;
   double w2 = 8.0 / 9.0;
   double w3 = 5.0 / 9.0;

   res[0] = new double[2] {( t1 + 1 ) / 2, 0.5 * w1};
   res[1] = new double[2] {( t2 + 1 ) / 2, 0.5 * w2};
   res[2] = new double[2] {( t3 + 1 ) / 2, 0.5 * w3};
   return res;
}
