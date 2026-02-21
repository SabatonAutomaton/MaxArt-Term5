#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

double func( double x )
{
   return ( x - 3 ) * ( x - 3 );
}

double dichotomy( double eps, double an, double bn )
{
   double x, x1, x2;
   double delta = 0.5 * eps;
   int n = 0;
   while ( bn - an >= eps )
   {
      n++;
      x = 0.5 * ( an + bn );
      x1 = x - 0.5 * delta;
      x2 = x + 0.5 * delta;
      double f1 = func( x1 );
      double f2 = func( x2 );
      if ( f1 < f2 )
      {
         bn = x1;
      }
      else if ( f1 > f2 )
      {
         an = x2;
      }
      else
      {
         break;
      }
      std::cout << n << '\t'<< std::setprecision(16) << x << std::endl;
   }
   return x;
}

int main( )
{
   double eps, an, bn;
   std::ifstream file( "in.txt" );
   file >> eps;
   file >> an;
   file >> bn;
   dichotomy( eps, an, bn );
   return 0;
}