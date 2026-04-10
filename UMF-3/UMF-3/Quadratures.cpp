#include "Quadratures.h"
#include <cmath>

std::vector<std::vector<double>> Quadratures::gauss3( )
{
   std::vector<std::vector<double>> res( 3 );
   double t1 = -std::sqrt( 15.0 ) / 5.0;
   double t2 = 0.0;
   double t3 = +std::sqrt( 15.0 ) / 5.0;
   double w1 = 5.0 / 9.0;
   double w2 = 8.0 / 9.0;
   double w3 = 5.0 / 9.0;

   res[0] = std::vector<double>{ ( t1 + 1 ) / 2, 0.5 * w1 };
   res[1] = std::vector<double>{ ( t2 + 1 ) / 2, 0.5 * w2 };
   res[2] = std::vector<double>{ ( t3 + 1 ) / 2, 0.5 * w3 };
   return res;
}

std::vector<std::vector<double>> Quadratures::gaussHex2( )
{
   const double a = 1.0 / std::sqrt( 3.0 );
   const double p0 = 0.5 * ( 1.0 - a );
   const double p1 = 0.5 * ( 1.0 + a );
   const double w = 0.125; // (1/2)^3 для перехода [-1,1]^3 -> [0,1]^3

   std::vector<std::vector<double>> q;
   q.reserve( 8 );
   for ( int ir = 0; ir < 2; ++ir )
   {
      for ( int is = 0; is < 2; ++is )
      {
         for ( int it = 0; it < 2; ++it )
         {
            q.push_back( { ir ? p1 : p0, is ? p1 : p0, it ? p1 : p0, w } );
         }
      }
   }
   return q;
}
