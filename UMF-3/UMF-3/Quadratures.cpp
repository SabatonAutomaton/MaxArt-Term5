#include "Quadratures.h"
#include <iostream>

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
