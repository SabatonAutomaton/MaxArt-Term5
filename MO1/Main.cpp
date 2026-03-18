#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

double func( double x )
{
   return x*x;
}

double dichotomy( double eps, double an, double bn )
{
   std::ofstream file( "dichotomy.txt" );
   double x, x1, x2;
   double delta = 0.5 * eps;
   int n = 0;
   file << an << '\n' << bn << '\n';
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
         bn = x2;
      }
      else if ( f1 > f2 )
      {
         an = x1;
      }
      else
      {
         break;
      }
      std::cout << n << '\t' << std::setprecision( 16 ) << x << std::endl;
      file << std::setprecision( 16 ) << x << std::endl;
   }
   file.close( );
   return x;
}

double goldenRatio( double eps, double an, double bn )
{
   std::ofstream file( "goldenRatio.txt" );
   double x, x1, x2;
   int n = 1;
   double const C = 0.38196601125010515;
   x1 = an + C * ( bn - an );
   x2 = bn - C * ( bn - an );
   double f1 = func( x1 );
   double f2 = func( x2 );
   file << an << '\n' << bn << '\n';
   while ( bn - an >= eps )
   {
      x = 0.5 * ( an + bn );
      std::cout << n << '\t' << std::setprecision( 16 ) << x << std::endl;
      file << std::setprecision( 16 ) << x << std::endl;
      n++;
      if ( f1 < f2 )
      {
         bn = x2;
         x2 = x1;
         f2 = f1;
         x1 = an + C * ( bn - an );
         f1 = func( x1 );
      }
      else if ( f1 > f2 )
      {
         an = x1;
         x1 = x2;
         f1 = f2;
         x2 = bn - C * ( bn - an );
         f2 = func( x2 );
      }
      else
      {
         break;
      }
   }
   file.close( );
   return x;
}

int FindN( double eps, double a0, double b0, std::vector<double> &res )
{
   int Fn = 1;      // F(1)
   int Fn1 = 1;     // F(2) 
   int Fn2 = 2;     // F(3)
   int n = 1;
   while ( b0 - a0 > Fn2 * eps )
   {
      int temp = Fn2;
      Fn2 = Fn2 + Fn1;
      Fn = Fn1;
      Fn1 = temp;
      n++;
      res.reserve( res.size( ) + 1 );
      res.push_back( Fn2 );
   }
   return n;
}

double Fibonacci( double eps, double a0, double b0 )
{
   std::ofstream file( "fibonacci.txt" );
   double x, x1, x2;
   double an = a0, bn = b0;
   int k = 2;
   std::vector<double> fibonacciFunc = { 1,1,2 };
   double f1, f2;
   double buff;
   int n = FindN( eps, a0, b0, fibonacciFunc );
   x1 = a0 + fibonacciFunc[n - 1] / fibonacciFunc[n + 1] * ( b0 - a0 );
   x2 = a0 + fibonacciFunc[n] / fibonacciFunc[n + 1] * ( b0 - a0 );
   f1 = func( x1 );
   f2 = func( x2 );
   file << a0 << '\n' << b0 << '\n';
   for ( ; k <= n; k++ )
   {

      if ( f1 < f2 )
      {
         x = x1;
         bn = x2;
         x2 = x;
         f2 = f1;
         x1 = an + fibonacciFunc[n - k] / fibonacciFunc[n + 1] * ( b0 - a0 );
         f1 = func( x1 );
      }
      else if ( f1 > f2 )
      {
         x = x2;
         an = x1;
         x1 = x;
         f1 = f2;
         x2 = an + fibonacciFunc[n - k + 1] / fibonacciFunc[n + 1] * ( b0 - a0 );
         f2 = func( x2 );
      }
      else
      {
         break;
      }
      std::cout << k - 1 << '\t' << std::setprecision( 16 ) << x << std::endl;
      file << std::setprecision( 16 ) << x << std::endl;
   }
   file.close( );
   return x;
}

double parabola( double eps, double an, double bn )
{
   double x1 = an, x2 = ( an + bn ) / 2.0, x3 = bn;
   double x = x2 + 10 *eps;
   double f, f1 = func( x1 ), f2 = func( x2 ), f3 = func( x3 );
   double numer, denom; // числитель, знаменатель
   int k = 0;
   while ( fabs( x2 - x ) >= eps )
   {
      numer = ( x2 - x1 ) * ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( x2 - x3 ) * ( f2 - f1 );
      denom = ( x2 - x1 ) * ( f2 - f3 ) - ( x2 - x3 ) * ( f2 - f1 );
      if ( denom == 0 )
         x = x2;
      else
         x = x2 - 0.5 * numer / denom;
      f = func( x );
      if ( x <= x2 )
         if ( f <= f2 )
         {
            x3 = x2;
            f3 = f2;
            x2 = x;
            f2 = f;
         }
         else
         {
            x1 = x;
            f1 = f;
         }
      else
         if ( f <= f2 )
         {
            x1 = x2;
            x2 = x;
            f1 = f2;
            f2 = f;
         }
         else
         {
            x3 = x;
            f3 = f;
         }
      k++;
      std::cout << k << '\t' << std::setprecision( 16 ) << x << std::endl;
   }
   return x2;
}

int main( )
{
   double eps, an, bn;
   std::ifstream file( "in.txt" );
   file >> eps;
   file >> an;
   file >> bn;
   std::cout << "Dichotomy" << std::endl;
   dichotomy( eps, an, bn );
   std::cout << "Golden" << std::endl;
   goldenRatio( eps, an, bn );
   std::cout << "Fibonacci" << std::endl;
   Fibonacci( eps, an, bn );
   std::cout << "parabola" << std::endl;
   parabola( eps, an, bn );
   return 0;
}