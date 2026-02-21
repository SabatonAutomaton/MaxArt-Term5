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
   }
   return x;
}

double goldenRatio( double eps, double an, double bn )
{
   double x, x1, x2;
   int n = 1;
   double const C = 0.38196601125010515;
   x1 = an + C * ( bn - an );
   x2 = bn - C * ( bn - an );
   double f1 = func( x1 );
   double f2 = func( x2 );
   while ( bn - an >= eps )
   {
      x = 0.5 * ( an + bn );
      std::cout << n << '\t' << std::setprecision( 16 ) << x << std::endl;
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
   return x;
}

//int FindN( double eps, double a0, double b0 )
//{
//   int Fn = 1, Fn1 = 1, Fn2 = 2;
//   int buff;
//   int n = 0;
//   for ( ; b0 - a0 > Fn2 * eps; n++ )
//   {
//      Fn = Fn1;
//      Fn1 = Fn2;
//      Fn2 = Fn + Fn1;
//   }
//   return n;
//}
//double FibonacciFunc( int n )
//{
//   int Fn = 1, Fn1 = 1, Fn2 = 2;
//   int buff;
//   for ( int i = 0; i < n - 3; i++ )
//   {
//      Fn = Fn1;
//      Fn1 = Fn2;
//      Fn2 = Fn + Fn1;
//      
//   }
//   return Fn2;
//}

int FindN( double eps, double a0, double b0 )
{
   int Fn = 1;      // F(1)
   int Fn1 = 1;     // F(2) 
   int Fn2 = 2;     // F(3)
   int n = 1;       // Начинаем с n=1

   while ( b0 - a0 > Fn2 * eps ) {
      int temp = Fn2;
      Fn2 = Fn2 + Fn1;
      Fn = Fn1;
      Fn1 = temp;
      n++;
   }
   return n;
}

double FibonacciFunc( int n )
{
   if ( n <= 0 ) return 0;
   if ( n == 1 || n == 2 ) return 1;

   int Fn = 1;      // F(1)
   int Fn1 = 1;     // F(2)
   int Fn2 = 2;     // F(3)

   for ( int i = 3; i < n; i++ ) {
      int temp = Fn2;
      Fn2 = Fn2 + Fn1;
      Fn = Fn1;
      Fn1 = temp;
   }

   return Fn2;
}
double Fibonacci( double eps, double a0, double b0 )
{
   double x, x1, x2;
   double an = a0, bn = b0;
   int k = 2;
   double f1, f2;
   double buff;
   int n = FindN( eps, a0, b0 );
   x1 = a0 + FibonacciFunc( n ) / FibonacciFunc( n + 2 ) * ( b0 - a0 );
   x2 = a0 + FibonacciFunc( n + 1 ) / FibonacciFunc( n + 2 ) * ( b0 - a0 );
   f1 = func( x1 );
   f2 = func( x2 );
   for ( ; k <= n; k++ )
   {

      if ( f1 < f2 )
      {
         x = x1;
         bn = x2;
         x2 = x;
         f2 = f1;
         x1 = an + FibonacciFunc( n - k + 1 ) / FibonacciFunc( n + 2 ) * ( b0 - a0 );
         f1 = func( x1 );
      }
      else if ( f1 > f2 )
      {
         x = x2;
         an = x1;
         x1 = x;
         f1 = f2;
         x2 = an + FibonacciFunc( n - k + 2 ) / FibonacciFunc( n + 2 ) * ( b0 - a0 );
         f2 = func( x2 );
      }
      else
      {
         break;
      }
      std::cout << k << '\t' << std::setprecision( 16 ) << x /*<< "\t" << x1 << "\t" << x2*/ << std::endl;
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
   //for ( int i = 3; i <= 15; i++ )
   //{
   //   std::cout << i - 1<< "\t" << FibonacciFunc( i ) << std::endl;
   //}
   std::cout << "Dichotomy" << std::endl;
   dichotomy( eps, an, bn );
   std::cout << "Golden" << std::endl;
   goldenRatio(eps,an,bn );
   std::cout << "Fibonacci" << std::endl;
   Fibonacci( eps, an, bn );
   return 0;
}