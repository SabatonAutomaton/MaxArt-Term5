#include "MSGFR.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>

void MSGFR::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { 6, 5};
   sk = ops.MultVecScal(f.gradFunc( x0 ), -1 );
   sk_1 = new double[n];
   xk = new double[n];
   xk_1 = new double[n];
}

double *MSGFR::Solver( )
{
   const double PI = 3.14159265359;
   std::ofstream file1( "out1.txt" );
   file1 << "i\t" << "x\t" << "y\t" << "f\t" << "s1\t" << "s2\t" << "lambda\t" << "|xi - xi-1|\t" << "|yi - yi-1|\t" << "|fi - fi-1|\t" << "angle\t" << "grad1\tgrad2" << std::endl;
   std::ofstream file2( "out2.txt" );
   file2 << 0 << ' ' << x0[0] << ' ' << x0[1] << std::endl;
   double a = 0, b = 10;
   lambdak = Min( eps1, a, b );
   xk = ops.AddVec( x0, ops.MultVecScal( sk, lambdak ) );
   xk_1 = x0;
   int k = 0;
   while ( ops.DotProduct( sk, sk ) >= eps )
   {
      if ( k % n1 == 0 )
      {
         sk = ops.MultVecScal( f.gradFunc( xk ), -1 );
      }
      lambdak = Min( eps1, a, b );
      for ( int i = 0; i < n; i++ ) xk_1[i] = xk[i];
      xk = ops.AddVec( xk_1, ops.MultVecScal( sk, lambdak ) );
      sk_1 = sk;
      double wk = CalcW( );
      sk = ops.AddVec( ops.MultVecScal( f.gradFunc( xk ), -1 ), ops.MultVecScal( sk_1, wk ) );
      k++;
      std::cout << k << "\t";
      file1 << k << "\t";
      file2 << k << " ";
      for ( int i = 0; i < n; i++ )
      {
         std::cout << std::setprecision( 16 ) << xk[i] << ' ';
         file1 << std::setprecision( 16 ) << xk[i] << '\t';
         file2 << std::setprecision( 16 ) << xk[i] << ' ';
      }
      file1 << std::setprecision( 16 ) << f.func( xk ) << '\t' << sk[0] << '\t' << sk[1] << '\t' << lambdak << '\t';
      file1 << abs( xk[0] - xk_1[0] ) << '\t' << abs( xk[1] - xk_1[1] ) << '\t' << abs( f.func( xk ) - f.func( xk_1 ) ) << '\t';
      file1 << std::acos( ( xk[0] * sk[0] + xk[1] * sk[1] ) / ( sqrt( ops.DotProduct( xk, xk ) ) * sqrt( ops.DotProduct( sk, sk ) ) ) ) * 180 / PI << '\t';
      file1 << f.gradFunc( xk )[0] << '\t' << f.gradFunc( xk )[1];
      std::cout << std::endl;
      file1 << std::endl;
      file2 << std::endl;
   }
   return xk;
}

double MSGFR::CalcW( )
{
   return ops.DotProduct( f.gradFunc( xk ), f.gradFunc( xk ) ) / ops.DotProduct( f.gradFunc( xk_1 ), f.gradFunc( xk_1 ) );
}

double MSGFR::Min( double eps, double an, double bn )
{
   double x, x1, x2;
   double delta = 0.5 * eps1;
   int n = 0;
   while ( bn - an >= eps )
   {
      n++;
      x = 0.5 * ( an + bn );
      x1 = x - 0.5 * delta;
      x2 = x + 0.5 * delta;
      double f1 = f.funcInDirection( xk, x1, sk );
      double f2 = f.funcInDirection( xk, x2, sk );
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
   }
   return x;
}


