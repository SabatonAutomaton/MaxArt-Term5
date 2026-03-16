#include "MSGFR.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>

void MSGFR::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { 0.5, 0.5 };
   sk = ops.MultVecScal( f.gradFunc( x0 ), -1 );
   sk_1 = new double[n];
   xk = new double[n];
   xk_1 = new double[n];
}

double *MSGFR::Solver( )
{
   const double PI = 3.14159265359;
   std::ofstream file1( "out1.txt" );
   std::ofstream file2("out2.txt" );
   double a = 0, b = 10;
   lambdak = Min( eps1, a, b );
   xk = ops.AddVec( x0, ops.MultVecScal( sk, lambdak ) );
   xk_1 = x0;
   int k = 0;
   while ( ops.DotProduct( sk ) >= eps)
   {
      if ( k % n1 == 0 )
      {
         sk = ops.MultVecScal( f.gradFunc( xk ), -1 );
      }
      lambdak = Min( eps1, a, b );
      xk_1 = xk;
      xk = ops.AddVec( xk_1, ops.MultVecScal( sk, lambdak ) );
      sk_1 = sk;
      double wk = CalcW( );
      sk = ops.AddVec( ops.MultVecScal( f.gradFunc( xk ), -1 ), ops.MultVecScal( sk_1, wk ) );
      k++;
      std::cout << k << "\t";
      file1 << k << " ";
      file2 << k << " ";
      for ( int i = 0; i < n; i++ )
      {
         std::cout << std::setprecision( 16 ) << xk[i] << ' ';
         file1 << std::setprecision( 16 ) << xk[i] << ' ';
         file2 << std::setprecision( 16 ) << xk[i] << ' ';
      }
      file1 << std::setprecision( 16 ) << f.func( xk ) << ' ' << sk[0]<< ' ' << sk[1] << ' ' << lambdak<< ' ';
      file1 << abs( xk[0] - xk_1[0] )<< ' ' << abs( xk[1] - xk_1[1] )<< ' ' << abs( f.func( xk ) - f.func(xk_1))<< ' ';
      file1 << std::acos(( xk[0] * sk[0] + xk[1] * sk[1] ) / ( sqrt( ops.DotProduct( xk ) ) * sqrt( ops.DotProduct( sk ) ) ))*180/PI;
      std::cout << std::endl;
      file1 << std::endl;
      file2 << std::endl;
   }
   return xk;
}

double MSGFR::CalcW( )
{
   return ops.DotProduct( f.gradFunc(xk) ) / ops.DotProduct( f.gradFunc( xk_1 ) );
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


