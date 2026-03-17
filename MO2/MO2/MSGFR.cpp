#include "MSGFR.h"
#include "iostream"
#include <iomanip>
#include <fstream>
#include <cmath>

void MSGFR::Init( int n1 )
{
   n = n1;
   x0 = new double[2] { 5, 6};
   sk = ops.MultVecScal(f.gradFunc( x0 ), -1 );
   sk_1 = new double[n];
   xk = new double[n];
   xk_1 = new double[n];
}

void MSGFR::FindInterval(double& a, double& b)
{
    double h = 0.1;
    double lambda0 = 0.5;

    double f0 = f.funcInDirection(x0, lambda0, sk);
    double f1 = f.funcInDirection(x0, lambda0 + h, sk);

    if (f1 > f0)
        h = -h;

    double lambda_prev = lambda0;
    double lambda_curr = lambda0 + h;

    int iter = 0;

    while (iter < 1000)   // защита от бесконечного цикла
    {
        double f_prev = f.funcInDirection(x0, lambda_prev, sk);
        double f_curr = f.funcInDirection(x0, lambda_curr, sk);

        if (f_curr > f_prev)
        {
            a = std::min(lambda_prev, lambda_curr);
            b = std::max(lambda_prev, lambda_curr);
            return;
        }

        lambda_prev = lambda_curr;
        lambda_curr += h;
        h *= 2;
        iter++;
    }

    a = 0;
    b = 1;
}

double *MSGFR::Solver( )
{
   const double PI = 3.14159265359;
   std::ofstream file1( "out1.txt" );
   file1 << "i\t" << "x\t" << "y\t" << "f\t" << "s1\t" << "s2\t" << "lambda\t" << "|xi - xi-1|\t" << "|yi - yi-1|\t" << "|fi - fi-1|\t" << "angle\t" << "grad1\tgrad2" << std::endl;
   std::ofstream file2( "out2.txt" );
   file2 << 0 << ' ' << x0[0] << ' ' << x0[1] << std::endl;
   double a, b;
   FindInterval(a, b);
   lambdak = GoldenRatio(eps1, a, b);
   xk = ops.AddVec( x0, ops.MultVecScal( sk, lambdak ) );
   xk_1 = x0;
   int k = 0;
   while ( ops.DotProduct( sk, sk ) >= eps )
   {
      if ( k % n1 == 0 )
      {
         sk = ops.MultVecScal( f.gradFunc( xk ), -1 );
      }
      FindInterval(a, b);
      lambdak = GoldenRatio( eps1, a, b );
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


double MSGFR::GoldenRatio( double eps, double an, double bn )
{
   double x, x1, x2;
   int n = 1;
   double const C = 0.38196601125010515;
   x1 = an + C * ( bn - an );
   x2 = bn - C * ( bn - an );
   double f1 = f.funcInDirection( xk, x1, sk );
   double f2 = f.funcInDirection( xk, x2, sk );
   while ( bn - an >= eps )
   {
      x = 0.5 * ( an + bn );
      n++;
      if ( f1 < f2 )
      {
         bn = x2;
         x2 = x1;
         f2 = f1;
         x1 = an + C * ( bn - an );
         f1 = f.funcInDirection( xk, x1, sk );
      }
      else if ( f1 > f2 )
      {
         an = x1;
         x1 = x2;
         f1 = f2;
         x2 = bn - C * ( bn - an );
         f2 = f.funcInDirection( xk, x2, sk );
      }
      else
      {
         break;
      }
   }
   return x;
}


